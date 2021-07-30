using PowerModels
using NetDecOPF
using DualDecomposition
using BundleMethod
using Ipopt
using OSQP
using JuMP
using LinearAlgebra

PowerModels.silence()

const DD = DualDecomposition
const parallel = DD.parallel
const BM = BundleMethod

function main(file, npartitions::Int, log_path)
    parallel.init()

    if parallel.is_root()
        println("Data file: $file")
    end
    data = parse_file(file)

    # Metis partitioning has been moved to NetDecOPF
    partitions = metis_cluster(data, npartitions)
    parallel.partition(length(partitions))

    has_empty_network = false
    for i in eachindex(partitions)
        if length(partitions[i]) == 0
            global has_empty_network = true
            break
        end
    end

    if has_empty_network == true
        if parallel.is_root()
            println("Network decomposition results in at least one empty sub-network.")
        end
        parallel.finalize()
        return
    end

    # Take care of my partitions only..
    my_sub_index = Int[]
    my_partitions = Vector{Vector{Int}}();
    for s in parallel.getpartition()
        push!(my_sub_index, s)
        push!(my_partitions, partitions[s])
    end

    dn_model = decompose(data, my_partitions, W_global_ACRModel, NetDecOPF.build_acopf_with_free_lines, extra_ref_extensions=[NetDecOPF.ref_add_global_bus!])
    # dn_model = decompose(data, my_partitions, W_ACRModel_V,      NetDecOPF.build_acopf_with_free_lines, extra_ref_extensions=[NetDecOPF.ref_add_global_bus!])
    models = Dict{Int, JuMP.Model}(my_sub_index[i] => dn_model.model_list[i].model for i in eachindex(my_sub_index))
    subproblem_dims = Dict{Int,Tuple{Int,Int}}()

    algo = DD.LagrangeDual()
    for s in my_sub_index
        DD.add_block_model!(algo, s, models[s])
        nconst = 0
        for t in JuMP.list_of_constraint_types(models[s])
            nconst += JuMP.num_constraints(models[s], t[1], t[2])
        end
        subproblem_dims[s] = (nconst,JuMP.num_variables(models[s]))
    end

    for i in 1:parallel.nprocs()
        if i - 1 == parallel.myid()
            for (s, dim) in subproblem_dims
                println("Proc $(i-1): sub $s, cons $(dim[1]), vars $(dim[2])")
            end
        end
        parallel.barrier()
    end

    coupling_variables = Vector{DD.CouplingVariableRef}()
    split_vars = get_split_vars(dn_model)
    for (i,s) in enumerate(my_sub_index)
        vars_dict = split_vars[i]
        for (varname, dict) in vars_dict
            for (idx, vref) in dict
                id = (varname, idx)
                push!(coupling_variables, DD.CouplingVariableRef(s, id, vref))
            end
        end
    end
    DD.set_coupling_variables!(algo, coupling_variables)

    suboptimizer = optimizer_with_attributes(
        Ipopt.Optimizer, 
        "print_level" => 0, 
        "warm_start_init_point" => "yes", 
        "tol" => 1e-4,
        # "linear_solver" => "ma27",
    )

    # Lagrange master method
    # LM = initialize_subgradient_method()
    LM = initialize_bundle_method()

    # Warmup solvers
    dummy_model = JuMP.Model()
    @variable(dummy_model, 0 <= x <= 1)
    @objective(dummy_model, Min, (x-0.5)^2)
    set_optimizer(dummy_model, suboptimizer)
    optimize!(dummy_model)

    for pm in dn_model.model_list
        m = pm.model
        # NOTE: This does not work with ipopt+mumps.
        # for v in JuMP.all_variables(m)
        #     if JuMP.has_lower_bound(v) && JuMP.has_upper_bound(v)
        #         JuMP.set_start_value(v, (JuMP.lower_bound(v) + JuMP.upper_bound(v)) / 2)
        #     end
        # end
        set_optimizer(m, suboptimizer)
    end

    parallel.barrier()
    DD.run!(algo, LM)

    # Write timing outputs to files
    if parallel.is_root()
        if !ispath(log_path)
            mkpath(log_path)
        end
    end
    DD.write_all(algo, dir = log_path)

    parallel.finalize()
end

function initialize_subgradient_method() 
    LM = DD.SubgradientMaster()
    LM.maxiter = 300
    LM.α = 0.0001
    # LM.step_size = (method) -> (LM.α / norm(method.∇f)) # constant step length
    LM.step_size = (method) -> (LM.α / method.iter) # square summable
    # LM.step_size = (method) -> (abs(1.1423e+06 - method.f) / norm(method.∇f)^2) # Polyak's
    return LM
end

function initialize_bundle_method()
    # optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_QPMETHOD" => 2)
    # optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_QPMETHOD" => 4)
    # optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "warm_start_init_point" => "yes")
    optimizer = optimizer_with_attributes(OSQP.Optimizer, "verbose" => false, "linsys_solver" => "mkl pardiso")
    # optimizer = optimizer_with_attributes(OSQP.MathOptInterfaceOSQP.Optimizer, "verbose" => false)

    # Warmup solvers
    dummy_model = JuMP.Model()
    @variable(dummy_model, 0 <= x <= 1)
    @objective(dummy_model, Min, (x-0.5)^2)
    set_optimizer(dummy_model, optimizer)
    optimize!(dummy_model)

    # Change parameters
    params = BM.Parameters()
    BM.set_parameter(params, "maxiter", 3000)
    BM.set_parameter(params, "ϵ_s", 1.e-4)
    BM.set_parameter(params, "ncuts_per_iter", npartitions)

    # Lagrange master method
    return DD.BundleMaster(BM.ProximalMethod, optimizer, params)
    # return DD.BundleMaster(BM.TrustRegionMethod, optimizer, params)
end

file = ARGS[1]
npartitions = parse(Int, ARGS[2])
log_path = ""
try 
    global log_path = ARGS[3]
catch BoundsError
    global log_path = pwd()
end

# file = "/home/kimk/REPOS/pglib-opf/pglib_opf_case5_pjm.m"
# npartitions = 2
# log_path = pwd()

main(file, npartitions, log_path)
