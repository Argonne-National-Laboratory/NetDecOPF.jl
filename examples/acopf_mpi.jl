using PowerModels
using NetDecOPF
using DualDecomposition
using BundleMethod
using Ipopt
using CPLEX
using OSQP
using JuMP
using LinearAlgebra
using SCS
using COSMO

PowerModels.silence()

const DD = DualDecomposition
const parallel = DD.parallel
const BM = BundleMethod

function main(file, npartitions::Int, log_path; max_iter = 3000)
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

    stime = time()

    # dn_model = decompose(data, my_partitions, W_global_ACRModel, NetDecOPF.build_acopf_with_free_lines, extra_ref_extensions=[NetDecOPF.ref_add_global_bus!])
    # dn_model = decompose(data, my_partitions, W_ACRModel_V, NetDecOPF.build_acopf_with_free_lines, extra_ref_extensions=[NetDecOPF.ref_add_global_bus!])
    dn_model = decompose(data, my_partitions, SparseSDPWRMPowerModel, NetDecOPF.build_acopf_with_free_lines)

    models = Dict{Int, JuMP.Model}(my_sub_index[i] => dn_model.model_list[i].model for i in eachindex(my_sub_index))
    subproblem_dims = Dict{Int,Tuple{Int,Int}}()

    if parallel.is_root()
        println("decompose(): $(time()-stime) sec.")
    end

    stime = time()

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

    if parallel.is_root()
        println("Initialize Lagrangian dual: $(time()-stime) sec.")
    end

    # suboptimizer = optimizer_with_attributes(
    #     Ipopt.Optimizer, 
    #     "print_level" => 0, 
    #     "warm_start_init_point" => "yes", 
    #     "tol" => 1e-4,
    #     "linear_solver" => "ma27",
    # )
    # suboptimizer = optimizer_with_attributes(
    #     SCS.Optimizer, 
    #     "verbose" => 0, 
    #     "eps" => 1e-3, 
    #     "alpha" => 1.8, 
    #     "warm_start" => true, 
    #     "linear_solver" => SCS.DirectSolver, 
    #     "acceleration_lookback" => 50
    # )
    suboptimizer = optimizer_with_attributes(
        COSMO.Optimizer, 
        "verbose" => true,
    )

    # Lagrange master method
    # LM = initialize_subgradient_method(;max_iter = max_iter)
    LM = initialize_bundle_method(;max_iter = max_iter)

    for pm in dn_model.model_list
        m = pm.model
        # NOTE: This does not work with ipopt.
        # for v in JuMP.all_variables(m)
        #     if JuMP.has_lower_bound(v) && JuMP.has_upper_bound(v)
        #         JuMP.set_start_value(v, (JuMP.lower_bound(v) + JuMP.upper_bound(v)) / 2)
        #     end
        # end
        set_optimizer(m, suboptimizer)
        # if parallel.is_root()
        #     set_optimizer_attribute(m, "print_level", 5)
        # end
    end

    parallel.barrier()

    # Warm up
    BM.set_parameter(LM.params, "maxiter", 2)
    DD.run!(algo, LM)
    algo.subsolve_time = []
    algo.subcomm_time = []
    algo.subobj_value = []
    algo.master_time = []

    # actual run
    BM.set_parameter(LM.params, "maxiter", max_iter)
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

function initialize_subgradient_method(;max_iter = max_iter) 
    LM = DD.SubgradientMaster()
    LM.maxiter = max_iter
    LM.α = 0.0001
    # LM.step_size = (method) -> (LM.α / norm(method.∇f)) # constant step length
    LM.step_size = (method) -> (LM.α / method.iter) # square summable
    # LM.step_size = (method) -> (abs(1.1423e+06 - method.f) / norm(method.∇f)^2) # Polyak's
    return LM
end

function initialize_bundle_method(;max_iter = max_iter)
    # optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "warm_start_init_point" => "yes")
    # optimizer = optimizer_with_attributes(
    #     CPLEX.Optimizer, 
    #     "CPX_PARAM_SCRIND" => 1, 
    #     "CPX_PARAM_QPMETHOD" => 2
    # )
    optimizer = optimizer_with_attributes(
        OSQP.Optimizer, 
        "verbose" => false, 
        "linsys_solver" => "mkl pardiso",
    )

    # Change parameters
    params = BM.Parameters()
    BM.set_parameter(params, "maxiter", max_iter)
    BM.set_parameter(params, "ϵ_s", 1.e-4)
    BM.set_parameter(params, "ncuts_per_iter", npartitions)

    # Lagrange master method
    return DD.BundleMaster(BM.ProximalMethod, optimizer, params)
    # return DD.BundleMaster(BM.TrustRegionMethod, optimizer, params)
end

# file = ARGS[1]
# npartitions = parse(Int, ARGS[2])
# log_path = ""
# try 
#     global log_path = ARGS[3]
# catch BoundsError
#     global log_path = pwd()
# end

file = "/home/kimk/REPOS/pglib-opf/pglib_opf_case300_ieee.m"
npartitions = 4
log_path = pwd()

main(file, npartitions, log_path,
    max_iter = 3000)
