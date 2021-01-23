# This file implements one solution method for the network decomposition
# using LagrangeDual from DualDecomposition.jl

function init_DD_algo(dn_model::NetDecModel; 
    method::Type{T} = BM.TrustRegionMethod, 
    maxiter::Int64 = 1000,
    tol::Float64 = 1e-6)::DD.LagrangeDual where T <: BM.AbstractMethod
    algo = DD.LagrangeDual(method, maxiter, tol)
    
    partition = get_partition(dn_model)
    models = get_models(dn_model)
    for i in eachindex(partition)
        DD.add_block_model!(algo, i, models[i].model)
    end
    add_split_vars_to_algo!(dn_model, algo)
    return algo
end

function add_split_vars_to_algo!(dn_model::NetDecModel, algo::DD.AbstractMethod)
    coupling_vars = Vector{DD.CouplingVariableRef}()
    partition = get_partition(dn_model)
    split_vars = get_split_vars(dn_model)
    for i in eachindex(partition)
        vars_dict = split_vars[i]
        for (varname, dict) in vars_dict
            for (idx, vref) in dict
                # id = varname * "_" * string(idx)
                id = (varname, idx)
                push!(coupling_vars, DD.CouplingVariableRef(i, id, vref))
            end
        end    
    end
    DD.set_coupling_variables!(algo, coupling_vars)
    return
end
