mutable struct NonconvexACRModel <: PM.AbstractWRModel PM.@pm_fields end

function PM.constraint_model_voltage(pm::NonconvexACRModel, n::Int)
    PM._check_missing_keys(PM.var(pm, n), [:w,:wr,:wi], typeof(pm))

    w  = PM.var(pm, n,  :w)
    wr = PM.var(pm, n, :wr)
    wi = PM.var(pm, n, :wi)

    for (i,j) in PM.ids(pm, n, :buspairs)
        JuMP.@constraint(pm.model, wr[(i,j)]^2 + wi[(i,j)]^2 == w[i] * w[j])
    end
end

function collect_split_vars(pm::NonconvexACRModel)
    wr = PM.var(pm, :wr)
    wi = PM.var(pm, :wi)
    # w  = PM.var(pm,  :w)
    p  = PM.var(pm,  :p)
    q  = PM.var(pm,  :q)

    shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    shared_vars_dict["wr"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["wi"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()

    for (i,j) in PM.ids(pm, :buspairs)
        shared_vars_dict["wr"][(i,j)] = wr[(i,j)]
        shared_vars_dict["wi"][(i,j)] = wi[(i,j)]
    end

    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end

    return shared_vars_dict
end
