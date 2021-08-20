function build_opf_mod(pm::PM.SparseSDPWRMPowerModel)
    PM.variable_bus_voltage(pm)
    PM.variable_gen_power(pm)
    PM.variable_branch_power(pm)
    PM.variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm) # allow for no generator case

    PM.constraint_model_voltage(pm)

    for i in setdiff(PM.ids(pm, :bus), PM.ids(pm, :cut_bus))
        PM.constraint_power_balance(pm, i)
    end

    for i in PM.ids(pm, :branch)
        PM.constraint_ohms_yt_from(pm, i)
        PM.constraint_ohms_yt_to(pm, i)

        PM.constraint_voltage_angle_difference(pm, i)

        PM.constraint_thermal_limit_from(pm, i)
        PM.constraint_thermal_limit_to(pm, i)
    end

    # for i in PM.ids(pm, :all_branch)
    #     constraint_voltage_angle_difference_all(pm, i)
    # end

    for i in PM.ids(pm, :dcline)
        PM.constraint_dcline_power_losses(pm, i)
    end
end

function collect_split_vars(pm::PM.SparseSDPWRMPowerModel)
    wr = PM.var(pm, :wr)
    wi = PM.var(pm, :wi)
    w  = PM.var(pm,  :w)
    p  = PM.var(pm,  :p)
    q  = PM.var(pm,  :q)

    shared_vars_dict = Dict{String, Dict{Any, JuMP.VariableRef}}()
    shared_vars_dict["wr"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["wi"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["w"] = Dict{Int64, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()

    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        if !((i,j) in keys(shared_vars_dict["wr"]))
            shared_vars_dict["wr"][(i,j)] = wr[(i,j)]
            shared_vars_dict["wi"][(i,j)] = wi[(i,j)]
        end
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
        shared_vars_dict["w"][i] = w[i]
        shared_vars_dict["w"][j] = w[j]
    end
    return shared_vars_dict
end
