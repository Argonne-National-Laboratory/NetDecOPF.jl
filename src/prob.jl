# This file includes functions that build subnetwork models
# with extensions from PowerModels.jl

function decompose(
    data::Dict{String, Any},
    N_gs::Vector{Vector{Int64}},
    modeltype::Type{T},
    build_function::Function;
    )::DeNetModel where T <: PM.AbstractPowerModel

    return decompose(data, [Set(i) for i in N_gs], modeltype, build_function)
end

function decompose(
    data::Dict{String, Any},
    N_gs::Vector{Set{Int64}},
    modeltype::Type{T},
    build_function::Function;
    )::DeNetModel where T <: PM.AbstractPowerModel

    models = T[]
    shared_vars_dict = Dict{Int64, Dict}()
    for i in eachindex(N_gs)
        N_g = N_gs[i]
        sub_data = generate_subnet_data(data, N_g)
        push!(models, instantiate_model(sub_data, modeltype, build_function, ref_extensions = [ref_add_cut_bus!, ref_add_cut_branch!]))
        shared_vars_dict[i] = collect_split_vars(models[i])
    end

    return DeNetModel(N_gs, models, shared_vars_dict)
end

# modified from original build_opf_bf
# get rid of ref buses and apply power balance to only nodes in partition
function build_opf_mod(pm::PM.AbstractPowerModel)
    variable_bus_voltage(pm, bounded=false)
    vr = var(pm, :vr)
    vi = var(pm, :vi)
    for (i, bus) in setdiff(ref(pm, :bus), ref(pm, :cut_bus))
        JuMP.set_lower_bound(vr[i], -bus["vmax"])
        JuMP.set_upper_bound(vr[i],  bus["vmax"])
        JuMP.set_lower_bound(vi[i], -bus["vmax"])
        JuMP.set_upper_bound(vi[i],  bus["vmax"])
        constraint_voltage_magnitude_bounds(pm, i)
    end
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm)

    constraint_model_voltage(pm)

    for i in setdiff(ids(pm, :bus), ids(pm, :cut_bus))
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

# modified from original build_opf_bf
# get rid of ref buses and apply power balance to only nodes in partition
function build_opf_bf_mod(pm::PM.AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_branch_current(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm)

    constraint_model_current(pm)

    for i in setdiff(ids(pm, :bus), ids(pm, :cut_bus))
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_power_losses(pm, i)
        constraint_voltage_magnitude_difference(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end


# modified from original objective_min_fuel_and_flow_cost to allow for
# the case with no generators
function objective_min_fuel_and_flow_cost_mod(pm::PM.AbstractPowerModel; kwargs...)
    model = check_cost_models(pm)
    if model == 1
        return objective_min_fuel_and_flow_cost_pwl(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    elseif model === nothing
        return JuMP.@objective(pm.model, Min, 0)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end
end

function variable_w_matrix(pm::PM.ACRPowerModel; nw::Int=pm.cnw)
    wrr = var(pm, nw)[:wrr] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus), j in ids(pm, nw, :bus)], base_name="$(nw)_wrr",
        start = comp_start_value(ref(pm, nw, :bus, i), "vr_start", 1.0) * comp_start_value(ref(pm, nw, :bus, j), "vr_start", 1.0)
        )
    wri = var(pm, nw)[:wri] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus), j in ids(pm, nw, :bus)], base_name="$(nw)_wri",
        start = comp_start_value(ref(pm, nw, :bus, i), "vr_start", 1.0) * comp_start_value(ref(pm, nw, :bus, j), "vi_start")
        )
    wii = var(pm, nw)[:wii] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus), j in ids(pm, nw, :bus)], base_name="$(nw)_wii",
        start = comp_start_value(ref(pm, nw, :bus, i), "vi_start") * comp_start_value(ref(pm, nw, :bus, j), "vi_start")
        )
end

#=
function build_acopf_with_free_lines(pm::PM.AbstractPowerModel)
    # this is for adding w variables to the model
    variable_bus_voltage_magnitude_sqr(pm)
    variable_buspair_voltage_product(pm)

    w  = var(pm,  :w)
    wr = var(pm, :wr)
    wi = var(pm, :wi)
    vr = var(pm, :vr)
    vi = var(pm, :vi)

    for (i, bus) in ref(pm, :bus)
        JuMP.@constraint(pm.model, w[i] == vr[i]^2 + vi[i]^2)
    end

    for (_, branch) in ref(pm, :branch)
        fbus = branch["f_bus"]
        tbus = branch["t_bus"]
        JuMP.@constraint(pm.model, wr[(fbus, tbus)] == vr[fbus] * vr[tbus] + vi[fbus] * vi[tbus])
        JuMP.@constraint(pm.model, wi[(fbus, tbus)] == vi[fbus] * vr[tbus] - vr[fbus] * vi[tbus])
    end
end
=#

function build_acopf_with_free_lines(pm::PM.AbstractPowerModel)
    build_opf_mod(pm)
    variable_w_matrix(pm)
    wrr = var(pm, :wrr)
    wri = var(pm, :wri)
    wii = var(pm, :wii)
    vr = var(pm, :vr)
    vi = var(pm, :vi)
    for (i, _) in ref(pm, :bus), (j, _) in ref(pm, :bus)
        JuMP.@constraint(pm.model, wrr[i,j] == vr[i] * vr[j])
        JuMP.@constraint(pm.model, wri[i,j] == vr[i] * vi[j])
        JuMP.@constraint(pm.model, wii[i,j] == vi[i] * vi[j])
    end
end

function build_socwr_with_free_lines(pm::PM.AbstractWRModel)
    build_opf_mod(pm)
end

function build_socbf_with_free_lines(pm::PM.AbstractBFModel)
    build_opf_bf_mod(pm)
end

function collect_split_vars(pm::PM.ACRPowerModel)
    wrr = var(pm, :wrr)
    wri = var(pm, :wri)
    wii = var(pm, :wii)
    p  = var(pm,  :p)
    q  = var(pm,  :q)

    shared_vars_dict = Dict{String, Dict{Tuple, VariableRef}}()
    shared_vars_dict["wrr"] = Dict{Tuple{Int64, Int64}, VariableRef}()
    shared_vars_dict["wri"] = Dict{Tuple{Int64, Int64}, VariableRef}()
    shared_vars_dict["wir"] = Dict{Tuple{Int64, Int64}, VariableRef}()
    shared_vars_dict["wii"] = Dict{Tuple{Int64, Int64}, VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()

    cut_arcs_from = ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        if !((i,j) in keys(shared_vars_dict["wrr"]))
            shared_vars_dict["wrr"][(i,j)] = wrr[i,j]
            shared_vars_dict["wri"][(i,j)] = wri[i,j]
            shared_vars_dict["wir"][(i,j)] = wri[j,i]
            shared_vars_dict["wii"][(i,j)] = wii[i,j]
        end
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end
    return shared_vars_dict
end

function collect_split_vars(pm::PM.AbstractWRModel)
    wr = var(pm, :wr)
    wi = var(pm, :wi)
    p  = var(pm,  :p)
    q  = var(pm,  :q)

    shared_vars_dict = Dict{String, Dict{Tuple, VariableRef}}()
    shared_vars_dict["wr"] = Dict{Tuple{Int64, Int64}, VariableRef}()
    shared_vars_dict["wi"] = Dict{Tuple{Int64, Int64}, VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()

    cut_arcs_from = ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        if !((i,j) in keys(shared_vars_dict["wr"]))
            shared_vars_dict["wr"][(i,j)] = wr[(i,j)]
            shared_vars_dict["wi"][(i,j)] = wi[(i,j)]
        end
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end
    return shared_vars_dict
end

function collect_split_vars(pm::PM.AbstractBFModel)
    p  = var(pm,  :p)
    q  = var(pm,  :q)
    ccm = var(pm, :ccm)

    shared_vars_dict = Dict{String, Dict}()
    shared_vars_dict["ccm"] = Dict{Int64, VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()

    cut_arcs_from = ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        shared_vars_dict["ccm"][l] = ccm[l]
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end
    return shared_vars_dict
end
