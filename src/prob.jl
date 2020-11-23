# This file includes functions that build subnetwork models
# with extensions from PowerModels.jl

function decompose(
    data::Dict{String, Any},
    N_gs::Vector{Vector{Int64}},
    modeltype::Type{T},
    build_function::Function;
    extra_ref_extensions=[]
    )::NetDecModel where T <: PM.AbstractPowerModel

    return decompose(data, [Set(i) for i in N_gs], modeltype, build_function, extra_ref_extensions = extra_ref_extensions)
end

function decompose(
    data::Dict{String, Any},
    N_gs::Vector{Set{Int64}},
    modeltype::Type{T},
    build_function::Function;
    extra_ref_extensions=[]
    )::NetDecModel where T <: PM.AbstractPowerModel

    models = T[]
    shared_vars_dict = Dict{Int64, Dict}()
    for i in eachindex(N_gs)
        N_g = N_gs[i]
        sub_data = generate_subnet_data(data, N_g)
        push!(models, PM.instantiate_model(sub_data, modeltype, build_function, ref_extensions = push!([ref_add_cut_bus!, ref_add_cut_branch!], extra_ref_extensions...)))
        shared_vars_dict[i] = collect_split_vars(models[i])
    end

    return NetDecModel(N_gs, models, shared_vars_dict)
end

# modified from original build_opf
# get rid of ref buses and apply power balance to only nodes in partition
function build_opf_mod(pm::PM.ACRPowerModel)
    # variable_bus_voltage(pm, bounded=false)
    # vr = var(pm, :vr)
    # vi = var(pm, :vi)
    # for (i, bus) in setdiff(ref(pm, :bus), ref(pm, :cut_bus))
    #     JuMP.set_lower_bound(vr[i], -bus["vmax"])
    #     JuMP.set_upper_bound(vr[i],  bus["vmax"])
    #     JuMP.set_lower_bound(vi[i], -bus["vmax"])
    #     JuMP.set_upper_bound(vi[i],  bus["vmax"])
    #     constraint_voltage_magnitude_bounds(pm, i)
    # end
    PM.variable_bus_voltage(pm)
    PM.variable_gen_power(pm)
    PM.variable_branch_power(pm)
    PM.variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm)

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

    for i in PM.ids(pm, :dcline)
        PM.constraint_dcline_power_losses(pm, i)
    end
end

function build_opf_mod(pm::PM.AbstractPowerModel)
    PM.variable_bus_voltage(pm)
    PM.variable_gen_power(pm)
    PM.variable_branch_power(pm)
    PM.variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm)

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

    for i in PM.ids(pm, :dcline)
        PM.constraint_dcline_power_losses(pm, i)
    end
end


#=
function PM.variable_bus_voltage_on_off(pm::PM.AbstractACRModel; kwargs...)
    variable_bus_voltage(pm; kwargs...)
end

function PM.constraint_model_voltage_on_off(pm::PM.AbstractACRModel; kwargs...)
end

function PM.constraint_ohms_yt_from_on_off(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    wrr_ff = var(pm, n, :wrr)[f_bus, f_bus]
    wii_ff = var(pm, n, :wii)[f_bus, f_bus]
    wrr_ft = var(pm, n, :wrr)[f_bus, t_bus]
    wii_ft = var(pm, n, :wii)[f_bus, t_bus]
    wri_ft = var(pm, n, :wri)[f_bus, t_bus]
    wri_tf = var(pm, n, :wri)[t_bus, f_bus]
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_fr == z * ((g+g_fr)/tm^2*(wrr_ff + wii_ff) + (-g*tr+b*ti)/tm^2*(wrr_ft + wii_ft) + (-b*tr-g*ti)/tm^2*(wri_tf - wri_ft)) )
    JuMP.@constraint(pm.model, q_fr == z * (-(b+b_fr)/tm^2*(wrr_ff + wii_ff) - (-b*tr-g*ti)/tm^2*(wrr_ft + wii_ft) + (-g*tr+b*ti)/tm^2*(wri_tf - wri_ft)) )
end
=#

# modified from original build_ots
# get rid of ref buses and apply power balance to only nodes in partition
function build_ots_mod(pm::PM.AbstractPowerModel)
    PM.variable_branch_indicator(pm)
    PM.variable_bus_voltage_on_off(pm)
    PM.variable_gen_power(pm)
    PM.variable_branch_power(pm)
    PM.variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm)

    PM.constraint_model_voltage_on_off(pm)

    for i in setdiff(PM.ids(pm, :bus), PM.ids(pm, :cut_bus))
        PM.constraint_power_balance(pm, i)
    end

    for i in PM.ids(pm, :branch)
        PM.constraint_ohms_yt_from_on_off(pm, i)
        PM.constraint_ohms_yt_to_on_off(pm, i)

        PM.constraint_voltage_angle_difference_on_off(pm, i)

        PM.constraint_thermal_limit_from_on_off(pm, i)
        PM.constraint_thermal_limit_to_on_off(pm, i)
    end

    for i in PM.ids(pm, :dcline)
        PM.constraint_dcline_power_losses(pm, i)
    end
end

# modified from original build_opf_bf
# get rid of ref buses and apply power balance to only nodes in partition
function build_opf_bf_mod(pm::PM.AbstractPowerModel)
    PM.variable_bus_voltage(pm)
    PM.variable_gen_power(pm)
    PM.variable_branch_power(pm)
    PM.variable_branch_current(pm)
    PM.variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm)

    PM.constraint_model_current(pm)

    for i in setdiff(PM.ids(pm, :bus), PM.ids(pm, :cut_bus))
        PM.constraint_power_balance(pm, i)
    end

    for i in PM.ids(pm, :branch)
        PM.constraint_power_losses(pm, i)
        PM.constraint_voltage_magnitude_difference(pm, i)

        PM.constraint_voltage_angle_difference(pm, i)

        PM.constraint_thermal_limit_from(pm, i)
        PM.constraint_thermal_limit_to(pm, i)
    end

    for i in PM.ids(pm, :dcline)
        PM.constraint_dcline_power_losses(pm, i)
    end
end


# modified from original objective_min_fuel_and_flow_cost to allow for
# the case with no generators
function objective_min_fuel_and_flow_cost_mod(pm::PM.AbstractPowerModel; kwargs...)
    model = PM.check_cost_models(pm)
    if model == 1
        return PM.objective_min_fuel_and_flow_cost_pwl(pm; kwargs...)
    elseif model == 2
        return PM.objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    elseif model === nothing
        return JuMP.@objective(pm.model, Min, 0)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end
end

function variable_w_matrix(pm::PM.ACRPowerModel; nw::Int=pm.cnw)
    wrr = PM.var(pm, nw)[:wrr] = JuMP.@variable(pm.model,
        [i in PM.ids(pm, nw, :bus), j in PM.ids(pm, nw, :bus)], base_name="$(nw)_wrr",
        start = PM.comp_start_value(PM.ref(pm, nw, :bus, i), "vr_start", 1.0) * PM.comp_start_value(PM.ref(pm, nw, :bus, j), "vr_start", 1.0)
        )
    wri = PM.var(pm, nw)[:wri] = JuMP.@variable(pm.model,
        [i in PM.ids(pm, nw, :bus), j in PM.ids(pm, nw, :bus)], base_name="$(nw)_wri",
        start = PM.comp_start_value(PM.ref(pm, nw, :bus, i), "vr_start", 1.0) * PM.comp_start_value(PM.ref(pm, nw, :bus, j), "vi_start")
        )
    wii = PM.var(pm, nw)[:wii] = JuMP.@variable(pm.model,
        [i in PM.ids(pm, nw, :bus), j in PM.ids(pm, nw, :bus)], base_name="$(nw)_wii",
        start = PM.comp_start_value(PM.ref(pm, nw, :bus, i), "vi_start") * PM.comp_start_value(PM.ref(pm, nw, :bus, j), "vi_start")
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
    wrr = PM.var(pm, :wrr)
    wri = PM.var(pm, :wri)
    wii = PM.var(pm, :wii)
    vr = PM.var(pm, :vr)
    vi = PM.var(pm, :vi)
    for (i, _) in PM.ref(pm, :bus), (j, _) in PM.ref(pm, :bus)
        JuMP.@constraint(pm.model, wrr[i,j] == vr[i] * vr[j])
        JuMP.@constraint(pm.model, wri[i,j] == vr[i] * vi[j])
        JuMP.@constraint(pm.model, wii[i,j] == vi[i] * vi[j])
    end
end

function build_acopf_with_free_lines(pm::PM.SDPWRMPowerModel)
    build_opf_mod(pm)
end

function build_socwr_with_free_lines(pm::PM.AbstractWRModel)
    build_opf_mod(pm)
end

function build_socbf_with_free_lines(pm::PM.AbstractBFModel)
    build_opf_bf_mod(pm)
end

function build_acots_with_free_lines(pm::PM.AbstractACPModel)
    build_ots_mod(pm)
end

# function build_acots_with_free_lines(pm::PM.AbstractACRModel)
#     variable_w_matrix(pm)
#     build_ots_mod(pm)
# end


function collect_split_vars(pm::PM.ACPPowerModel)
    p = PM.var(pm, :p)
    q = PM.var(pm, :q)
    z = PM.var(pm, :z_branch)
    shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    shared_vars_dict["z"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        shared_vars_dict["z"][(l,i,j)] = z[l]
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end
    return shared_vars_dict
end

function collect_split_vars(pm::PM.ACRPowerModel)
    wrr = PM.var(pm, :wrr)
    wri = PM.var(pm, :wri)
    wii = PM.var(pm, :wii)
    p  = PM.var(pm,  :p)
    q  = PM.var(pm,  :q)

    shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    shared_vars_dict["wrr"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["wri"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["wir"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["wii"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()

    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
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
    wr = PM.var(pm, :wr)
    wi = PM.var(pm, :wi)
    p  = PM.var(pm,  :p)
    q  = PM.var(pm,  :q)

    shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    shared_vars_dict["wr"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["wi"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
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
    end
    return shared_vars_dict
end

function collect_split_vars(pm::PM.AbstractBFModel)
    p  = PM.var(pm,  :p)
    q  = PM.var(pm,  :q)
    ccm = PM.var(pm, :ccm)

    shared_vars_dict = Dict{String, Dict}()
    shared_vars_dict["ccm"] = Dict{Int64, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()

    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        shared_vars_dict["ccm"][l] = ccm[l]
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end
    return shared_vars_dict
end

function collect_split_vars(pm::PM.SDPWRMPowerModel)
    p = PM.var(pm, :p)
    q = PM.var(pm, :q)
    WR = PM.var(pm, :WR)
    WI = PM.var(pm, :WI)
    bus_ids = PM.ids(pm, :bus)
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))

    shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["WR"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["WI"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
        shared_vars_dict["WR"][(i,j)] = WR[lookup_w_index[i], lookup_w_index[j]]
        shared_vars_dict["WI"][(i,j)] = WI[lookup_w_index[i], lookup_w_index[j]]
        shared_vars_dict["WR"][(j,i)] = WR[lookup_w_index[j], lookup_w_index[i]]
        shared_vars_dict["WI"][(j,i)] = WI[lookup_w_index[j], lookup_w_index[i]]
    end
    return shared_vars_dict
end