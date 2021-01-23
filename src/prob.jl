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
    denet_refs = [ref_add_cut_bus!, ref_add_cut_branch!]
    if T <: PM.AbstractWRMModel
        push!(denet_refs, ref_add_global_bus!)
    end
    for i in eachindex(N_gs)
        N_g = N_gs[i]
        sub_data = generate_subnet_data(data, N_g)
        push!(models, PM.instantiate_model(sub_data, modeltype, build_function, ref_extensions = push!(denet_refs, extra_ref_extensions...)))
        shared_vars_dict[i] = collect_split_vars(models[i])
    end

    return NetDecModel(N_gs, models, shared_vars_dict)
end

# modified from original build_opf
# get rid of ref buses and apply power balance to only nodes in partition
function build_opf_mod(pm::PM.AbstractPowerModel)
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

function build_opf_mod(pm::PM.AbstractWRMModel)
    function variable_bus_voltage(pm::PM.AbstractWRMModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
        wr_min, wr_max, wi_min, wi_max = PM.ref_calc_voltage_product_bounds(PM.ref(pm, nw, :all_buspairs))
        bus_ids = PM.ids(pm, nw, :all_bus)
    
        w_index = 1:length(bus_ids)
        lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))
    
        WR_start = zeros(length(bus_ids), length(bus_ids)) + I
    
        WR = PM.var(pm, nw)[:WR] = JuMP.@variable(pm.model,
            [i=1:length(bus_ids), j=1:length(bus_ids)], Symmetric, base_name="$(nw)_WR", start=WR_start[i,j]
        )
        if report
            PM.sol(pm, nw)[:WR] = WR
        end
    
        WI = PM.var(pm, nw)[:WI] = JuMP.@variable(pm.model,
            [1:length(bus_ids), 1:length(bus_ids)], base_name="$(nw)_WI", start=0.0
        )
        if report
            PM.sol(pm, nw)[:WI] = WI
        end
    
        # bounds on diagonal
        for (i, bus) in PM.ref(pm, nw, :all_bus)
            w_idx = lookup_w_index[i]
            wr_ii = WR[w_idx,w_idx]
            wi_ii = WR[w_idx,w_idx]
    
            if bounded
                JuMP.set_lower_bound(wr_ii, (bus["vmin"])^2)
                JuMP.set_upper_bound(wr_ii, (bus["vmax"])^2)    
            else
                JuMP.set_lower_bound(wr_ii, 0)
            end
        end
    
        # bounds on off-diagonal
        for (i,j) in PM.ids(pm, nw, :all_buspairs)
            wi_idx = lookup_w_index[i]
            wj_idx = lookup_w_index[j]
    
            if bounded
                JuMP.set_upper_bound(WR[wi_idx, wj_idx], wr_max[(i,j)])
                JuMP.set_lower_bound(WR[wi_idx, wj_idx], wr_min[(i,j)])
    
                JuMP.set_upper_bound(WI[wi_idx, wj_idx], wi_max[(i,j)])
                JuMP.set_lower_bound(WI[wi_idx, wj_idx], wi_min[(i,j)])
            end
        end
    
        PM.var(pm, nw)[:w] = Dict{Int,Any}()
        for (i, bus) in PM.ref(pm, nw, :all_bus)
            w_idx = lookup_w_index[i]
            PM.var(pm, nw, :w)[i] = WR[w_idx,w_idx]
        end
        report && PM._IM.sol_component_value(pm, nw, :all_bus, :w, PM.ids(pm, nw, :all_bus), PM.var(pm, nw)[:w])
    
        PM.var(pm, nw)[:wr] = Dict{Tuple{Int,Int},Any}()
        PM.var(pm, nw)[:wi] = Dict{Tuple{Int,Int},Any}()
        for (i,j) in PM.ids(pm, nw, :all_buspairs)
            w_fr_index = lookup_w_index[i]
            w_to_index = lookup_w_index[j]
    
            PM.var(pm, nw, :wr)[(i,j)] = WR[w_fr_index, w_to_index]
            PM.var(pm, nw, :wi)[(i,j)] = WI[w_fr_index, w_to_index]
        end
    end
    function constraint_voltage_angle_difference_all(pm::PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
        branch = PM.ref(pm, nw, :all_branch, i)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        pair = (f_bus, t_bus)
        buspair = PM.ref(pm, nw, :all_buspairs, pair)
    
        if buspair["branch"] == i
            PM.constraint_voltage_angle_difference(pm, nw, f_idx, buspair["angmin"], buspair["angmax"])
        end
    end
    variable_bus_voltage(pm) # re-define W matrix
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

        # PM.constraint_voltage_angle_difference(pm, i)

        PM.constraint_thermal_limit_from(pm, i)
        PM.constraint_thermal_limit_to(pm, i)
    end

    for i in PM.ids(pm, :all_branch)
        constraint_voltage_angle_difference_all(pm, i)
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
    build_opf_mod(pm)
    PM.variable_bus_voltage_magnitude_sqr(pm)
    PM.variable_buspair_voltage_product(pm)

    w  = PM.var(pm,  :w)
    wr = PM.var(pm, :wr)
    wi = PM.var(pm, :wi)
    vr = PM.var(pm, :vr)
    vi = PM.var(pm, :vi)

    for (i, bus) in PM.ref(pm, :bus)
        JuMP.@constraint(pm.model, w[i] == vr[i]^2 + vi[i]^2)
    end

    for (_, branch) in PM.ref(pm, :branch)
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
    # shared_vars_dict["wr"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    # shared_vars_dict["wi"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()

    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        # if !((i,j) in keys(shared_vars_dict["wr"]))
        #     shared_vars_dict["wr"][(i,j)] = wr[(i,j)]
        #     shared_vars_dict["wi"][(i,j)] = wi[(i,j)]
        # end
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

function collect_split_vars(pm::PM.AbstractWRMModel)
    p = PM.var(pm, :p)
    q = PM.var(pm, :q)
    WR = PM.var(pm, :WR)
    WI = PM.var(pm, :WI)
    bus_ids = PM.ids(pm, :all_bus) # all buses from the full network 
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))

    shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["WR"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["WI"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from # penalizing cut lines
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end
    bus = PM.ids(pm, :bus)
    for i in bus_ids, j in bus_ids
        shared_vars_dict["WR"][(i,j)] = WR[lookup_w_index[i], lookup_w_index[j]]
        shared_vars_dict["WI"][(i,j)] = WI[lookup_w_index[i], lookup_w_index[j]]
        if i != j
            shared_vars_dict["WR"][(j,i)] = WR[lookup_w_index[j], lookup_w_index[i]]
            shared_vars_dict["WI"][(j,i)] = WI[lookup_w_index[j], lookup_w_index[i]]
        end
    end
    return shared_vars_dict
end

######################################################
###             EXPERIMENTAL FUNCTIONS             ###
######################################################

function constraint_model_voltage_mod(pm::PM.AbstractPowerModel)
    PM.constraint_model_voltage(pm)
end

# function constraint_model_voltage_mod(pm::PM.AbstractWRMModel)
#     PM._check_missing_keys(PM.var(pm), [:WR,:WI], typeof(pm))

#     WR = PM.var(pm)[:WR]
#     WI = PM.var(pm)[:WI]

#     bus_ids = PM.ids(pm, :bus)
#     lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))

#     local_buses = setdiff(PM.ids(pm, :bus), PM.ids(pm, :cut_bus))
#     local_w_idx = [lookup_w_index[i] for i in local_buses]

#     WR_local = WR[local_w_idx, local_w_idx]
#     WI_local = WI[local_w_idx, local_w_idx]
#     JuMP.@constraint(pm.model, [WR_local WI_local; -WI_local WR_local] in JuMP.PSDCone())
# end