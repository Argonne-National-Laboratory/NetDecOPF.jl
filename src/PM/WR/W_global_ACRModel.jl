#############################################################################################
# W_global_ACRModel
#############################################################################################

mutable struct W_global_ACRModel <: PM.AbstractWRModel PM.@pm_fields end

"""
    build_opf_mod

This creates additional variables:
- `WR[i,j]` for each bus i, j
- `WI[i,j]` for each bus i, j
"""
function PM.variable_bus_voltage(pm::W_global_ACRModel; nw::Int=PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    wr_min, wr_max, wi_min, wi_max = PM.ref_calc_voltage_product_bounds(PM.ref(pm, nw, :all_buspairs))
    bus_ids = PM.ids(pm, nw, :all_bus)

    w_index = 1:length(bus_ids)
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))
    pm.data["bus_to_widx"] = lookup_w_index
    pm.data["widx_to_bus"] = Dict((i,bi) for (i,bi) in enumerate(bus_ids))

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
            
            ##
            JuMP.set_upper_bound(WI[wj_idx, wi_idx], -wi_min[(i,j)])
            JuMP.set_lower_bound(WI[wj_idx, wi_idx], -wi_max[(i,j)])
        end
    end

    PM.var(pm, nw)[:w] = Dict{Int,Any}()
    for (i, bus) in PM.ref(pm, nw, :all_bus)
        w_idx = lookup_w_index[i]
        PM.var(pm, nw, :w)[i] = WR[w_idx,w_idx]
    end
    report && PM.sol_component_value(pm, nw, :all_bus, :w, PM.ids(pm, nw, :all_bus), PM.var(pm, nw)[:w])

    PM.var(pm, nw)[:wr] = Dict{Tuple{Int,Int},Any}()
    PM.var(pm, nw)[:wi] = Dict{Tuple{Int,Int},Any}()
    for (i,j) in PM.ids(pm, nw, :all_buspairs)
        w_fr_index = lookup_w_index[i]
        w_to_index = lookup_w_index[j]

        PM.var(pm, nw, :wr)[(i,j)] = WR[w_fr_index, w_to_index]
        PM.var(pm, nw, :wi)[(i,j)] = WI[w_fr_index, w_to_index]
    end
end

# function PM.variable_bus_voltage_magnitude_sqr(pm::W_global_ACRModel; nw::Int=PM.nw_id_default, bounded::Bool=true, report::Bool=true)
#     w = PM.var(pm, nw)[:w] = JuMP.@variable(pm.model,
#         [i in PM.ids(pm, nw, :all_bus)], base_name="$(nw)_w",
#         lower_bound = 0.0,
#         start = PM.comp_start_value(PM.ref(pm, nw, :all_bus, i), "w_start", 1.001)
#     )

#     if bounded
#         for (i, bus) in PM.ref(pm, nw, :all_bus)
#             JuMP.set_lower_bound(w[i], bus["vmin"]^2)
#             JuMP.set_upper_bound(w[i], bus["vmax"]^2)
#         end
#     end

#     report && PM.sol_component_value(pm, nw, :all_bus, :w, PM.ids(pm, nw, :all_bus), w)
# end

# function PM.variable_buspair_voltage_product(pm::W_global_ACRModel; nw::Int=PM.nw_id_default, bounded::Bool=true, report::Bool=true)
#     wr = PM.var(pm, nw)[:wr] = JuMP.@variable(pm.model,
#         [bp in PM.ids(pm, nw, :all_buspairs)], base_name="$(nw)_wr",
#         start = PM.comp_start_value(PM.ref(pm, nw, :all_buspairs, bp), "wr_start", 1.0)
#     )
#     wi = PM.var(pm, nw)[:wi] = JuMP.@variable(pm.model,
#         [bp in PM.ids(pm, nw, :all_buspairs)], base_name="$(nw)_wi",
#         start = PM.comp_start_value(PM.ref(pm, nw, :all_buspairs, bp), "wi_start")
#     )

#     if bounded
#         wr_min, wr_max, wi_min, wi_max = PM.ref_calc_voltage_product_bounds(PM.ref(pm, nw, :all_buspairs))

#         for bp in PM.ids(pm, nw, :all_buspairs)
#             JuMP.set_lower_bound(wr[bp], wr_min[bp])
#             JuMP.set_upper_bound(wr[bp], wr_max[bp])

#             JuMP.set_lower_bound(wi[bp], wi_min[bp])
#             JuMP.set_upper_bound(wi[bp], wi_max[bp])
#         end
#     end

#     report && PM.sol_component_value_buspair(pm, nw, :all_buspairs, :wr, PM.ids(pm, nw, :all_buspairs), wr)
#     report && PM.sol_component_value_buspair(pm, nw, :all_buspairs, :wi, PM.ids(pm, nw, :all_buspairs), wi)
# end

"""
    build_acopf_with_free_lines

This is the main function to create the model. 
    
The additional variables are added:
- `vr[i]` for each bus `i`
- `vi[i]` for each bus `i`

The following constraints are added:
- `WI[i,j] == vi[i] * vr[j] - vr[i] * vi[j]` for each pair `i,j`
- `WR[i,j] == vr[i] * vr[j] + vi[i] * vi[j]` for each pair `i,j` such that `i <= j`
"""
function build_acopf_with_free_lines(pm::W_global_ACRModel)
    build_opf_mod(pm)
    nw = PM.nw_id_default
    vr = PM.var(pm, nw)[:vr] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :all_bus)], base_name="$(nw)_vr", start = PM.comp_start_value(PM.ref(pm, nw, :all_bus, i), "vr_start", 1.0))
    vi = PM.var(pm, nw)[:vi] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :all_bus)], base_name="$(nw)_vi", start = PM.comp_start_value(PM.ref(pm, nw, :all_bus, i), "vi_start"))
    for (i, bus) in PM.ref(pm, nw, :all_bus)
        JuMP.set_lower_bound(vr[i], -bus["vmax"])
        JuMP.set_upper_bound(vr[i],  bus["vmax"])
        JuMP.set_lower_bound(vi[i], -bus["vmax"])
        JuMP.set_upper_bound(vi[i],  bus["vmax"])
    end

    # w = PM.var(pm, :w)
    # for (i, bus) in PM.ref(pm, nw, :all_bus)
    #     JuMP.@constraint(pm.model, w[i] == vr[i]^2 + vi[i]^2)
    # end

    # wr = PM.var(pm, :wr)
    # wi = PM.var(pm, :wi)
    # for (i,j) in PM.ids(pm, :all_buspairs)
    #     JuMP.@constraint(pm.model, wr[(i,j)] == vr[i] * vr[j] + vi[i] * vi[j])
    #     JuMP.@constraint(pm.model, wi[(i,j)] == vi[i] * vr[j] - vr[i] * vi[j])
    # end
    WR = PM.var(pm, :WR)
    WI = PM.var(pm, :WI)
    for i in keys(pm.data["widx_to_bus"]), j in keys(pm.data["widx_to_bus"])
        b_i_index = pm.data["widx_to_bus"][i]
        b_j_index = pm.data["widx_to_bus"][j]
        JuMP.@constraint(pm.model, WI[i,j] == vi[b_i_index] * vr[b_j_index] - vr[b_i_index] * vi[b_j_index])
        if i <= j
            JuMP.@constraint(pm.model, WR[i,j] == vr[b_i_index] * vr[b_j_index] + vi[b_i_index] * vi[b_j_index])
        end
    end
end

# function collect_split_vars(pm::W_global_ACRModel)
#     wr = PM.var(pm, :wr)
#     wi = PM.var(pm, :wi)
#     w  = PM.var(pm,  :w)
#     p  = PM.var(pm,  :p)
#     q  = PM.var(pm,  :q)

#     # shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
#     shared_vars_dict = Dict{String, Dict{Any, JuMP.VariableRef}}()
#     shared_vars_dict["wr"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
#     shared_vars_dict["wi"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
#     shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
#     shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
#     shared_vars_dict["vr"] = Dict{Int64, JuMP.VariableRef}()
#     shared_vars_dict["vi"] = Dict{Int64, JuMP.VariableRef}()
#     shared_vars_dict["w"] = Dict{Int64, JuMP.VariableRef}()

#     cut_arcs_from = PM.ref(pm, :cut_arcs_from)
#     for (l,i,j) in cut_arcs_from
#         if !((i,j) in keys(shared_vars_dict["wr"]))
#             shared_vars_dict["wr"][(i,j)] = wr[(i,j)]
#             shared_vars_dict["wi"][(i,j)] = wi[(i,j)]
#         end
#         shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
#         shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
#         shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
#         shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
#     end

#     all_arcs_from = PM.ref(pm, :all_arcs_from)
#     for (l,i,j) in all_arcs_from
#         if !((i,j) in keys(shared_vars_dict["wr"]))
#             shared_vars_dict["wr"][(i,j)] = wr[(i,j)]
#             shared_vars_dict["wi"][(i,j)] = wi[(i,j)]
#         end
#     end

#     vr = PM.var(pm, :vr)
#     vi = PM.var(pm, :vi)
#     buses = PM.ref(pm, :all_bus)
#     for (i,_) in buses
#         shared_vars_dict["vr"][i] = vr[i]
#         shared_vars_dict["vi"][i] = vi[i]
#         shared_vars_dict["w"][i] = w[i]
#     end
    
#     return shared_vars_dict
# end

function collect_split_vars(pm::W_global_ACRModel)
    WR = PM.var(pm, :WR)
    WI = PM.var(pm, :WI)
    p  = PM.var(pm,  :p)
    q  = PM.var(pm,  :q)

    # shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    shared_vars_dict = Dict{String, Dict{Any, JuMP.VariableRef}}()
    shared_vars_dict["WR"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["WI"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()

    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end

    # coupling the complete W matrix
    for i in keys(pm.data["widx_to_bus"]), j in keys(pm.data["widx_to_bus"])
        shared_vars_dict["WI"][(i,j)] = WI[i,j]
        if i <= j
            shared_vars_dict["WR"][(i,j)] = WR[i,j]
        end
    end
    
    # for (_,i,j) in PM.ref(pm, :all_arcs_from)
    #     wi = pm.data["bus_to_widx"][i]
    #     wj = pm.data["bus_to_widx"][j]
    #     shared_vars_dict["WI"][(wi,wj)] = WI[wi,wj]
    #     shared_vars_dict["WI"][(wj,wi)] = WI[wj,wi]
    #     # if (wi,wi) in keys(shared_vars_dict["WR"])
    #     #     shared_vars_dict["WR"][(wi,wi)] = WR[wi,wi]
    #     # end
    #     # if (wj,wj) in keys(shared_vars_dict["WR"])
    #     #     shared_vars_dict["WR"][(wj,wj)] = WR[wj,wj]
    #     # end
    #     if wi <= wj
    #         shared_vars_dict["WR"][(wi,wj)] = WR[wi,wj]
    #     else
    #         shared_vars_dict["WR"][(wj,wi)] = WR[wj,wi]
    #     end
    # end

    return shared_vars_dict
end
