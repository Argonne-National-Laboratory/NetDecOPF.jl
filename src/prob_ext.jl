mutable struct W_ACRModel <: PM.AbstractWRModel PM.@pm_fields end
mutable struct WU_ACRModel <: PM.AbstractWRModel PM.@pm_fields end

function build_acopf_with_free_lines(pm::W_ACRModel)
    build_opf_mod(pm)
    nw = PM.nw_id_default
    vr = PM.var(pm, nw)[:vr] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :bus)], base_name="$(nw)_vr", start = PM.comp_start_value(PM.ref(pm, nw, :bus, i), "vr_start", 1.0))
    vi = PM.var(pm, nw)[:vi] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :bus)], base_name="$(nw)_vi", start = PM.comp_start_value(PM.ref(pm, nw, :bus, i), "vi_start"))
    for (i, bus) in PM.ref(pm, nw, :bus)
        JuMP.set_lower_bound(vr[i], -bus["vmax"])
        JuMP.set_upper_bound(vr[i],  bus["vmax"])
        JuMP.set_lower_bound(vi[i], -bus["vmax"])
        JuMP.set_upper_bound(vi[i],  bus["vmax"])
    end

    w = PM.var(pm, :w)
    for (i, bus) in PM.ref(pm, nw, :bus)
        JuMP.@constraint(pm.model, w[i] == vr[i]^2 + vi[i]^2)
    end

    wr = PM.var(pm, :wr)
    wi = PM.var(pm, :wi)
    for (i,j) in PM.ids(pm, :buspairs)
        JuMP.@constraint(pm.model, wr[(i,j)] == vr[i] * vr[j] + vi[i] * vi[j])
        JuMP.@constraint(pm.model, wi[(i,j)] == vi[i] * vr[j] - vr[i] * vi[j])
    end
end

function build_acopf_with_free_lines(pm::WU_ACRModel)
    build_opf_mod(pm)
    nw = PM.nw_id_default
    vr = PM.var(pm, nw)[:vr] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :bus)], base_name="$(nw)_vr", start = PM.comp_start_value(PM.ref(pm, nw, :bus, i), "vr_start", 1.0))
    vi = PM.var(pm, nw)[:vi] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :bus)], base_name="$(nw)_vi", start = PM.comp_start_value(PM.ref(pm, nw, :bus, i), "vi_start"))
    for (i, bus) in PM.ref(pm, nw, :bus)
        JuMP.set_lower_bound(vr[i], -bus["vmax"])
        JuMP.set_upper_bound(vr[i],  bus["vmax"])
        JuMP.set_lower_bound(vi[i], -bus["vmax"])
        JuMP.set_upper_bound(vi[i],  bus["vmax"])
    end

    w = PM.var(pm, :w)
    ur = PM.var(pm, nw)[:ur] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :bus)], base_name="$(nw)_ur", start = 1.0)
    ui = PM.var(pm, nw)[:ui] = JuMP.@variable(pm.model, [i in PM.ids(pm, nw, :bus)], base_name="$(nw)_ui", start = 1.0)
    for (i, bus) in PM.ref(pm, nw, :bus)
        JuMP.@constraint(pm.model, w[i] == ur[i] + ui[i])
        JuMP.@constraint(pm.model, ur[i] == vr[i]^2)
        JuMP.@constraint(pm.model, ui[i] == vi[i]^2)
        JuMP.@constraint(pm.model, ur[i] <= +bus["vmax"] * vr[i])
        JuMP.@constraint(pm.model, ur[i] >= -bus["vmax"] * vr[i])
        JuMP.@constraint(pm.model, ui[i] <= +bus["vmax"] * vi[i])
        JuMP.@constraint(pm.model, ui[i] >= -bus["vmax"] * vi[i])

        # SOC relaxation constraints
        # JuMP.@constraint(pm.model, ur[i] >= vr[i]^2)
        # JuMP.@constraint(pm.model, ui[i] >= vi[i]^2)
    end

    urr = PM.var(pm, nw)[:urr] = JuMP.@variable(pm.model, [(i,j) in PM.ids(pm, :buspairs)], base_name="$(nw)_urr", start = 1.0)
    uii = PM.var(pm, nw)[:uii] = JuMP.@variable(pm.model, [(i,j) in PM.ids(pm, :buspairs)], base_name="$(nw)_uii", start = 1.0)
    uir = PM.var(pm, nw)[:uir] = JuMP.@variable(pm.model, [(i,j) in PM.ids(pm, :buspairs)], base_name="$(nw)_uir", start = 1.0)
    uri = PM.var(pm, nw)[:uri] = JuMP.@variable(pm.model, [(i,j) in PM.ids(pm, :buspairs)], base_name="$(nw)_uri", start = 1.0)
    wr = PM.var(pm, :wr)
    wi = PM.var(pm, :wi)
    for (i,j) in PM.ids(pm, :buspairs)
        JuMP.@constraint(pm.model, wr[(i,j)] == urr[(i,j)] + uii[(i,j)])
        JuMP.@constraint(pm.model, wi[(i,j)] == uir[(i,j)] - uri[(i,j)])

        # Gurobi will branch on these guys...
        JuMP.@constraint(pm.model, urr[(i,j)] == vr[i] * vr[j])
        JuMP.@constraint(pm.model, uii[(i,j)] == vi[i] * vi[j])
        JuMP.@constraint(pm.model, uir[(i,j)] == vi[i] * vr[j])
        JuMP.@constraint(pm.model, uri[(i,j)] == vr[i] * vi[j])

        busi = PM.ref(pm, :bus)[i]
        busj = PM.ref(pm, :bus)[j]

        # McCormick relaxation of the bilinear constraints
        # JuMP.@constraint(pm.model, -urr[(i,j)] - busi["vmax"] * vr[j] - busj["vmax"] * vr[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -urr[(i,j)] + busi["vmax"] * vr[j] + busj["vmax"] * vr[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -urr[(i,j)] + busi["vmax"] * vr[j] - busj["vmax"] * vr[i] >= -busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -urr[(i,j)] - busi["vmax"] * vr[j] + busj["vmax"] * vr[i] >= -busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uii[(i,j)] - busi["vmax"] * vi[j] - busj["vmax"] * vi[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uii[(i,j)] + busi["vmax"] * vi[j] + busj["vmax"] * vi[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uii[(i,j)] + busi["vmax"] * vi[j] - busj["vmax"] * vi[i] >= -busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uii[(i,j)] - busi["vmax"] * vi[j] + busj["vmax"] * vi[i] >= -busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uir[(i,j)] - busi["vmax"] * vr[j] - busj["vmax"] * vi[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uir[(i,j)] + busi["vmax"] * vr[j] + busj["vmax"] * vi[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uir[(i,j)] + busi["vmax"] * vr[j] - busj["vmax"] * vi[i] >= -busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uir[(i,j)] - busi["vmax"] * vr[j] + busj["vmax"] * vi[i] >= -busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uri[(i,j)] - busi["vmax"] * vi[j] - busj["vmax"] * vr[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uri[(i,j)] + busi["vmax"] * vi[j] + busj["vmax"] * vr[i] <= +busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uri[(i,j)] + busi["vmax"] * vi[j] - busj["vmax"] * vr[i] >= -busi["vmax"] * busj["vmax"])
        # JuMP.@constraint(pm.model, -uri[(i,j)] - busi["vmax"] * vi[j] + busj["vmax"] * vr[i] >= -busi["vmax"] * busj["vmax"])

        # Simple bounds
        JuMP.@constraint(pm.model, urr[(i,j)] >= -busj["vmax"] * vr[i])
        JuMP.@constraint(pm.model, urr[(i,j)] >= -busi["vmax"] * vr[j])
        JuMP.@constraint(pm.model, urr[(i,j)] <= busj["vmax"] * vr[i])
        JuMP.@constraint(pm.model, urr[(i,j)] <= busi["vmax"] * vr[j])
        JuMP.@constraint(pm.model, uii[(i,j)] >= -busj["vmax"] * vi[i])
        JuMP.@constraint(pm.model, uii[(i,j)] >= -busi["vmax"] * vi[j])
        JuMP.@constraint(pm.model, uii[(i,j)] <= busj["vmax"] * vi[i])
        JuMP.@constraint(pm.model, uii[(i,j)] <= busi["vmax"] * vi[j])
        JuMP.@constraint(pm.model, uir[(i,j)] >= -busj["vmax"] * vi[i])
        JuMP.@constraint(pm.model, uir[(i,j)] >= -busi["vmax"] * vr[j])
        JuMP.@constraint(pm.model, uir[(i,j)] <= busj["vmax"] * vi[i])
        JuMP.@constraint(pm.model, uir[(i,j)] <= busi["vmax"] * vr[j])
        JuMP.@constraint(pm.model, uri[(i,j)] >= -busj["vmax"] * vr[i])
        JuMP.@constraint(pm.model, uri[(i,j)] >= -busi["vmax"] * vi[j])
        JuMP.@constraint(pm.model, uri[(i,j)] <= busj["vmax"] * vr[i])
        JuMP.@constraint(pm.model, uri[(i,j)] <= busi["vmax"] * vi[j])

        # SOC relaxation constraints
        # JuMP.@constraint(pm.model, wr[(i,j)]^2 + wi[(i,j)]^2 <= w[i] * w[j])
    end
end
#############################################################################################
# Legacy Formulation for ACRPowerModel with W matrices of size 2N-by-2N (Wrr & Wii)
#############################################################################################

mutable struct Large_W_ACRModel <: PM.AbstractACRModel PM.@pm_fields end

function build_acopf_with_free_lines(pm::Large_W_ACRModel)
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

function collect_split_vars(pm::Large_W_ACRModel)
    # wrr = PM.var(pm, :wrr)
    # wri = PM.var(pm, :wri)
    # wii = PM.var(pm, :wii)
    p  = PM.var(pm,  :p)
    q  = PM.var(pm,  :q)

    shared_vars_dict = Dict{String, Dict{Tuple, JuMP.VariableRef}}()
    # shared_vars_dict["wrr"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    # shared_vars_dict["wri"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    # shared_vars_dict["wir"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    # shared_vars_dict["wii"] = Dict{Tuple{Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["p"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()
    shared_vars_dict["q"] = Dict{Tuple{Int64, Int64, Int64}, JuMP.VariableRef}()

    cut_arcs_from = PM.ref(pm, :cut_arcs_from)
    for (l,i,j) in cut_arcs_from
        # if !((i,j) in keys(shared_vars_dict["wrr"]))
        #     shared_vars_dict["wrr"][(i,j)] = wrr[i,j]
        #     shared_vars_dict["wri"][(i,j)] = wri[i,j]
        #     shared_vars_dict["wir"][(i,j)] = wri[j,i]
        #     shared_vars_dict["wii"][(i,j)] = wii[i,j]
        # end
        shared_vars_dict["p"][(l,i,j)] = p[(l,i,j)]
        shared_vars_dict["p"][(l,j,i)] = p[(l,j,i)]
        shared_vars_dict["q"][(l,i,j)] = q[(l,i,j)]
        shared_vars_dict["q"][(l,j,i)] = q[(l,j,i)]
    end
    return shared_vars_dict
end

function variable_w_matrix(pm::Large_W_ACRModel; nw::Int=PM.nw_id_default)
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

#############################################################################################
# W_global_ACRModel
#############################################################################################

mutable struct W_global_ACRModel <: PM.AbstractWRModel PM.@pm_fields end

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

    w = PM.var(pm, :w)
    for (i, bus) in PM.ref(pm, nw, :all_bus)
        JuMP.@constraint(pm.model, w[i] == vr[i]^2 + vi[i]^2)
    end

    wr = PM.var(pm, :wr)
    wi = PM.var(pm, :wi)
    for (i,j) in PM.ids(pm, :all_buspairs)
        JuMP.@constraint(pm.model, wr[(i,j)] == vr[i] * vr[j] + vi[i] * vi[j])
        JuMP.@constraint(pm.model, wi[(i,j)] == vi[i] * vr[j] - vr[i] * vi[j])
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

    for i in keys(pm.data["widx_to_bus"]), j in keys(pm.data["widx_to_bus"])
        shared_vars_dict["WI"][(i,j)] = WI[i,j]
        if i <= j
            shared_vars_dict["WR"][(i,j)] = WR[i,j]
        end
    end
    
    return shared_vars_dict
end