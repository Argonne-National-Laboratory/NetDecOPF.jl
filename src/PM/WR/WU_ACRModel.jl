mutable struct WU_ACRModel <: PM.AbstractWRModel PM.@pm_fields end

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
