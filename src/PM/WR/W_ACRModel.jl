mutable struct W_ACRModel <: PM.AbstractWRModel PM.@pm_fields end

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
