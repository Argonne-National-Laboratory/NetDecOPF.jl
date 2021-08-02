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
