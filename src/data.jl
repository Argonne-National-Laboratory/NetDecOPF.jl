# This file includes all functions that serve to create data structures that 
# record partition information and relevant network parameters

function find_neighbor_buses(data::Dict{String, Any}, N_g::Set{Int64})::Set{Int64}
    neighbors = Set{Int64}()
    for line in keys(data["branch"])
        if data["branch"][line]["f_bus"] in N_g || data["branch"][line]["t_bus"] in N_g
            push!(neighbors, data["branch"][line]["f_bus"], data["branch"][line]["t_bus"])
        end
    end
    return setdiff!(neighbors, N_g)
end

function generate_subnet_data(data::Dict{String, Any}, N_g::Vector{Int64})
    return generate_subnet_data(data, Set(N_g))
end

function generate_subnet_data(data::Dict{String, Any}, N_g::Set{Int64})
    sub_data = deepcopy(data)
    nbus = find_neighbor_buses(data, N_g)
    N = Set(parse.(Int, collect(keys(data["bus"]))))
    # deactivate unrelated buses
    for i in setdiff!(N, N_g, nbus)
        sub_data["bus"]["$(i)"]["bus_type"] = PM.pm_component_status_inactive["bus"]
    end
    # deactivate lines between neighbor buses
    for line in keys(data["branch"])
        if sub_data["branch"][line]["f_bus"] in nbus && data["branch"][line]["t_bus"] in nbus
            sub_data["branch"][line]["br_status"] = PM.pm_component_status_inactive["branch"]
        end
    end
    # deactivate unrelated generators (for those at neighbor nodes)
    for (_, gen) in sub_data["gen"]
        if !(gen["gen_bus"] in N_g)
            gen["gen_status"] = PM.pm_component_status_inactive["gen"]
        end
    end
    sub_data["cut_bus"] = nbus
    return sub_data
end

function ref_add_cut_bus!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
    for (nw, nw_ref) in ref[:nw]
        nw_ref[:cut_bus] = Dict(i => data["bus"]["$(i)"] for i in data["cut_bus"])
    end
end

function ref_add_cut_branch!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
    for (nw, nw_ref) in ref[:nw]
        # set up cut branches and arcs
        nw_ref[:cut_branch] = Dict(parse(Int, x.first) => x.second for x in data["branch"] if
            (x.second["f_bus"] in keys(nw_ref[:bus])) + (x.second["t_bus"] in keys(nw_ref[:bus])) +
            (x.second["f_bus"] in keys(nw_ref[:cut_bus])) + (x.second["t_bus"] in keys(nw_ref[:cut_bus])) == 3)
        nw_ref[:cut_arcs_from] = [(i,branch["f_bus"],branch["t_bus"]) for (i,branch) in nw_ref[:cut_branch]]
        nw_ref[:cut_arcs_to]   = [(i,branch["t_bus"],branch["f_bus"]) for (i,branch) in nw_ref[:cut_branch]]
        nw_ref[:cut_arcs] = [nw_ref[:cut_arcs_from]; nw_ref[:cut_arcs_to]]
    end
end

#=
function ref_add_cut_bus_arcs_refs!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
    for (nw, nw_ref) in ref[:nw]
        cut_bus_arcs = Dict((i, Tuple{Int,Int,Int}[]) for (i,bus) in merge(nw_ref[:bus], nw_ref[:cut_bus]))
        for (l,i,j) in nw_ref[:cut_arcs]
            push!(cut_bus_arcs[i], (l,i,j))
        end
        nw_ref[:cut_bus_arcs] = cut_bus_arcs
    end
end
=#