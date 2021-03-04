# This file includes all functions that serve to create data structures that 
# record partition information and relevant network parameters

# NOTE: this file is updated to incorporate the case when soem lines are inactivated while parsed by PowerModels.jl
#       however, we are still assuming that all buses are active (which might cause errors if not true)

function find_neighbor_buses(data::Dict{String, Any}, N_g::Set{Int64})::Set{Int64}
    neighbors = Set{Int64}()
    for line in keys(data["branch"])
        if data["branch"][line]["br_status"] != PM.pm_component_status_inactive["branch"] && 
                (data["branch"][line]["f_bus"] in N_g || data["branch"][line]["t_bus"] in N_g)
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
            x.second["br_status"] != PM.pm_component_status_inactive["branch"] &&
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

function ref_add_global_bus!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
    function calc_all_buspair_parameters(buses, branches, conductor_ids, ismulticondcutor)
        bus_lookup = Dict(bus["index"] => bus for (i,bus) in buses)    
        branch_lookup = Dict(branch["index"] => branch for (i,branch) in branches)
        buspair_indexes = Set((branch["f_bus"], branch["t_bus"]) for (i,branch) in branch_lookup)
        bp_branch = Dict((bp, typemax(Int64)) for bp in buspair_indexes)
    
        if ismulticondcutor
            bp_angmin = Dict((bp, [-Inf for c in conductor_ids]) for bp in buspair_indexes)
            bp_angmax = Dict((bp, [ Inf for c in conductor_ids]) for bp in buspair_indexes)
        else
            @assert(length(conductor_ids) == 1)
            bp_angmin = Dict((bp, -Inf) for bp in buspair_indexes)
            bp_angmax = Dict((bp,  Inf) for bp in buspair_indexes)
        end
    
        for (l,branch) in branch_lookup
            i = branch["f_bus"]
            j = branch["t_bus"]
    
            if ismulticondcutor
                for c in conductor_ids
                    bp_angmin[(i,j)][c] = max(bp_angmin[(i,j)][c], branch["angmin"][c])
                    bp_angmax[(i,j)][c] = min(bp_angmax[(i,j)][c], branch["angmax"][c])
                end
            else
                bp_angmin[(i,j)] = max(bp_angmin[(i,j)], branch["angmin"])
                bp_angmax[(i,j)] = min(bp_angmax[(i,j)], branch["angmax"])
            end
    
            bp_branch[(i,j)] = min(bp_branch[(i,j)], l)
        end
    
        buspairs = Dict((i,j) => Dict(
            "branch"=>bp_branch[(i,j)],
            "angmin"=>bp_angmin[(i,j)],
            "angmax"=>bp_angmax[(i,j)],
            "tap"=>branch_lookup[bp_branch[(i,j)]]["tap"],
            "vm_fr_min"=>bus_lookup[i]["vmin"],
            "vm_fr_max"=>bus_lookup[i]["vmax"],
            "vm_to_min"=>bus_lookup[j]["vmin"],
            "vm_to_max"=>bus_lookup[j]["vmax"]
            ) for (i,j) in buspair_indexes
        )
    
        # add optional parameters
        for bp in buspair_indexes
            branch = branch_lookup[bp_branch[bp]]
            if haskey(branch, "rate_a")
                buspairs[bp]["rate_a"] = branch["rate_a"]
            end
            if haskey(branch, "c_rating_a")
                buspairs[bp]["c_rating_a"] = branch["c_rating_a"]
            end
        end
    
        return buspairs
    end
    for (nw, nw_ref) in ref[:nw]
        nw_ref[:all_bus] = Dict(j["bus_i"] => j for (_, j) in data["bus"])
        nw_ref[:all_branch] = Dict{Int, Any}(parse(Int, k) => v for (k,v) in data["branch"] if v["br_status"] != PM.pm_component_status_inactive["branch"])
        nw_ref[:all_buspairs] = calc_all_buspair_parameters(nw_ref[:all_bus], nw_ref[:all_branch], nw_ref[:conductor_ids], haskey(nw_ref, :conductors))
    end
end