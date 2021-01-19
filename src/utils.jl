
function metis_cluster(file, N_partition)
    data = PM.parse_matpower(file)
    N = length(data["bus"])
    buses = collect(keys(data["bus"]))
    buses = sort([parse(Int, i) for i in buses])
    bus_idx_dict = Dict(buses[i] => i for i in eachindex(buses))
    bus_idx_dict_rev = Dict(i => buses[i] for i in eachindex(buses))
    L = length(data["branch"])
    lines = [(data["branch"]["$(i)"]["f_bus"], data["branch"]["$(i)"]["t_bus"]) for i in 1:L]

    g = SimpleGraph(N)
    for (i,j) in lines
        LightGraphs.add_edge!(g, bus_idx_dict[i], bus_idx_dict[j])
    end

    part = Metis.partition(g, N_partition)

    partitions = [Int64[] for i in 1:N_partition]

    for i in buses
        part_id = part[bus_idx_dict[i]]
        push!(partitions[part_id], i)
    end

    return partitions
end