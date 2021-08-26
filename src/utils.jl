function metis_cluster(data::Dict, N_partition, alg=:KWAY)
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

    part = Metis.partition(g, N_partition, alg=alg)

    metis_options = copy(Metis.options)
    metis_options[Metis.METIS_OPTION_CONTIG] = 1

    # The following code is copied from Metis library
    G = Metis.graph(g)
    part = Vector{Metis.idx_t}(undef, G.nvtxs)
    vwgt = isdefined(G, :vwgt) ? G.vwgt : C_NULL
    edgecut = fill(Metis.idx_t(0), 1)
    if alg === :RECURSIVE
        Metis.METIS_PartGraphRecursive(G.nvtxs, Metis.idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, C_NULL,
                                 Metis.idx_t(N_partition), C_NULL, C_NULL, metis_options, edgecut, part)
    elseif alg === :KWAY
        Metis.METIS_PartGraphKway(G.nvtxs, Metis.idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, C_NULL,
                            Metis.idx_t(N_partition), C_NULL, C_NULL, metis_options, edgecut, part)
    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end

    partitions = [Int64[] for i in 1:N_partition]

    for i in buses
        part_id = part[bus_idx_dict[i]]
        push!(partitions[part_id], i)
    end

    return partitions
end

metis_cluster(file::String, N_partition, alg=:KWAY) = metis_cluster(PM.parse_matpower(file), N_partition, alg)