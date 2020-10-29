mutable struct NetDecModel{T <: PM.AbstractPowerModel}
    node_partition::Vector{Set{Int64}}
    model_list::Vector{T}
    split_line_vars::Dict{Int64, Dict}
    function NetDecModel(
        node_partition::Vector{Set{Int64}},
        model_list::Vector{T} = T[],
        split_line_vars::Dict{Int64, Dict} = Dict{Int64, Dict}()
        ) where T <: PM.AbstractPowerModel
        return new{T}(node_partition, model_list, split_line_vars)
    end
end

get_partition(dm_model::NetDecModel)::Vector{Set{Int64}} = dm_model.node_partition

get_models(dm_model::NetDecModel)::Vector{<:PM.AbstractPowerModel} = dm_model.model_list

get_split_vars(dm_model::NetDecModel)::Dict{Int64, Dict} = dm_model.split_line_vars

function set_subnet_optimizer!(dm_model::NetDecModel, optimizer)
    for pm in get_models(dm_model)
        JuMP.set_optimizer(pm.model, optimizer)
    end
    return
end