module NetDecOPF

# Import necessary packages
import JuMP
import PowerModels
import DualDecomposition
import BundleMethod
import LinearAlgebra: I
const PM = PowerModels
const DD = DualDecomposition
const BM = BundleMethod

using Metis
using LightGraphs, SimpleWeightedGraphs

include("model.jl")
include("data.jl")
include("prob.jl")
include("solution.jl")
include("utils.jl")

export NetDecModel, get_partition, get_models, get_split_vars, set_subnet_optimizer!

export find_neighbor_buses, generate_subnet_data, ref_add_cut_bus!, ref_add_cut_branch!

export decompose, build_acopf_with_free_lines, build_socbf_with_free_lines, build_acots_with_free_lines, collect_split_vars

export init_DD_algo, add_split_vars_to_algo!

export metis_cluster

include("PM/ACR/Large_W_ACRModel.jl")
include("PM/WR/W_ACRModel.jl")
include("PM/WR/WU_ACRModel.jl")
include("PM/WR/W_global_ACRModel.jl")
include("PM/WR/W_ACRModel_V.jl")
include("PM/WR/NonconvexWRModel.jl")
include("PM/WRM/SparseSDPWRMPowerModel.jl")
include("PM/WRM/RelaxedSDPWRMPowerModel.jl")

export W_ACRModel, WU_ACRModel, Large_W_ACRModel, W_global_ACRModel, W_ACRModel_V, NonconvexACRModel, RelaxedSDPWRMPowerModel

end # module
