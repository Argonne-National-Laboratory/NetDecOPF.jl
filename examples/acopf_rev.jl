using PowerModels
using NetDecOPF
using DualDecomposition
using Ipopt, KNITRO
using Mosek, MosekTools
using JuMP
using BundleMethod
const DD = DualDecomposition
const PM = PowerModels
using JLD, HDF5
using COSMO
using OSQP

PowerModels.silence()
sub_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 1, "warm_start_init_point" => "yes", "linear_solver" => "ma27")
# sub_optimizer = optimizer_with_attributes(KNITRO.Optimizer, "outlev" => 0, "algorithm" => 4, "par_numthreads" => 1, "strat_warm_start" => 1)
# sub_optimizer = optimizer_with_attributes(KNITRO.Optimizer, "outlev" => 1, "strat_warm_start" => 1)
optimizer = optimizer_with_attributes(OSQP.Optimizer, "verbose" => false)

dir_to_pglib = "" # CHANGE THIS TO DIRECTORY TO PGLIB-OPF
file = dir_to_pglib * "pglib_opf_case5_pjm.m"
data = parse_file(file)

# i = 320
partitions = [[bus["bus_i"]] for (_,bus) in data["bus"]]
# partitions = metis_cluster(file, 2)

# This is decomposing with full W matrix as consensus
dn_model = decompose(data, partitions, W_global_ACRModel, NetDecOPF.build_acopf_with_free_lines, extra_ref_extensions=[NetDecOPF.ref_add_global_bus!])

# This is decomposing with vr and vi as consensus
dn_model = decompose(data, partitions, W_ACRModel_V, NetDecOPF.build_acopf_with_free_lines, extra_ref_extensions=[NetDecOPF.ref_add_global_bus!])

algo = init_DD_algo(dn_model)
set_subnet_optimizer!(dn_model, sub_optimizer)
LM = DD.BundleMaster(BM.ProximalMethod, optimizer)
DD.run!(algo, LM)