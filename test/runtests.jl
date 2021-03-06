using Test
using PowerModels
using NetDecOPF
using DualDecomposition
using Ipopt
const DD = DualDecomposition

sub_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "warm_start_init_point" => "yes")
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "warm_start_init_point" => "yes")

file = "../examples/case5.m"
data = parse_file(file)

# Partition network
partitions = metis_cluster(file, 2)

dn_model = decompose(data, partitions, ACRPowerModel, NetDecOPF.build_acopf_with_free_lines)
# dn_model = decompose(data, partitions, SOCBFPowerModel, NetDecOPF.build_socbf_with_free_lines)

algo = init_DD_algo(dn_model)
set_subnet_optimizer!(dn_model, sub_optimizer)

# Change parameters (for example)
params = BM.Parameters()
BM.set_parameter(params, "maxiter", 3000)
BM.set_parameter(params, "ϵ_s", 1.e-5)

# Lagrange master method
LM = DD.BundleMaster(BM.ProximalMethod, optimizer, params)

DD.run!(algo, LM)

@test isapprox(DD.dual_objective_value(algo), 17751, rtol = 0.1)
# @show DD.dual_solution(algo)
