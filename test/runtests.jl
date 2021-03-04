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
# partitions = [[2, 3, 4], [1, 5]]
partitions = metis_cluster(file, 2)
dn_model = decompose(data, partitions, ACRPowerModel, NetDecOPF.build_acopf_with_free_lines)
# dn_model = decompose(data, partitions, SOCBFPowerModel, NetDecOPF.build_socbf_with_free_lines)
algo = init_DD_algo(dn_model, method = BM.ProximalMethod, maxiter=1000)
# algo = init_DD_algo(dn_model, method = BM.TrustRegionMethod, maxiter=5000)
set_subnet_optimizer!(dn_model, sub_optimizer)
DD.run!(algo, optimizer)
@show DD.dual_objective_value(algo)
@show DD.dual_solution(algo)

@test true
