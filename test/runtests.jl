using Test
using PowerModels
using DeNet
using DualDecomposition
using Ipopt
const DD = DualDecomposition

sub_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

file = "../examples/case5.m"
data = parse_file(file)
partitions = [[2, 3, 4], [1, 5]]
dn_model = decompose(data, partitions, ACRPowerModel, DeNet.build_acopf_with_free_lines)
# dn_model = decompose(data, partitions, SOCBFPowerModel, DeNet.build_socbf_with_free_lines)
algo = init_DD_algo(dn_model)
set_subnet_optimizer!(dn_model, sub_optimizer)
DD.run!(algo, optimizer)