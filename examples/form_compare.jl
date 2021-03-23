using Test
using PowerModels
using NetDecOPF
using DualDecomposition
using Ipopt, OSQP
# using BenchmarkTools
const DD = DualDecomposition

sub_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "warm_start_init_point" => "yes", "linear_solver" => "ma27")
optimizer = optimizer_with_attributes(OSQP.Optimizer, "warm_start" => true, "verbose" => false)

file = "/home/wzhang483/pglib-opf/api/pglib_opf_case89_pegase__api.m"
data = parse_file(file)
partitions = metis_cluster(data, 10)

function run_model(modeltype, partitions, data)
    dn_model = decompose(data, partitions, modeltype, NetDecOPF.build_acopf_with_free_lines)

    algo = init_DD_algo(dn_model)
    set_subnet_optimizer!(dn_model, sub_optimizer)

    # Change parameters (for example)
    params = BM.Parameters()
    BM.set_parameter(params, "maxiter", 1000)
    BM.set_parameter(params, "ϵ_s", 1.e-5)

    # Lagrange master method
    LM = DD.BundleMaster(BM.ProximalMethod, optimizer, params)

    DD.run!(algo, LM)
end

# @benchmark run_model(ACRPowerModel, partitions, data)
# @benchmark run_model(W_ACRModel, partitions, data)
run_model(ACRPowerModel, partitions, data) # without W variables
# run_model(W_ACRModel, partitions, data) # with W variables
