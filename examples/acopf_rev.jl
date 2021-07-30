using Revise

using PowerModels
using NetDecOPF
using DualDecomposition
using BundleMethod

const DD = DualDecomposition

using Ipopt
# using KNITRO
using filterSQP
using AmplNLWriter
# using JuMP
using OSQP

PowerModels.silence()
sub_optimizer = optimizer_with_attributes(filterSQP.Optimizer, "use_warm_start" => false, "eps" => 1e-4)
# sub_optimizer = optimizer_with_attributes(
#     () -> AmplNLWriter.Optimizer("/Users/kibaekkim/Documents/REPOS/filterSQP/build/bin/amplfilter"), 
#     "iprint" => 1,
#     "eps" => 1e-4)
# sub_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 1, "warm_start_init_point" => "yes", "linear_solver" => "ma27")
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
# dn_model = decompose(data, partitions, W_ACRModel_V, NetDecOPF.build_acopf_with_free_lines, extra_ref_extensions=[NetDecOPF.ref_add_global_bus!])

# initialization
for pm in dn_model.model_list
    m = pm.model
    for v in JuMP.all_variables(m)
        if JuMP.has_lower_bound(v) && JuMP.has_upper_bound(v)
            JuMP.set_start_value(v, (JuMP.lower_bound(v) + JuMP.upper_bound(v)) / 2)
        end
    end
end

algo = init_DD_algo(dn_model)
set_subnet_optimizer!(dn_model, sub_optimizer)
LM = DD.BundleMaster(BM.ProximalMethod, optimizer)
DD.run!(algo, LM)