# NetDecOPF.jl

[![Run tests](https://github.com/Argonne-National-Laboratory/NetDecOPF.jl/actions/workflows/test.yml/badge.svg)](https://github.com/Argonne-National-Laboratory/NetDecOPF.jl/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/Argonne-National-Laboratory/NetDecOPF.jl/branch/master/graph/badge.svg?token=NBSOXPS06D)](https://codecov.io/gh/Argonne-National-Laboratory/NetDecOPF.jl)
[![DOI](https://zenodo.org/badge/285039510.svg)](https://zenodo.org/badge/latestdoi/285039510)

This package implements a network decomposition formulation to solve large-scale optimization problems over power networks, e.g. ACOPF, OTS etc.

## Installation

This package can be installed by

```julia
] add NetDecOPF
```

## Basic Usage

This package interfaces with power system data generated from [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl). Given the directory `file` to a power system data file
```julia
julia> using NecDecOPF, PowerModels

julia> file = "../examples/case5.m";

julia> data = parse_file(file);
```

A graph partition of the power system network is needed for the network decomposition. This package provides a function to generate partitions using METIS:
```julia
julia> partitions = metis_cluster(data, 2);
2-element Array{Array{Int64,1},1}:
 [2, 3, 4]
 [1, 5]
```

Alternatively, users can define custom partition to the network in the format of `Array` of `Array`s, where each `Array` collects all nodes of one partition:
```julia
julia> partitions = [[2, 3], [1, 4, 5]]
2-element Array{Array{Int64,1},1}:
 [2, 3, 4]
 [1, 5]
```

Now we are ready to call the `decompose` function to generated network decomposition:
```julia
julia> dn_model = decompose(data, partitions, ACRPowerModel, build_acopf_with_free_lines);
```
This code generates 2 subproblems based on the partitions specified by `partitions`. `ACRPowerModel` and `build_acopf_with_free_lines` specifies that the generated subproblems are rectangular ACOPF. This package also provides some relaxations for network decomposition:
```julia
dn_model = decompose(data, partitions, SDPWRMPowerModel, NetDecOPF.build_acopf_with_free_lines) # SDP relaxation for each subproblem
dn_model = decompose(data, partitions, SOCBFPowerModel, NetDecOPF.build_socbf_with_free_lines) # SOC relaxation of the branch flow model for each subproblem
dn_model = decompose(data, partitions, SOCWRPowerModel, NetDecOPF.build_socwr_with_free_lines) # SOC relaxation of the bus injection model for each subproblem
```

This package uses [DualDecomposition.jl](https://github.com/kibaekkim/DualDecomposition.jl) to solve the network decomposition problem, where the network decomposition is formulated as a dual decomposition problem and can be solved by methods such as bundle methods. This can be done as follows:
```julia
julia> using DualDecomposition; const DD = DualDecomposition

julia> sub_optimizer = optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "warm_start_init_point" => "yes")

julia> algo = init_DD_algo(dn_model); # initialize the dual decomposition algorithm

julia> set_subnet_optimizer!(dn_model, sub_optimizer); # set the optimizer for each network subproblem

julia> LM = DD.BundleMaster(BM.ProximalMethod, optimizer); # initialize the bundle method with optimizer
```
Here line 3 initializes the dual decomposition based on our decomposed network model `dn_model`. Line 4 sets the optimizer (`sub_optimizer`) for each network subproblem. `sub_optimizer` needs to be set as an optimizer capable of solving the subproblem. For instance, if each subproblem is formulated as ACOPF, then `sub_optimizer` has to be able to solve nonconvex QCQP. Line 5 sets the optimizer (`optimizer`) for the master problem in the bundle method.

In the end we solve the network decomposition by calling the `DD.run!` function:
```julia
julia> DD.run!(algo, LM)
```

## Citing this package

```
@misc{NetDecOPF.jl.0.1.0,
  author       = {Zhang, Weiqi and Kim, Kibaek},
  title        = {{NetDecOPF.jl: Implementation of network decomposition in Julia}},
  month        = Mar,
  year         = 2021,
  doi          = {10.5281/zenodo.4592258},
  version      = {0.1.0},
  publisher    = {Zenodo},
  url          = {https://doi.org/10.5281/zenodo.4592258}
}
```

## Acknowledgement

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.
