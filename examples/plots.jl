using Plots
using DelimitedFiles

# t = readdlm("output/$(casename)/subsolve_time.txt_0.txt", ',', header=true)[1]
# plot(t)

##
nprocs=32
casename = "pglib_opf_case300_ieee_p$(nprocs)"
subtime = Dict{Int,Vector{Float64}}()
for i=1:nprocs
    t = readdlm("output/$(casename)/subsolve_time.txt_$(i-1).txt", ',', header=true)[1]
    for j=1:1
        subtime[(j-1)*4+i] = t[:,j]
    end
end

##
plot(
    [v for (i,v) in subtime],
    legend=:none,
    xlabel="Iteration",
    ylabel="Subproblem time",
    title=casename
)

##
num_iters = length(subtime[1])
accum_time = zeros(num_iters)
for n=1:num_iters
    if n > 1
        accum_time[n] += accum_time[n-1]
    end
    accum_time[n] += maximum([v[n] for (i,v) in subtime])
end

##
plot(
    accum_time, 
    legend=:none, 
    xlabel="Iteration", 
    ylabel="Accumulated subproblem time",
    title=casename
)

##
subobj = readdlm("output/$(casename)/subobj_value.txt")[:] * -1.0
bestlb = deepcopy(subobj)
for i=2:length(bestlb)
    bestlb[i] = max(bestlb[i],bestlb[i-1])
end
plot(
    bestlb,
    xlabel="Iteration",
    ylabel="Lagrangian bound",
    title=casename,
    legend=:none,
)