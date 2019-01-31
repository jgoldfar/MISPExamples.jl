# README

![Bitbucket Pipeline Status](https://img.shields.io/bitbucket/pipelines/jgoldfar/mispexamples.jl.svg?style=flat)
[![Travis Build Status](https://travis-ci.org/jgoldfar/MISPExamples.jl.svg?branch=master)](https://travis-ci.org/jgoldfar/MISPExamples.jl)
[![codecov](https://codecov.io/gh/jgoldfar/MISPExamples.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jgoldfar/MISPExamples.jl)
[![Coverage Status](https://coveralls.io/repos/github/jgoldfar/MISPExamples.jl/badge.svg?branch=master)](https://coveralls.io/github/jgoldfar/MISPExamples.jl?branch=master)

`MISPExamples` defines a "general" version of the self-similar solution to the Stefan problem from Tikhonov & Samarskii (with all parameters free to choose) as well as a "default" version with pre-set parameters.

The README is generated automatically from the notebook [doc/TikhonovSamarskiiPlots.ipynb](/doc/TikhonovSamarskiiPlots.ipynb), so the commands below will work without modification from `/doc/`

It has no significant external dependencies except [`SpecialFunctions.jl`](https://github.com/JuliaMath/SpecialFunctions.jl); the first time you run the package, you'll need to `instantiate` that dependency:

```julia
using Pkg
Pkg.activate(@__DIR__)

# # If you've already run the command below, all you need to do is execute
Pkg.instantiate()

# # The first time you use this notebook, you'll have to run
# Pkg.develop(PackageSpec(url=joinpath(@__DIR__, "..")))
```

In most cases, the values of the parameters before and after applying the Kirchoff (integral) transformation are available, and the names of the functions available in `MISPExamples` should match that in e.g. Abdulla & Poggi (2019)


```julia
using MISPExamples
const TSE = MISPExamples.TikhonovSamarskiiDefault

@static if VERSION < v"0.7-"
    using Compat: stdout, range
end
```

If you'd like to examine the problem data, it is available as `t1` in that module.
For the most part, we don't need to work with these parameters ourselves; the main thing this is used for is to hold the solution of the rootfinding method needed to calculate `alpha`:


```julia
print(TSE.t1) # Verbose output by default
```

    ProblemData(c1=-1.0, c2=1.0, k1=1.0, k2=1.0, gamma=1.0, alpha=0.7555195764068221)

We should be able to plot the corresponding data using our preferred plotting package.


```julia
using PyPlot
```

The length of the space domain is "free" since the example in Tikhonov & Samarskii is a solution in the half-space, but we fix a right-boundary location to be a multiple of the final location of the "true" boundary to avoid degenerate phases appearing in our solution.


```julia
lBdy = TSE.spaceLength
xGrid = range(0, stop=lBdy, length=100)
tGrid = TSE.tGrid(100)
```




    0.0:0.010101010101010102:1.0



The transformed initial data is available as `Phi`; evidently the data is continuous, but note the lack of differentiability at the level `Phi=0`: this is the initial interface location


```julia
withfig(figure()) do
    plot(xGrid, map(TSE.Phi, xGrid))
    title(L"\Phi(x)")
    xlabel("x")
end
```




![png](/doc/figures/README_files/README_10_0.png?raw=true)




```julia
outputDir = joinpath(pwd(), "figures")
mkpath(outputDir)
```

The function `v(x,t)` is the analytic (true) transformed state vector


```julia
withfig(figure()) do
    surf(xGrid, tGrid, [TSE.v(x,t) for x in xGrid, t in tGrid])
    title(L"v(x,t)")
    xlabel("t")
    ylabel("x")
    # # Uncomment the line below to save this plot to `outputDir`
    # savefig(joinpath(outputDir, "TikhonovSamarskiiSolveV.png"))
end
```




![png](/doc/figures/README_files/README_13_0.png?raw=true)



We could also plot the solution as a series of slices for more convenient inclusion into e.g. a report.


```julia
withfig(figure()) do
    plot(xGrid, map(TSE.Phi, xGrid))
    for t in TSE.tGrid(4)
        plot(xGrid, map(x->TSE.v(x,t), xGrid))
    end
    legend(append!(["\\Phi"], ["t=$t" for t in TSE.tGrid(4)]))
    title(L"v(x,t)")
    xlabel("x")
    # # Uncomment the line below to save this plot to `outputDir`
    # savefig(joinpath(outputDir, "TikhonovSamarskiiSolveVSlice.png"))
end
```




![png](/doc/figures/README_files/README_15_0.png)



The functions `p` and `g` represent the heat flux on the right- and left-hand sides of the fixed domain, respectively.


```julia
withfig(figure()) do
    plot(tGrid, map(TSE.p, tGrid))
    title(L"p(t)")
    xlabel("t")
    legend(L"p(t) = k_2 u_x(\ell, t)", loc=8, fontsize="x-small")
    # # Uncomment the line below to save this plot to `outputDir`
    # savefig(joinpath(outputDir, "TikhonovSamarskiiSolveP.png"))
end
```




![png](/doc/figures/README_files/README_17_0.png)




```julia
withfig(figure()) do
    plot(tGrid, map(TSE.g, tGrid))
    title(L"g(t)")
    xlabel("t")
    legend("g(t) = k_1 u_x(0, t)", loc=8, fontsize="x-small")
    # # Uncomment the line below to save this plot to `outputDir`
    # savefig(joinpath(outputDir, "TikhonovSamarskiiSolveG.png"))
end
```




![png](/doc/figures/README_files/README_18_0.png?raw=true)



`nu(t)` is the Dirichlet data on the right-hand side of the domain


```julia
withfig(figure()) do
    plot(tGrid, map(TSE.nu, tGrid))
    legend("nu(t) = u(\ell, t)", loc=8, fontsize="x-small")
end
```




![png](/doc/figures/README_files/README_20_0.png?raw=true)



`Gamma` is the transformed Dirichlet data on the right-hand side of the domain.


```julia
withfig(figure()) do
    plot(tGrid, map(TSE.Gamma, tGrid))
    title(L"\Gamma(t)")
    xlabel("t")
    legend(L"\Gamma(t) = F(u(\ell, t))", loc=8, fontsize="x-small")
    # # Uncomment the line below to save this plot to `outputDir`
    # savefig(joinpath(outputDir, "TikhonovSamarskiiSolveGamma.png"))
end
```




![png](/doc/figures/README_files/README_22_0.png?raw=true)



Having `b(v)` and the function `v(x,t)`, we can show the composition, which should have a jump along the boundary curve:


```julia
withfig(figure()) do
    surf(xGrid, tGrid, [TSE.b(TSE.v(x,t)) for x in xGrid, t in tGrid])
    title(L"b(v(x,t))")
    xlabel("t")
    ylabel("x")
end
```




![png](/doc/figures/README_files/README_24_0.png?raw=true)



Having the analytic data, could (if we were so inclined) calculate the "optimal" error possible from a finite difference scheme; this in particular would serve as a check against any algebra issues.


```julia
using LinearAlgebra
NVals = 8:10:1000
L2ErrorVals = zeros(length(NVals))
for (Ni, N) in enumerate(NVals)

    lBdy = TSE.spaceLength
    xGrid = range(0, stop=lBdy, length=N)
    tGrid = TSE.tGrid(N)

    vapprox = [TSE.v(x,t) for x in xGrid, t in tGrid]
    @static if VERSION >= v"0.7-"
        btapprox = diff([TSE.b(vapprox[i, j]) for i in eachindex(xGrid), j in eachindex(tGrid)], dims=2) / step(tGrid)
        vxxapprox = diff(diff(vapprox, dims=1), dims=1)/(step(xGrid)^2)
    else
        btapprox = diff([TSE.b(vapprox[i, j]) for i in eachindex(xGrid), j in eachindex(tGrid)], 2) / step(tGrid)
        vxxapprox = diff(diff(vapprox, 1), 1)/(step(xGrid)^2)
    end
    L2ErrorVals[Ni] = vecnorm(btapprox[1:(end-2), :] - vxxapprox[:, 1:(end-1)])*step(xGrid)*step(tGrid)
end

withfig(figure()) do
    plot(NVals, L2ErrorVals)
    title("Optimal Error in Finite Difference Scheme")
    xlabel(L"N_x \equiv N_t")
    ylabel(L"\Vert b(v)_\bar{t} - v_{\bar{x}x}\Vert_{\ell_2, h, \tau}")
    # # Uncomment the line below to save this plot to `outputDir`
    # savefig(joinpath(outputDir, "TikhonovSamarskiiSolveApproxError.png"))
end
```




![png](/doc/figures/README_files/README_27_0.png?raw=true)



This example shows the spatial distribution of such a finite difference approximation; we can see that the errors all appear near the interface curve, as one might expect.


```julia
N = 64

lBdy = TSE.spaceLength
xGrid = range(0, stop=lBdy, length=N)
tGrid = TSE.tGrid(N)
vapprox = [TSE.v(x,t) for x in xGrid, t in tGrid]
@static if VERSION >= v"0.7-"
    btapprox = diff([TSE.b(vapprox[i, j]) for i in eachindex(xGrid), j in eachindex(tGrid)], dims=2) / step(tGrid)
    vxxapprox = diff(diff(vapprox, dims=1), dims=1)/(step(xGrid)^2)
else
    btapprox = diff([TSE.b(vapprox[i, j]) for i in eachindex(xGrid), j in eachindex(tGrid)], 2) / step(tGrid)
    vxxapprox = diff(diff(vapprox, 1), 1)/(step(xGrid)^2)
end
fapprox = btapprox[1:(end-2), :] - vxxapprox[:, 1:(end-1)]
withfig(figure()) do
    surf(tGrid[1:(end-1)], xGrid[1:(end-2)], fapprox)
    title(L"b(v)_\bar{t} - v_{\bar{x}x} \equiv f(x,t) \approx 0")
    xlabel("t")
    ylabel("x")
end
```




![png](/doc/figures/README_files/README_29_0.png?raw=true)


