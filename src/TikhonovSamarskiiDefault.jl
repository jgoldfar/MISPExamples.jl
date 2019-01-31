"""
TikhonovSamarskiiDefault: Provides a example problem for MISP with some "default" values chosen for
various parameters.
"""
module TikhonovSamarskiiDefault
using ..TikhonovSamarskii

# Note: Minimal interface is exported below.
export g, p, f, Gamma, bcoeff, Phi

# The Tikhonov & Samarskii example is parameterized by the coefficients appearing in
# the various equations.
function initProblem(c1=-1.0, c2=1.0, k1=1.0, k2=1.0, gammaVal=1.0)
    # The module provides a method to find the last coefficient, alpha, approximately using
    # bisection method; the next two lines solve for and extract this value of alpha.
    TikhonovSamarskii.withAlpha(
        TikhonovSamarskii.ProblemData(
            c1, #c1
            c2, #c2
            k1, #k1
            k2, #k2
            gammaVal), # Gamma (latent heat)
        1e-2, # minimum value for alpha search
        10)# maximum value for alpha search
end

const t1 = initProblem()

# Length of time domain
const finalMoment = TikhonovSamarskii.finalMoment

# Length of space domain
const spaceLength = TikhonovSamarskii.lBdy(t1)

# tShift is the initial moment we chose for this simulation.
const tShift = 1e-1

# Analytic solution v(x,t)
v(x,t) = TikhonovSamarskii.v(x, t + tShift, t1)

# Our interface does not support shifting the initial data in time, otherwise we would use the defined function Phi.
Phi(x) = v(x, 0) # Initial data
p(t) = TikhonovSamarskii.p(t + tShift, t1) # Flux on RHS
g(t) = TikhonovSamarskii.g(t + tShift, t1) # Flux on LHS
Gamma(t) = TikhonovSamarskii.Gamma(t + tShift, t1) # Temperature on RHS
b(v) = TikhonovSamarskii.bcoeff(0, 0, v, t1) # b(x,t,v)
f(x,t) = TikhonovSamarskii.f(x,t) # Heat sources
const gAverageValue = TikhonovSamarskii.gAverageValue(t1)
nu(t) = TikhonovSamarskii.nu(t + tShift, t1)

tGrid(NT) = TikhonovSamarskii.tGrid(NT)
xGrid(NX) = TikhonovSamarskii.xGrid(NX, t1)
end
