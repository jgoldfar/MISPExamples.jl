module TikhonovSamarskiiExample
if VERSION >= v"0.6-"
    import SpecialFunctions: erf
end

import Base.FastMath: sqrt_fast

#TODO: Update alphaEstimate inplace using a Ref type
# on first access.
struct ProblemData{T}
    c1::T
    c2::T
    k1::T
    k2::T
    gammaVal::T
    alphaEstimate::Float64
    function ProblemData(c1::T, c2::T, k1::T, k2::T, gammaVal::T) where {T<:Real}
        @assert k1>0 && k2>0 && c1<0
        new{T}(c1, c2, k1, k2, gammaVal, NaN)
    end
    function ProblemData(c1::T, c2::T, k1::T, k2::T, gammaVal::T, alphaEstimate) where {T<:Real}
        @assert k1>0 && k2>0 && c1<0
        new{T}(c1, c2, k1, k2, gammaVal, T(alphaEstimate))
    end
end
ProblemData(c1::Real, c2=zero(c1), k1=one(c2), k2 = k1+1, gammaVal = one(c1)) = ProblemData(promote(c1, c2, k1, k2, gammaVal)...)
ProblemData(;c1::Real=-1.0, c2=zero(c1), k1=one(c2), k2 = k1+1, gammaVal = one(c1)) = ProblemData(c1, c2, k1, k2, gammaVal)

alpha(d::ProblemData) = d.alphaEstimate
hasAlpha(d::ProblemData) = !isnan(alpha(d))

function Base.show(io::IO, d::ProblemData)
    if get(io, :compact, false)
        print(io, "ProblemData(", d.c1,
        ", ", d.c2,
        ", ", d.k1,
        ", ", d.k2,
        ", ", d.gammaVal,
        ", ", (hasAlpha(d) ? d.alphaEstimate : "unknown"),
        ")")
    else
        print(io, "ProblemData(c1=", d.c1,
        ", c2=", d.c2,
        ", k1=", d.k1,
        ", k2=", d.k2,
        ", gamma=", d.gammaVal,
        ", alpha=", (hasAlpha(d) ? d.alphaEstimate : "unknown"),
        ")")
    end
end

B1(d::ProblemData) = hasAlpha(d) ? -(d.c1) / erf(alpha(d)/(2*sqrt(d.k1))) : NaN
B2(d::ProblemData) = hasAlpha(d) ? (d.c2) / (1-erf(alpha(d)/(2*sqrt(d.k2)))) : NaN
A1(d::ProblemData) = (d.c1)
A2(d::ProblemData) = hasAlpha(d) ? -erf(alpha(d)/(2*sqrt(d.k2)))*B2(d) : NaN

function alphaResidual(alphaTest, d::ProblemData)
    k1s = sqrt_fast(d.k1)
    k2s = sqrt_fast(d.k2)
    alphaOverTwo = alphaTest / 2

    V1 = k1s * (d.c1) * exp(-(alphaOverTwo/k1s)^2)/erf(alphaOverTwo/k1s)
    V2 = k2s * (d.c2) * exp(-(alphaOverTwo/k2s)^2)/(1-erf(alphaOverTwo/k2s))
    V3 = (d.gammaVal)*sqrt(pi)*alphaOverTwo
    return V1 + V2 + V3
end

function _bisect(f, a::Real, fa, b::Real, fb, numits::Int = 100)
    if isapprox(fb, 0)
        return b, fb
    elseif  isapprox(fa, 0) || (fa * fb > 0)
        return a, fa
    end

    c = (a+b)/2
    fc = f(c)
    if numits < 1 || isapprox(fc, 0)
        return c, fc
    end

    if fa * fc < 0 # Sign change on [a,c]
        return _bisect(f, a, fa, c, fc, numits-1)
    elseif fb * fc < 0 # Sign change on [c,b]
        return _bisect(f, c, fc, b, fb, numits-1)
    else # No sign change (how did we get here?)
        return c, fc
    end
end

function withAlpha(d::ProblemData, a, b, tolerance::Real = 1e-6)
    residual(alp) = alphaResidual(alp, d)
    ra = residual(a)
    rb = residual(b)

    alphaEstimate, errorEstimate = _bisect(residual, a, ra, b, rb)

    if errorEstimate > tolerance
        error("No alpha found with given tolerance and initial interval.")
        return d
    end

    ProblemData(d.c1, d.c2, d.k1, d.k2, d.gammaVal, alphaEstimate)
end

"""
    finalMoment

Last time we will use in simulation.
"""
const finalMoment = 1

## Actual model functions.
xi(t, d::ProblemData) = hasAlpha(d) ? alpha(d) * sqrt(t) : NaN

"""
    lBdy(d)

Left boundary point chosen as twice the magnitude of the "known" final position for problem with data `d` at time `finalMoment`.
"""
lBdy(d::ProblemData) = 2*xi(finalMoment, d)

## Functions below all assume hasAlpha(data)
u1(x, t, d::ProblemData) = d.c1 + B1(d)*erf(x/(2 * sqrt_fast(t * d.k1)))
u2(x, t, d::ProblemData) = A2(d) + B2(d)*erf(x/(2 * sqrt_fast(t * d.k2)))

"""
    u(x,t,d)

Analytic solution to untransformed problem with data `d` at `(x,t)`.
"""
function u(x, t, d::ProblemData{T}) where {T}
    if isapprox(t, 0)
        if x < 0
            return first(promote(d.c1, zero(x), zero(t)))
        else
            return first(promote(d.c2, zero(x), zero(t)))
        end
    end
    xiVal = xi(t, d)
    if x < xiVal
        return u1(x, t, d)
    elseif x > xiVal
        return u2(x, t, d)
    end
    # Fallthrough definition
    return first(promote(zero(x), zero(t), zero(T)))
end


# Note: Minimal interface is exported below.
export g, p, f, Gamma, bcoeff, Phi

"""
    g(t,d)

Analytic heat flux for solution to transformed or untransformed problem with data `d` at `(0,t)`.
"""
g(t, d::ProblemData) = (sqrt_fast(d.k1)*B1(d)) / sqrt_fast(pi * t)

"""
    gAverageValue(d)

Average value of function `g(t)` when problem data is `d`.
"""
gAverageValue(d::ProblemData) = (sqrt_fast(d.k1)*9*B1(d)) / (5*sqrt_fast(pi))

"""
    p(t,d)

Analytic heat flux for solution to transformed or untransformed problem with data `d` at `(lBdy(d),t)`.
"""
p(t, d::ProblemData) = (sqrt_fast(d.k1) * B2(d) / sqrt_fast(pi * t)) * exp(-(alpha(d)^2)/(d.k2*t))

"""
    nu(t,d)

Analytic temperature function for solution to untransformed problem with data `d` at `(lBdy(d),t)`.
"""
nu(t, d::ProblemData) = A2(d) + B2(d)*erf(alpha(d)/sqrt_fast(d.k2*t))

"""
    F(u,d)

Kirchoff transformation for PDE problem with data `d` at value `u`.

The output values should be considered valid "`v`" values.
"""
function F(u, d::ProblemData)
    if u>0
        d.k2*u
    else
        d.k1*u
    end
end
## Below interface is all that examples _must_ define

"""
    Gamma(t,d)

Analytic transformed temperature of solution to transformed problem with data `d` at `(lBdy(d),t)`.
"""
Gamma(t, d::ProblemData) = F(nu(t, d), d)

"""
    bcoeff(x,t,v,d)

Function appearing in parabolic problem after transformation,

    (bcoeff(v))_t - v_{xx} = f

when problem data is given according to Tikhonov & Samarskii, pp. 283--288, with data `d` at `(x,t,v)`.
"""
function bcoeff(x,t,v,d::ProblemData)
    if v>0
        d.gammaVal + v/d.k2
    else
        v/d.k1
    end
end

"""
    g(t,d)

Analytic heat flux for solution to transformed or untransformed problem with data `d` at `(0,t)`.
"""
Phi(x, d::ProblemData) = F(u(x, 0, d), d)

"""
    f(x,t)

Analytic heat sources for solution to transformed or untransformed problem with data `d` at `(x,t)`.
"""
f(x,t) = first(promote(zero(x), zero(t)))

"""
    v(x,t,d)

Analytic solution to transformed problem with data `d` at `(x,t)`.
"""
v(x, t, d::ProblemData) = F(u(x, t, d), d)

end
