using MISPExamples

@static if VERSION >= v"0.7-"
    using Test
else
    using Base.Test
    import Compat: @info
end

const TSE = MISPExamples.TikhonovSamarskii
@testset "TikhonovSamarskii Example" begin
    # Simple test that _bisect works as expected.
    approxRoot, errorEstimate = @inferred TSE._bisect(x->x^2-2, 0.0, -2.0, 2.0, 2.0, 30)
    @test isapprox(approxRoot, sqrt(2), rtol = abs(errorEstimate))

    @testset "eltype T=$T, c2=$c2, k2=$(k2ratio)" for T in [Float32, Float64, BigFloat], c2 in [0, 1, 2], k2ratio in [1, 2]
        t1 = TSE.ProblemData(-one(T), T(c2), one(T), k2ratio*one(T))

        # Test that alphaResidual function is inferrable
        @inferred TSE.alphaResidual(1.0, t1)
        @inferred TSE.alphaResidual(1, t1)

        # Test that we can find alpha
        t1WithAlpha = TSE.withAlpha(t1, 1e-5, 10)
        @test TSE.hasAlpha(t1WithAlpha)
        alphaValue = @inferred TSE.alpha(t1WithAlpha)
        @test alphaValue >= 0

        # For a SimpleProblem (simplification applied as in Tikhonov and Samarskii) A2=B2=0.
        @test !isapprox(c2, 0) || isapprox(TSE.A2(t1WithAlpha), 0)
        @test !isapprox(c2, 0) || isapprox(TSE.B2(t1WithAlpha), 0)

        # Test that parameter functions are all inferrable.
        xi1 = @inferred TSE.xi(1.0, t1WithAlpha)
        @test isapprox(xi1, alphaValue)
        u1v = @inferred TSE.u1(0.1, 0.1, t1WithAlpha)
        uv = @inferred TSE.u(0.1, 0.1, t1WithAlpha)
        @test isapprox(u1v, uv)
        u2v = @inferred TSE.u2(xi1, 0.1, t1WithAlpha)
        @test u2v >= 0
        gv = @inferred TSE.g(1.0, t1WithAlpha)
        @test gv >= 0
        pv = @inferred TSE.p(1.0, t1WithAlpha)
        @test pv >= 0
        nuv = @inferred TSE.nu(1.0, t1WithAlpha)
        @test nuv >= 0


        smallAccum = 0.0
        for t in TSE.tGrid(5)[2:end]
            xit = TSE.xi(t, t1WithAlpha)
            smallAccum += TSE.u1(xit, t, t1WithAlpha) - TSE.u2(xit, t, t1WithAlpha)
        end
        @test smallAccum < 1e-6
    end
end

# Tests for "Default" Tikhonov & Samarskii Example
const TSD = MISPExamples.TikhonovSamarskiiDefault
@testset "TikhonovSamarskiiDefault" begin
    # Test that parameter functions are all inferrable.
    xi1 = @inferred TSD.xi(1.0)
    @test xi1 > TSD.t1.alphaEstimate # since tShift > 0

    vv = @inferred TSD.v(0, 0.1)
    vv2 = @inferred TSD.v(TSD.spaceLength, 0.1)
    @test vv * vv2 < 0 # Change in phase across domain

    gv = @inferred TSD.g(1.0)
    @test gv >= 0

    pv = @inferred TSD.p(1.0)
    @test pv >= 0

    nuv = @inferred TSD.nu(1.0)
    @test nuv >= 0

    PhiLeft = @inferred TSD.Phi(0)
    PhiRight = @inferred TSD.Phi(TSD.spaceLength)
    @test PhiLeft * PhiRight < 0 # Change in phase across domain initially

    GammaV = @inferred TSD.Gamma(0)
    @test GammaV > 0

    bv = @inferred TSD.b(vv)
    bv2 = @inferred TSD.b(vv2)
    @test bv * bv2 < 0 # Change in phase across domain results in change in sign for b

    fv = @inferred TSD.f(0, 0)
    @test isapprox(fv, 0) # This example only

    tgv = TSD.tGrid(6)
    @test length(tgv) == 6
    @test first(tgv) == 0
    @test isapprox(last(tgv), TSD.finalMoment)

    xgv = TSD.xGrid(7)
    @test length(xgv) == 7
    @test first(xgv) == 0
    @test isapprox(last(xgv), TSD.spaceLength)
end

@static if VERSION >= v"0.7-"
    @testset "exampleFiles" begin
        expectedNames = [:bcoeff, :f, :Phi, :g, :p, :Gamma]
        @testset "Example module $M satisfies correct interface" for M in filter(!isequal(:MISPExamples), names(MISPExamples))
            MNames = names(getproperty(MISPExamples, M))
            @test all(n in MNames for n in expectedNames)
        end
    end
end

@static if VERSION >= v"0.7-"

@testset "exampleFiles" begin
    expectedNames = [:bcoeff, :f, :Phi, :g, :p, :Gamma]
    @testset "Example module $M satisfies correct interface" for M in filter(!isequal(:MISPExamples), names(MISPExamples))
        MNames = names(getproperty(MISPExamples, M))
        @test all(n in MNames for n in expectedNames)
    end
end

end
