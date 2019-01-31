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

        # TODO: add test that u1 and u2 are continuous by checking values at xi(t)
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
