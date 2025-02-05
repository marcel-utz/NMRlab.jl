using NMRlab
using Test

@testset "NMRlab.jl FileIO" begin
    # Write your tests here.
   
    ex=NMRlab.Examples.Data["HCC cell culture media spectra"]["files"][1]
    params,d=NMRlab.load(ex,:Bruker)

    @test size(d) == size(d.coord[1])
    @test isapprox(2.616357072520625e11 + 1.335611666755e11im,sum(d))
end
