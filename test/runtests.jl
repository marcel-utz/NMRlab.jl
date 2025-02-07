using NMRlab
using Test


@testset "SpectData basics" begin
    
    # construction
    d = NMRlab.SpectData(ones(45,45,45),(range(0.0,1.0,45),range(1.0,2.0,45),range(2.0,3.0,45)))
    
    # slicing
    @test typeof(d[:,1,1]) == NMRlab.SpectData{Float64,1}

end

@testset "NMRlab.jl FileIO" begin
    # Write your tests here.
   
    ex=NMRlab.Examples.Data["HCC cell culture media spectra"]["files"][1]
    params,d=NMRlab.load(ex,:Bruker)

    @test size(d) == size(d.coord[1])
    @test isapprox(2.616357072520625e11 + 1.335611666755e11im,sum(d))
end
