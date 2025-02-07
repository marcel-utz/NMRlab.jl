module NMRlab

    export SpectData,coords,load
    export NMRProcessor, Chain, FourierTransform, Apodize, ZeroFill

    include("DataSet.jl")
    include("Examples.jl")
    include("NMRProcessor.jl")
    include("FileIO.jl")

end
