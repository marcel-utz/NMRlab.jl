module NMRlab

    export SpectData,coords,load
    export NMRProcessor, Chain, FourierTransform, Apodize, ZeroFill
    export SpinSim
    export GISSMO

    include("DataSet.jl")
    include("Examples.jl")
    include("NMRProcessor.jl")
    include("FileIO.jl")
    include("SpinSim.jl")
    include("GISSMO.jl")

end
