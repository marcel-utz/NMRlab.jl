"""
module **GISSMO**

Interface routines with GISSMO database.

Hesam Dashti, Jonathan R. Wedell, William M. Westler, Marco Tonelli, David Aceti, Gaya K. Amarasinghe, John L. Markley, and Hamid R. Eghbalnia, Applications of Parametrized NMR Spin Systems of Small Molecules, *Anal. Chem.*, 2018, **90** (18), pp 10646–10649, DOI: 10.1021/acs.analchem.8b02660
Hesam Dashti, William M. Westler, Marco Tonelli, Jonathan R. Wedell, John L. Markley, and Hamid R. Eghbalnia, Spin System Modeling of Nuclear Magnetic Resonance Spectra for Applications in Metabolomics and Small Molecule Screening, *Analytical Chemistry*, 2017, **89** (22), pp 12201–12208, doi: 10.1021/acs.analchem.7b02884
""" 
module GISSMO

using HTTP
using LightXML
import JSON
using  NMRlab.SpinSim


export Hamiltonian, search


"""
    function Hamiltonian(fn::String;freq=600.0,ctr=4.8)
    function Hamiltonian(system::XMLElement;freq=600.0,ctr=4.8)

computes the Hamiltonian of an entry in the GISSMO database and returns it in the natural basis as a sparse matrix.
`fn` is the GISSMO reference name of the compound. Keyword parameters are used to indicate spectrometer
base frequency (in MHz) and the spectral zero point (carrier position) in ppm. 
`maxSpin` represents the maximum size of the spin system that will be simulated.
If the spin system size exceeds `maxSpin`, the function returns `(-1,nothing)`.
"""
function Hamiltonian(id::String;freq=600.0,ctr=4.8)
    xs=String(HTTP.request("GET","https://gissmo.bmrb.io/entry/$(id)/simulation_1/spin_simulation.xml").body)
    xdoc=parse_string(xs)
    nspin,H = Hamiltonian(root(xdoc);freq=freq,ctr=ctr)
end

function Hamiltonian(system::XMLElement;freq=600.0,ctr=4.8,maxSpin=25)
    compound = system["name"][1]
    xspin = system["coupling_matrix"][1]
    chem_shifts = Array{Float64,1}([])
    xcs = xspin["chemical_shifts_ppm"][1]["cs"]
    chem_shifts = [parse(Float64,attribute(c,"ppm")) for c in xcs]
    nspin=length(chem_shifts)
    
    if nspin>maxSpin
        return (-1,nothing)
    end
    
    H=sum(j->2pi*freq*(chem_shifts[j]-ctr)*SpinOp(nspin,Sz,j),1:nspin)
    
    xJs=xspin["couplings_Hz"][1]["coupling"]
    for c in xJs
        k=parse(Int64,attribute(c,"from_index"))
        l=parse(Int64,attribute(c,"to_index"))
        J=parse(Float64,attribute(c,"value"))
        
        H.+=2pi*J*OpJstrong(nspin,k,l)
    end
    nspin,H
end

"""
    function search(term::String)

searches the online GISSMO database for `term` and returns a `JSON` object
with the search result.
"""
function search(term::String)
    res=String(HTTP.request("GET","https://gissmo.bmrb.io/search?term=$(term)").body)
    return JSON.parse(res)
end

end