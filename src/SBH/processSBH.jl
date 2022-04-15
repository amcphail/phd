using Statistics

using FileIO, JLD2

using Unzip

include("StimulatedBoseHubbard.jl")

data_dir = "data/SBH"

fn = "results"

function processNumber(fn,ann)

    tout = []
    number = []
    phase = []
    number_var = []
    phase_var = []

    for i in 1:10
        o_fn = string(data_dir,"/",fn,ann,"-",i,".jld2")
        jldopen(o_fn,"r")

        t = load(o_fn,"tout")
        psi = load(o_fn,"Ψt")

        b = StimulatedBoseHubbard.LatticeBasis(8,8,8)

        func1 = a -> unzip(a)[1]
        func2 = a -> unzip(a)[2]
        func3 = a -> unzip(a)[3]
        func4 = a -> unzip(a)[4]

        ns1 = [func1(StimulatedBoseHubbard.averages(b,p)) for p in psi]
        ns2 = [func2(StimulatedBoseHubbard.averages(b,p)) for p in psi]
        ns3 = [func3(StimulatedBoseHubbard.averages(b,p)) for p in psi]
        ns4 = [func4(StimulatedBoseHubbard.averages(b,p)) for p in psi]

        push!(number,ns1)
        push!(phase,ns2)
        push!(number_var,ns3)
        push!(phase_var,ns4)

        tout = t
    end

    jldsave(string(data_dir,"/","variances",ann,".jld2");number,phase,number_var,phase_var,tout)
end

processNumber(fn,"")
processNumber(fn,"-stim")
processNumber(fn,"-phase")
