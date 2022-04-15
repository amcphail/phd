module MutualInformation

using LinearAlgebra

using StatsBase

export mutual_information_dataset

function mutual_information(p1,p2,p12,bins)
    mi = 0

    for i in 1:bins
        for j in 1:bins
            if p12.weights[i,j] == 0 || p1.weights[i] == 0 || p2.weights[j] == 0
                Nothing
            else
                mi += p12.weights[i,j]*log2(p12.weights[i,j]/(p1.weights[i]*p2.weights[j]))
            end
        end
    end

    return mi
end

function mutual_information_dataset(data,lowest,highest,bins)
    (t,num_nodes) = size(data)

    edges = lowest:((highest-lowest)/bins):highest

    hists_1D = []
    hists_2D = []

    for i in 1:num_nodes
        push!(hists_1D,normalize(fit(Histogram,data[:,i],edges),mode=:probability))

        for j in 1:num_nodes
            push!(hists_2D,normalize(fit(Histogram,(data[:,i],data[:,j]),(edges,edges)),mode=:probability))
        end
    end

    hists_2D = reshape(hists_2D,(num_nodes,num_nodes))

    mi = []

    for i in 1:num_nodes
        for j in 1:num_nodes
            if j > i
                push!(mi,0)
            else
                push!(mi,mutual_information(hists_1D[i],hists_1D[j],hists_2D[i,j],bins))
            end
        end
    end

    mi = reshape(mi,(num_nodes,num_nodes))

    return mi
end

end
