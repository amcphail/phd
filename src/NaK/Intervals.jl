module Intervals

using Statistics, StatsBase

export measure_intervals

function create_hist(series,bins)
    intervals = measure_intervals(series)

    @show typeof(intervals)

    m = maximum(intervals)

    edges = 1:floor(m/bins):m

    h = fit(Histogram,Float64.(intervals),edges)

    return (h.edges,h.weights)
end

function measure_intervals(series)
    len = size(series)[1]

    intervals = []

    start = 0

    for i = 1:len
        if series[i] == 1
            start = i
            break
        end
    end

    if start == 0
        return intervals
    end

    for i = (start+1):len
        if series[i] == 1
            interval = i - start
            push!(intervals,interval)
            start = i
        end
    end

    return intervals
end

end
