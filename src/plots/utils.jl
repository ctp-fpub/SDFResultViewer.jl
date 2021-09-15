struct MinMaxTicks end

function MakieLayout.get_tickvalues(::MinMaxTicks, min, max)
    ticks = MakieLayout.get_tickvalues(WilkinsonTicks(5, k_min = 3), min, max)
    pushfirst!(ticks, min)
    push!(ticks, max)
    ticks
end
