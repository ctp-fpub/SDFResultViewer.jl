function npart_plot(sim; species="electron")
    ts = get_time.(sim) .|> unit_t
    npart = get_npart.(sim, (species,))

    Plots.plot(ts, npart, xlabel="t", ylabel="N")
end
