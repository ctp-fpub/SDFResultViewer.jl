function phase_space_summary_plot(sim, dir; species="electron", save=false)
    r_mean, p_mean,
    r_mean_plus, p_mean_plus,
    r_mean_minus, p_mean_minus = @withprogress name="Computing phase space mean" begin
        r_mean, p_mean = phase_space_mean(sim, dir; species)
        @logprogress 1/3
        r_mean_plus, p_mean_plus = phase_space_mean(sim, dir; species, rcond=(i,r)->r>r_mean[i])
        @logprogress 2/3
        r_mean_minus, p_mean_minus = phase_space_mean(sim, dir; species, rcond=(i,r)->r<r_mean[i])
        @logprogress 1

        r_mean, p_mean, r_mean_plus, p_mean_plus, r_mean_minus, p_mean_minus
    end

    ts = ustrip.(unit_t, get_time.(sim))
    r = uconvert.(unit_l, r_mean)
    r_plus = uconvert.(unit_l, r_mean_plus)
    r_minus = uconvert.(unit_l, r_mean_minus)
    rm, rmp, rmm = mean(r), mean(r_plus), mean(r_minus)
    p = uconvert.(pₑ, p_mean)
    p_plus = uconvert.(pₑ, p_mean_plus)
    p_minus = uconvert.(pₑ, p_mean_minus)
    pm, pmp, pmm = mean(p), mean(p_plus), mean(p_minus)

    r_label = "$dir - $(dir)ₘ"
    plt1 = Plots.plot(ts, r .- rm, xlabel="t", ylabel=r_label, label="all")
    Plots.plot!(plt1, ts, r_plus .- rmp, xlabel="t", ylabel=r_label, label="r > r̅")
    Plots.plot!(plt1, ts, r_minus .- rmm, xlabel="t", ylabel=r_label, label="r < r̅")

    plt2 = Plots.plot(ts, p, xlabel="t", ylabel="p$dir", label="all")
    Plots.plot!(plt2, ts, p_plus, xlabel="t", ylabel="p$dir", label="r > r̅")
    Plots.plot!(plt2, ts, p_minus, xlabel="t", ylabel="p$dir", label="r < r̅")

    plt3 = Plots.scatter(r, p, zcolor=ts, xlabel=string(dir), ylabel="p$dir", legend=false)

    plt4 = Plots.scatter(r .- rm, p, xlabel=r_label, ylabel="p$dir", label="all")
    Plots.scatter!(r_plus .- rmp, p_plus, xlabel=r_label, ylabel="p$dir", label="r > r̅")
    Plots.scatter!(r_minus .- rmm, p_minus, xlabel=r_label, ylabel="p$dir", label="r < r̅")

    plt = Plots.plot(plt1, plt2, plt3, plt4, layout=4, size=(900,700))

    save && Plots.savefig(plt, joinpath(sim.dir, "phase_space_summary.png"))

    return plt
end
