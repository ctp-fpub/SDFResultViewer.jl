import Plots
using UnitfulRecipes

function mean_Lx_plot(sim; species="electron")
    Lx(file) = compute_Lx(file, species) |> unit_L
    @time m_Lx = mean(Lx, sim)
    @time m_Lx_plus = mean(Lx, sim, cond=lx->lx>0)
    @time m_Lx_minus = mean(Lx, sim, cond=lx->lx<0)

    ts = get_time.(sim) .|> unit_t

    plt = Plots.plot(ts, m_Lx,
        title="Average Lx for l = $l",
        xlabel="t",
        ylabel="Lx",
        label="Total Lx")
    Plots.plot!(plt, ts, m_Lx_plus, label="Lx > 0")
    Plots.plot!(plt, ts, m_Lx_minus, label="Lx < 0")

    Plots.savefig(plt, joinpath(dir, "mean_Lx_l_$l.png"))
    plt
end

function Lx_section_plots(sim, slice_location = 54.97unit_l; ϵ = 2e-3unit_l, species="electron")
    if !ispath(joinpath(dir, "Lx_section"))
        mkdir(joinpath(dir, "Lx_section"))
    end

    for file in sim
        Lx = compute_Lx(file, species)
        Lx_section = sclice(Lx, :x, slice_location, ϵ)
        if all(iszero.(Lx))
            continue
        end
        cl = max(abs.(extrema(Lx))...)

        bg = Plots.cgrad(:RdYlBu_11, categorical=true)[6]
        t = sprint(show, get_time(file), context=:compact => true)

        plt = Plots.scatter(y, z, zcolor=Lx,
            xlabel = "y",
            ylabel = "z",
            clims = (-cl,cl),
            cb_title = "Lx",
            title = "Lx at x = $slice_location and t = $t",
            markersize = 2.5, mswidth = 0,
            color = :RdYlBu_11,
            background_inside = bg,
            framestyle = :box,
            label = "l = $l",
            aspect_ratio=1)

        Plots.savefig(plt, joinpath(dir, "Lx_section", "Lx_x$(slice_location)_eps$(ϵ)_t$t.png"))
    end
end

function npart_plot(sim; species="electron")
    ts = get_time.(sim) |> unit_t
    npart = get_npart.(sim, (species,))

    plt = Plots.plot(ts, npart, xlabel="t", ylabel="N")
    Plots.savefig(plt, joinpath(dir, "npart.png"))
    plt
end

# function center_of_mass_plot(dir, λ, d, species)
#     i = i_from_dir(d)
#     cm = ustrip.(u"μm", center_of_mass(dir, i, species))
#     cm_plus = ustrip.(u"μm", center_of_mass(dir, i, species, cond=ri->ri>0u"m"))
#     cm_minus = ustrip.(u"μm", center_of_mass(dir, i, species, cond=ri->ri<0u"m"))
#     ts = get_times(dir, λ)

#     plt = plot(ts, cm, xlabel = "t (1/ω₀)", ylabel = "$d-cm (micron)", title = species)
#     savefig(plt, joinpath(dir, "$d-cm.png"))

#     plt = plot(ts, cm_plus, xlabel = "t (1/ω₀)", ylabel = "$d-cm (micron)", title = species)
#     savefig(plt, joinpath(dir, "$d-cm_plus.png"))

#     plt = plot(ts, cm_minus, xlabel = "t (1/ω₀)", ylabel = "$d-cm (micron)", title = species)
#     savefig(plt, joinpath(dir, "$d-cm_minus.png"))

# end

# function center_of_momenta_plot(dir, λ, d, species)
#     unit_p = species == "electron" ? m_e * c_0 : 1836 * m_e * c_0

#     cm = uconvert.(NoUnits, center_of_momenta(dir, d, species)/unit_p)
#     cm_plus = uconvert.(NoUnits, center_of_momenta(dir, d, species, cond=p_i->p_i>zero(p_i))/unit_p)
#     cm_minus = uconvert.(NoUnits, center_of_momenta(dir, d, species, cond=p_i->p_i<zero(p_i))/unit_p)
#     ts = get_times(dir, λ)
#     ylab = species == "electron" ? "P$(d)cm (mₑ c)" : "P$(d)cm (mₚ c)"

#     plt = plot(ts, cm, xlabel = "t (1/ω₀)", ylabel = ylab, title = species)
#     savefig(plt, joinpath(dir, "cp$(d).png"))

#     plt = plot(ts, cm_plus, xlabel = "t (1/ω₀)", ylabel = ylab, title = species)
#     savefig(plt, joinpath(dir, "cp$(d)_plus.png"))

#     plt = plot(ts, cm_minus, xlabel = "t (1/ω₀)", ylabel = ylab, title = species)
#     savefig(plt, joinpath(dir, "cp$(d)_minus.png"))

# end

# function average_linear_density_x_plot(dir, λ)
#     files = filter(f->endswith(f, ".sdf"), readdir(dir))
#     ts = get_times(dir, λ)
#     ω = 2π*c_0/λ
#     unit_length = c_0 / ω

#     plt = plot(xlabel="x (c/ω)", ylabel="∫∫n (c/ω)⁻³")
#     for (i, file) in enumerate(files)
#         file = joinpath(dir, file)

#         x, nlx = average_linear_density_x(file, λ)
#         t = sprint(show, ts[i], context=:compact => true)
#         p = plot!(plt, x, nlx, label="t=$t")
#         savefig(p, joinpath(dir, "nlx_$i.png"))
#     end
#     savefig(plt, joinpath(dir, "nlx.png"))
#     plt
# end

# function average_linear_density_x_start_end_plot(dir, λ)
#     files = filter(f->endswith(f, ".sdf"), readdir(dir))
#     ts = get_times(dir, λ)
#     ti, i_i = findmin(ts)
#     te, i_e = findmax(ts)

#     ω = 2π*c_0/λ
#     unit_length = c_0 / ω

#     plt = plot(xlabel="x (c/ω)", ylabel="∫∫n (c/ω)⁻³")
#     ts = sprint.((show,), (ti,te), context=:compact => true)
#     for (file, t) in zip(files[[i_i, i_e]], ts)
#         file = joinpath(dir, file)

#         x, nlx = average_linear_density_x(file, λ)
#         plt = plot!(plt, x, nlx, label="t=$t")
#     end
#     savefig(plt, joinpath(dir, "nlx_start_end.png"))
#     plt
# end
