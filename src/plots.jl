@doc """
    powertostring(base, number)

Converts a number b^p into a string "bᵖ".

#Example
```julia-repl
julia> powertostring(2, 2^BigInt(1345))
"2¹³⁴⁵"
```
"""
function powertostring(base, number)
    power = Int(round(log(base,number))) #scale = 10^power
    d = Dict(0 => "⁰",
             1 => "¹",
             2 => "²",
             3 => "³",
             4 => "⁴",
             5 => "⁵",
             6 => "⁶",
             7 => "⁷",
             8 => "⁸",
             9 => "⁹"
    ) # tried Dict([(i, "\^$i") for i in 0:1:9]), but didn't work
    s = "$base"
    for i in reverse(digits(power))
        s = s*d[i]
    end
    return s
end

@doc """
    photon_dens_integrated(sim)

Plot of photon number density along a given direction (direction::Symbol argument; default is :x).
The distribution is normalized to the number density of the electrons at t=0.
"""
function photon_dens_integrated(sim; species = "photon", direction = :x, plot_title = "Photon number density integrated along x-axis")
    file = sim[1]
    n₀ = file["number_density/electron"]

    file = sim[end]
    nᵧ = file["number_density/$species"]

    kwargs = (color=:rainbow,
              colorbar_title="nᵧ/n₀",
              tick_direction = :out,
              minorticks = true,
              title = plot_title
    )
    dim = dir_to_idx(direction)
    nᵧ_int = dropdims(sum(nᵧ*cell_volume(file), dims=dim), dims=dim)
    n₀_int = dropdims(sum(n₀*cell_volume(file), dims=dim), dims=dim)

    if direction == :x
        a1 = cell_length(file, :y)
        a2 = cell_length(file, :z)
        plt = Plots.heatmap([1:1:size(nᵧ_int,1)]*a1, [1:1:size(nᵧ_int,2)]*a2, transpose(nᵧ_int ./ n₀_int); kwargs...) #damn heatmap is dumbly implemented
        Plots.plot!(plt,
              xlabel = "y (μm)",
              ylabel = "z (μm)",
              size = (600,500)
        )
        return plt
    elseif direction == :y
        a1 = cell_length(file, :x)
        a2 = cell_length(file, :z)
        plt = Plots.heatmap([1:1:size(nᵧ_int,1)]*a1, [1:1:size(nᵧ_int,2)]*a2, transpose(nᵧ_int ./ n₀_int); kwargs...)
        Plots.plot!(plt,
              xlabel = "x (μm)",
              ylabel = "z (μm)",
              size = (665,500)
        )
        return plt
    elseif direction == :z
        a1 = cell_length(file, :x)
        a2 = cell_length(file, :y)
        plt = Plots.heatmap([1:1:size(nᵧ_int,1)]*a1, [1:1:size(nᵧ_int,2)]*a2, transpose(nᵧ_int ./ n₀_int); kwargs...)
        Plots.plot!(plt,
              xlabel = "x (μm)",
              ylabel = "y (μm)",
              size = (665,500)
        )
        return plt
    end
end

function photon_en_dens_integrated(sim; species = "photon", direction = :x, endens_scale = 10^12 , plot_title = "Photon energy density integrated along x-axis")
    file = sim[end]
    nᵧ, enᵧ= file["number_density/$species", "ekbar/$species"]
    enᵧ = ustrip.(u"MeV",enᵧ)

    scale = powertostring(10, endens_scale)
    kwargs = (color=:rainbow,
             colorbar_title="ϵᵧ ($scale MeV/μm²)",
             size = (665,600),
             tick_direction = :out,
             minorticks = true,
             title = plot_title
    )
    dim = dir_to_idx(direction)
    n = uconvert(NoUnits, nᵧ*cell_volume(file).*enᵧ)

    if direction == :x
        a1 = ustrip(u"μm",cell_length(file, :y))
        a2 = ustrip(u"μm",cell_length(file, :z))
        nᵧ_int = dropdims(sum(n, dims=dim), dims=dim) / (a1*a2) / endens_scale
        plt = Plots.heatmap([1:1:size(nᵧ_int,1)]*a1, [1:1:size(nᵧ_int,2)]*a2, transpose(nᵧ_int); kwargs...)
        Plots.plot!(plt,
              xlabel = "y (μm)",
              ylabel = "z (μm)",
        )
        return plt
    elseif direction == :y
        a1 = ustrip(u"μm",cell_length(file, :x))
        a2 = ustrip(u"μm",cell_length(file, :z))
        nᵧ_int = dropdims(sum(n, dims=dim), dims=dim) / (a1*a2) / endens_scale
        plt = Plots.heatmap([1:1:size(nᵧ_int,1)]*a1, [1:1:size(nᵧ_int,2)]*a2, transpose(nᵧ_int); kwargs...)
        Plots.plot!(plt,
              xlabel = "x (μm)",
              ylabel = "z (μm)",
        )
        return plt
    elseif direction == :z
        a1 = ustrip(u"μm",cell_length(file, :x))
        a2 = ustrip(u"μm",cell_length(file, :y))
        nᵧ_int = dropdims(sum(n, dims=dim), dims=dim) / (a1*a2) / endens_scale
        plt = Plots.heatmap([1:1:size(nᵧ_int,1)]*a1, [1:1:size(nᵧ_int,2)]*a2, transpose(nᵧ_int); kwargs...)
        Plots.plot!(plt,
              xlabel = "x (μm)",
              ylabel = "y (μm)",
        )
        return plt
    end
end

@doc """
Photon energy spectra sᵧ² dN/dsᵧ = f(sᵧ = log₁₀(E/keV)). The photon species can be selected.
One can filter photons in cone. default angle is π (180⁰), i.e. all photons going to the right
"""
function photon_energy_spectrum(sim; label = "label", mark_max = false, angle = pi, species = "photon", plot_title = "Compton scattering spectrum")
    file = sim[end]
    px, py, pz = file["px/$species", "py/$species", "pz/$species"]
    # filtering data: keep only photons moving in the direction of the laser
    idx = findall(p -> p > zero(eltype(px)), px)
    px = px[idx]
    py = py[idx]
    pz = pz[idx]

    # θ is measured from the x axis
    # filtering data: keep photons that propagate within 2θ
    pₚₑᵣₚ = sqrt.(py.^2 + pz.^2)   #transversal momentum of each photon
    θ = atan.(pₚₑᵣₚ./px)           #emission angles for each photon
    idx = findall(θ₀ -> θ₀ < angle/2, θ)
    px = px[idx]
    py = py[idx]
    pz = pz[idx]

    # make plot
    Eᵧ = sqrt.(px.^2 + py.^2 + pz.^2) * c_0
    sᵧ = log.(ustrip(u"keV", Eᵧ))

    h = fit(Histogram, sᵧ; nbins = 200)
    s = midpoints(collect(h.edges[1]))
    N = h.weights
    dNds = N./abs(s[1]-s[2])

    idx = findall(x -> x > 0, s)
    s = s[idx]
    dNds = dNds[idx]

    kwargs = (xlabel = "sᵧ = log₁₀(Eᵧ/keV)",
              ylabel = "sᵧ² dN/dsᵧ",
              label = label,
              tickfontsize = 10,
              guidefontsize = 15,
              title = plot_title,
              legendfontsize = 15,
              legendtitlefont = 15,
              titlefont = 16
    )

    plt = Plots.plot(s, dNds.*s.*s, linewidth = 4; kwargs...)

    if mark_max
        dNds = dNds.*s.*s
        p0 = findfirst(smax -> smax == maximum(dNds), dNds)
        Plots.plot!(plt, [s[p0]], [dNds[p0]], seriestype = :scatter,
                                    markersize = 7,
                                    markercolor = :black,
                                    label = false
        )
    end

    return plt
end

# an alternative version of photon_energy_spectrum() that plots multiple spectra on the same plot
function photon_energy_spectrum(sims, labels; legend_title = "labels", mark_max = false, angle = pi, species = "photon", plot_title = "Compton scattering spectrum")
    kwargs = (xlabel = "sᵧ = log₁₀(Eᵧ/keV)",
            ylabel = "sᵧ² dN/dsᵧ",
            tickfontsize = 10,
            guidefontsize = 15,
            title = plot_title,
            legendfontsize = 15,
            legendtitlefont = 15,
            titlefont = 16,
            legendtitle = legend_title
    )
    plt = Plots.plot(;kwargs...)
    for (i, sim) in enumerate(sims)
        file = sim[end]
        px, py, pz = file["px/$species", "py/$species", "pz/$species"]
        # filtering data: keep only photons moving in the direction of the laser
        idx = findall(p -> p > zero(eltype(px)), px)
        px = px[idx]
        py = py[idx]
        pz = pz[idx]

        # θ is measured from the x axis
        # filtering data: keep photons that propagate within 2θ
        pₚₑᵣₚ = sqrt.(py.^2 + pz.^2)   #transversal momentum of each photon
        θ = atan.(pₚₑᵣₚ./px)           #emission angles for each photon
        idx = findall(θ₀ -> θ₀ < angle/2, θ)
        px = px[idx]
        py = py[idx]
        pz = pz[idx]

        # make plot
        Eᵧ = sqrt.(px.^2 + py.^2 + pz.^2) * c_0
        sᵧ = log.(ustrip(u"keV", Eᵧ))

        h = fit(Histogram, sᵧ; nbins = 200)
        s = midpoints(collect(h.edges[1]))
        N = h.weights
        dNds = N./abs(s[1]-s[2])

        idx = findall(x -> x > 0, s)
        s = s[idx]
        dNds = dNds[idx]
        plt = Plots.plot!(plt, s, dNds.*s.*s, linewidth = 4; label = labels[i])

        if mark_max
            dNds = dNds.*s.*s
            p0 = findfirst(smax -> smax == maximum(dNds), dNds)
            Plots.plot!(plt, [s[p0]], [dNds[p0]], seriestype = :scatter,
                                    markersize = 7,
                                    markercolor = :black,
                                    label = false
            )
        end
    end
    return plt
end

@doc """
Energy spectra s² dN/ds = f(s = log₁₀(E/MeV)) for particles with mass.
The mass of the particle has to be added by hand with the mass argument (default is the electron mass).
"""
function energy_spectrum(sim; mass = m_e, species = "electron", plot_title = "Electron energy spectrum")
    file = sim[end]
    px, py, pz = file["px/$species", "py/$species", "pz/$species"]

    #should be implemented: read species mass from input.deck
    E₀ = sqrt.(mass^2*c_0^4 .+ (px.^2 + py.^2 + pz.^2) * c_0^2)
    s₀ = log.(ustrip(u"MeV", E₀))

    h = fit(Histogram, s₀; nbins = 200)
    s = midpoints(collect(h.edges[1]))
    N = h.weights
    dNds = N./abs(s[1]-s[2])

    idx = findall(x -> x > 0, s)
    s = s[idx]
    dNds = dNds[idx]

    kwargs = (xlabel = "s = log₁₀(E/MeV)",
              ylabel = "s² dN/ds",
              tickfontsize = 10,
              guidefontsize = 15,
              title = plot_title,
              legendfontsize = 15,
              legendtitlefont = 15,
              titlefont = 16
    )

    plt = Plots.plot(s, dNds.*s.*s, linewidth = 4; kwargs...)

    return plt
end

# an alternative version of energy_spectrum() that plots multiple spectra on the same plot
function energy_spectrum(sims, labels; legend_title = "labels", mass = m_e, species = "electron", plot_title = "Electron energy spectrum")
    kwargs = (xlabel = "s = log₁₀(E/MeV)",
              ylabel = "s² dN/ds",
              tickfontsize = 10,
              guidefontsize = 15,
              title = plot_title,
              legendfontsize = 15,
              legendtitlefont = 15,
              titlefont = 16,
              legendtitle = legend_title
    )
    plt = Plots.plot(;kwargs...)
    for (i, sim) in enumerate(sims)
        file = sim[end]
        px, py, pz = file["px/$species", "py/$species", "pz/$species"]

        #should be implemented: read species mass from input.deck
        E₀ = sqrt.(mass^2*c_0^4 .+ (px.^2 + py.^2 + pz.^2) * c_0^2)
        s₀ = log.(ustrip(u"MeV", E₀))

        h = fit(Histogram, s₀; nbins = 200)
        s = midpoints(collect(h.edges[1]))
        N = h.weights
        dNds = N./abs(s[1]-s[2])

        idx = findall(x -> x > 0, s)
        s = s[idx]
        dNds = dNds[idx]

        plt = Plots.plot!(plt, s, dNds.*s.*s, linewidth = 4; label = labels[i])
    end
    return plt
end

@doc """
Angular distribution of gamma energy emission (MeV/srad). The angles are given by the direction of the momentum
(the reference axis is Ox). One can filter photons in cone. default angle is π (180⁰), i.e. all photons going
to the right.
"""
function photon_solid_angle_emission(sim;  angle = pi, species = "photon", plot_title = "Solid angle emission")
    file = sim[end]

    px, py, pz= file["px/$species",
                     "py/$species",
                     "pz/$species"
    ]

    idx = findall(p -> p > zero(eltype(px)), px)
    px = px[idx]
    py = py[idx]
    pz = pz[idx]

    pₚₑᵣₚ = sqrt.(py.^2 + pz.^2)   #transversal momentum of each photon
    θ = atan.(pₚₑᵣₚ./px)           #emission angles for each photon
    ϕ = atan.(py, pz)

    # filtering data: keep photons that propagate within 2θ
    idx = findall(θ₀ -> θ₀ < angle/2, θ)
    px = px[idx]
    py = py[idx]
    pz = pz[idx]
    θ = θ[idx]
    ϕ = ϕ[idx]

    Eᵧ = sqrt.(px.^2 + py.^2 + pz.^2) * c_0

    data = (ϕ .+ π, θ)
    w = weights(ustrip.(u"MeV", Eᵧ))

    h = fit(Histogram, data, w, nbins = (360,90))
    xedges = collect(h.edges[1])
    yedges = rad2deg.(h.edges[2])
    H = h.weights

    ax = PyPlot.Axes3D(PyPlot.figure())
    PyPlot.subplot(projection="polar")
    PyPlot.pcolormesh(xedges[2:end], yedges[2:end], transpose(H))
    PyPlot.colorbar(label = "dEᵧ/dΩ (MeV/srad)")
    PyPlot.title(plot_title)
    ax[:grid](true)
    # PyPlot.display(PyPlot.gcf())
    return PyPlot.gcf()
end

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
    ts = get_time.(sim) .|> unit_t
    npart = get_npart.(sim, (species,))

    Plots.plot(ts, npart, xlabel="t", ylabel="N")
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
