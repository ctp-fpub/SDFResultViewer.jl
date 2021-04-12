@doc """
    powertostring(base, number)
Converts a number b^p into a string "bᵖ".
#Example
```julia-repl
julia> powertostring(2, 2^BigInt(1345))
"2¹³⁴⁵"
```
""" powertostring!
function powertostring(base, number)
    power = Int(round(log(base,number))) #scale = base^power
    s = "$base"*superscript(power)
    return s
end

@doc """
    photon_dens_integrated(sim; species = "photon", direction = direction, plot_title = title)

Plot of photon number density along a given direction (direction::Symbol argument; default is :x).
The distribution is normalized to the number density of the electrons at t=0.
""" photon_dens_integrated!
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

@doc """
    photon_en_dens_integrated(sim; species = "photon", direction = direction, endens_scale = scale , plot_title = title)

Plot of photon energy density along a given direction (direction::Symbol argument; default is :x) in MeV/μm².
One can choose to add a 10ˣ scale to the units with the endens_scale argument (pro tip: use BigInt for large powers).
""" photon_en_dens_integrated!
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
        a1 = ustrip(u"μm", cell_length(file, :y))
        a2 = ustrip(u"μm", cell_length(file, :z))
        nᵧ_int = dropdims(sum(n, dims=dim), dims=dim) / (a1*a2) / endens_scale
        plt = Plots.heatmap([1:1:size(nᵧ_int,1)]*a1, [1:1:size(nᵧ_int,2)]*a2, transpose(nᵧ_int); kwargs...)
        Plots.plot!(plt,
              xlabel = "y (μm)",
              ylabel = "z (μm)",
        )
        return plt
    elseif direction == :y
        a1 = ustrip(u"μm", cell_length(file, :x))
        a2 = ustrip(u"μm", cell_length(file, :z))
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
    photon_energy_spectrum(sim; species = photon_species, plot_title = title, angle = pi)

Photon energy spectra sᵧ² dN/dsᵧ = f(sᵧ = log₁₀(E/keV)). The photon species can be selected.
One can filter photons in cone. default angle is π (180⁰), i.e. all photons going to the right
""" photon_energy_spectrum!
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
    energy_spectrum(sim; mass = m_e, species = "electron", plot_title = title)

Energy spectra s² dN/ds = f(s = log₁₀(E/MeV)) for particles with mass.
The mass of the particle has to be added by hand with the mass argument (default is the electron mass).
""" energy_spectrum!
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
    photon_solid_angle_emission(sim;  angle = pi, species = photon_species, plot_title = title)

Angular distribution of gamma energy emission (MeV/srad). The angles are given by the direction of the momentum
(the reference axis is Ox). One can filter photons in cone. default angle is π (180⁰), i.e. all photons going
to the right.
""" photon_solid_angle_emission!
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
