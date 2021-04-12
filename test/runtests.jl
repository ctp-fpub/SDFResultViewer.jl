using SDFResultViewer
using Unitful
using SDFResults
using Test

@testset "SDFResultViewer.jl" begin
    if isempty(ENV["CI"])
        dir = "/mnt/storage/epoch/Wang-Arefiev2020mod/12mtt_0.8w0_l0_MDO"
    else
        dir = pwd()
    end
    sim = read_simulation(dir)

    test_plot1 = SDFResultViewer.photon_energy_spectrum(sim; mark_max = true)

    test_plot21 = SDFResultViewer.photon_dens_integrated(sim; direction = :y)
    test_plot22 = SDFResultViewer.photon_dens_integrated(sim; direction = :z)

    test_plot31 = SDFResultViewer.photon_en_dens_integrated(sim; direction = :y, endens_scale = 10^3)
    test_plot32 = SDFResultViewer.photon_en_dens_integrated(sim; direction = :z, endens_scale = 10^3)

    test_plot4 = SDFResultViewer.energy_spectrum(sim)

    test_plot5 = SDFResultViewer.photon_solid_angle_emission(sim; angle = pi/4)

    @test typeof(test_plot1) == Plots.Plot{Plots.GRBackend}
    @test typeof(test_plot21) == Plots.Plot{Plots.GRBackend}
    @test typeof(test_plot22) == Plots.Plot{Plots.GRBackend}
    @test typeof(test_plot32) == Plots.Plot{Plots.GRBackend}
    @test typeof(test_plot31) == Plots.Plot{Plots.GRBackend}
    @test typeof(test_plot4) == Plots.Plot{Plots.GRBackend}

    @test typeof(test_plot5) == PyPlot.Figure

    @test SDFResultViewer.powertostring(2, 2^BigInt(1345)) == "2¹³⁴⁵"
    @test SDFResultViewer.powertostring(10, 1000) == "10³"

    @test SDFResultViewer.cell_length(sim, :x) == 15u"μm"/450

end
