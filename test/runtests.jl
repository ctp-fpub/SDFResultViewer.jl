using SDFResultViewer
using Unitful
using Test

@testset "SDFResultViewer.jl" begin
    λ = 800u"nm"
    dir = pwd()
    file = "0002.sdf"
    times = get_times(dir, λ)
    @test length(times) == 1

    m_Lx = mean_Lx(dir, λ)
    @test length(m_Lx) == 1

    Lx, y, z = Lx_slice(joinpath(dir, file), λ, 0)
    @test length(Lx) > 0
    @test length(Lx) == length(y) == length(z)
end
