function add_λ_units(λ)
    expr = quote

        using Unitful
        using Unitful: MW
        using PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e

        local λ = $λ
        local ω = 2π*c_0/λ

        @unit unit_l "c/ω" UnitLength c_0/ω false
        @unit unit_t "1/ω" UnitTime 1/ω false
        @unit unit_B "mₑ ω/e" UnitMagneticField m_e * ω / e false
        @unit unit_E "mₑ c ω/e" UnitElectricField m_e * ω / e * c_0 false
        @unit unit_L "mₑ c²/ω" UnitAngularMomenta m_e * c_0^2 / ω false
        @unit unit_S "MW ω²/c²" UnitPoyntingFlux MW * ω^2 / c_0^2 false
        @unit pₑ "mₑ c" ElectronMomenta m_e * c_0 false

        Unitful.register(@__MODULE__)
    end
    eval(expr)
end

function add_λ_units(sim::EPOCHSimulation)
    λ = get_parameter(sim, :laser, :lambda) |> u"nm"
    add_λ_units(λ)
end
