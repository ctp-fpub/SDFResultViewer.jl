module SDFResultViewer

export compute_L, compute_Lx,
    mean_quantity, approximate_field,
    add_λ_units, unit_l, unit_t, unit_E, unit_B, unit_L, pₑ,
    average_linear_density_x, npart_plot,
    scalar_field_widget

using SDFResults
using PICDataStructures
using PICDataStructures: dir_to_idx
using RecursiveArrayTools
using LinearAlgebra
using Unitful
using PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
using StaticArrays
using FileTrees
using NumericalIntegration
using Statistics
using Transducers, ThreadsX

include("units.jl")
include("analysis.jl")
include("number_density.jl")
include("plots.jl")
include("widgets.jl")

end
