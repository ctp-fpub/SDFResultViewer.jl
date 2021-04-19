module SDFResultViewer

export add_λ_units, unit_l, unit_t, unit_E, unit_B, unit_L, unit_S, pₑ,
    batch_plot,
    mean_Lx_plot, Lx_section_plot,
    phase_space_summary_plot,
    average_linear_density_x, npart_plot,
    section_widget, approx_plot,
    compare_E_slice_with_analytic, compare_S_slice_with_analytic,
    poynting_summary,
    apply_analytic, analytic_laser

using SDFResults
using PICDataStructures
using PICDataStructures: dir_to_idx
using PICAnalysisTools
using LaserTypes
using StaticArrays
using Rotations
using RecursiveArrayTools: recursive_bottom_eltype
using Unitful
using Unitful: superscript
using PhysicalConstants.CODATA2018: c_0, ε_0, μ_0, m_e, e
using Statistics, StatsBase
using LinearAlgebra
using ProgressLogging
# plots
import Plots
import PyPlot
# widgets
using WGLMakie
using JSServe
using JSServe: DOM, Slider
using Markdown

Unitful.register(@__MODULE__)

include("units.jl")
include("batch.jl")
include("plots/angular_momentum.jl")
include("plots/fields.jl")
include("plots/phase_space.jl")
include("plots/qed.jl")
include("plots/statistics.jl")
include("widgets.jl")
include("utils.jl")

end
