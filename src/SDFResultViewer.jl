module SDFResultViewer

export add_λ_units, unit_l, unit_t, unit_E, unit_B, unit_L, pₑ,
    average_linear_density_x, npart_plot,
    scalar_field_widget

using SDFResults
using PICDataStructures
using PICDataStructures: dir_to_idx
using PICAnalysisTools
using Unitful
using PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
using Statistics, StatsBase
# plots
import Plots
import PyPlot
using UnitfulRecipes
# widgets
using WGLMakie
using JSServe
using JSServe: DOM, Slider
using Markdown

include("units.jl")
include("plots.jl")
include("widgets.jl")

end
