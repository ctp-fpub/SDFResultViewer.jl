module SDFResultViewer

export add_λ_units, unit_l, unit_t, unit_E, unit_B, unit_L, pₑ,
    batch_plot,
    mean_Lx_plot, Lx_section_plot,
    phase_space_summary_plot,
    average_linear_density_x, npart_plot,
    section_widget

using SDFResults
using PICDataStructures
using PICDataStructures: dir_to_idx
using PICAnalysisTools
using Unitful
using PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
using Statistics, StatsBase
using ProgressLogging
# plots
import Plots
import PyPlot
# widgets
using WGLMakie
using JSServe
using JSServe: DOM, Slider
using Markdown

include("units.jl")
include("batch.jl")
include("plots.jl")
include("widgets.jl")

end
