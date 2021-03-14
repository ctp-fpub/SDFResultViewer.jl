# SDFResultViewer

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SebastianM-C.github.io/SDFResultViewer.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SebastianM-C.github.io/SDFResultViewer.jl/dev)
[![Build Status](https://github.com/SebastianM-C/SDFResultViewer.jl/workflows/CI/badge.svg)](https://github.com/SebastianM-C/SDFResultViewer.jl/actions)
[![Coverage](https://codecov.io/gh/SebastianM-C/SDFResultViewer.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SebastianM-C/SDFResultViewer.jl)

This package aims to provide visualization tools for [EPOCH](https://cfsa-pmw.warwick.ac.uk/mediawiki/index.php/EPOCH:FAQ) simulations.

## Installation

In order to use this package you need to add the [ctp-fpub Julia registry](https://github.com/ctp-fpub/JuliaRegistry)
```
]registry add https://github.com:ctp-fpub/JuliaRegistry.git
```

After that, you can install the package with
```
]add SDFResultViewer
```

## Quck start

```julia
using SDFResults
using SDFResultViewer
using WGLMakie
using Unitful
using PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e

dir = "epoch_simulation"
sim = read_simulation(dir)
file = sim[4]

λ = get_parameter(file, :lambda)
ω = 2π * c_0 / λ
nc = (ε_0 * m_e) / e^2 * ω^2

add_λ_units(λ)
Ex = file[:ex]

scalar_field_widget(uconvert(unit_E, unit_l, Ex))
```

![widget](assets/widget.png)
