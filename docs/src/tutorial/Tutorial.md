# Modia Tutorial

This tutorial gives an overview of package [Modia](https://github.com/ModiaSim/Modia.jl)
to construct component-based and equation-based models with the **Modia language**
on a high level, symbolically transforming these models into ODEs
(Ordinary Differential Equations in state space form), simulating them and plotting result variables.

Note, all examples in this tutorial can be executed with\
`using Modia; include("$(Modia.path)/examples/Tutorial.jl")`

Modeling of 3D components is explained in the [Modia3D Tutorial](https://modiasim.github.io/Modia3D.jl/stable/)

!!! info
    Modia has an interface to various plot packages. A plot package can be
    either selected by setting `ENV["MODIA_PLOT"] = XXX`, for example in the `config/startup.jl`
    file of Julia or by command `Modia.usePlotPackage(XXX)`. Possible values for `XXX`: 
    - "[PyPlot](https://github.com/JuliaPy/PyPlot.jl)" (plots with Matplotlib from Python), 
    - "[GLMakie](https://github.com/JuliaPlots/GLMakie.jl)" (interactive plots in an OpenGL window),
    - "[WGLMakie](https://github.com/JuliaPlots/WGLMakie.jl)" (interactive plots in a browser window),
    - "[CairoMakie](https://github.com/JuliaPlots/CairoMakie.jl)" (static plots on file with publication quality),
    - "NoPlot" (= all `plot(...)` calls are ignored), or
    - "SilentNoPlot" (= NoPlot without messages).
    
