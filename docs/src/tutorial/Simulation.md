# 3 Simulation

A particular model is instantiated, simulated and results plotted with the commands:

```julia
using Modia
using DifferentialEquations
@usingModiaPlot

filter = @instantiateModel(Filter)
simulate!(filter, stopTime=10.0)
plot(filter, "y", figure=1)
```


## 3.1 Instantiating

The `@instantiateModel` macro takes additional arguments:

```julia
modelInstance = @instantiateModel(model;
                    FloatType = Float64, aliasReduction=true, unitless=false,
                    log=false, logModel=false, logDetails=false, logStateSelection=false,
                    logCode=false, logExecution=false, logTiming=false)
```

The macro performs structural and symbolic transformations, generates a function for
calculation of derivatives suitable for use with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
and returns [`SimulationModel`](@ref) that can be used in other functions,
for example to simulate or plot results:

* `model`: model (declarations and equations).
* `FloatType`: Variable type for floating point numbers (see below).
* `aliasReduction`: Perform alias elimination and remove singularities.
* `unitless`: Remove units (useful while debugging models and needed for MonteCarloMeasurements).
* `log`: Log the different phases of translation.
* `logModel`: Log the variables and equations of the model.
* `logDetails`: Log internal data during the different phases of translation.
* `logStateSelection`: Log details during state selection.
* `logCode`: Log the generated code.
* `logExecution`: Log the execution of the generated code (useful for timing compilation)
* `logCalculations`: Log the calculations of the generated code (useful for finding unit bugs)
* `logTiming`: Log timing of different phases.
* `return modelInstance prepared for simulation`

## 3.2 Simulating

The [`simulate!`](@ref) function performs one simulation with
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) using the default integrator
that this package automatically selects and stores the result in `modelInstance`.
 It is also possible to specify the integrator as second argument of `simulate!`:

```julia
using Modia
using DifferentialEquations
@usingModiaPlot

filter = @instantiateModel(Filter)
sol = simulate!(filter, Tsit5(), stopTime=10.0, merge=Map(T=0.5, x=0.8))
plot(filter, ["y", "x"], figure=1)
```

Integrator `DifferentialEquations.Tsit5` is an
[adaptive Runge-Kutta method of order 5/4 from Tsitouras](https://www.sciencedirect.com/science/article/pii/S0898122111004706).
There are > 100 ODE integrators provided. For details, see [here](https://docs.sciml.ai/stable/solvers/ode_solve/).

Parameters and init/start values can be changed with the `merge` keyword.
The effect is the same, as if the filter would have been instantiated with:

```julia
filter = @instantiateModel(Filter | Map(T=0.5, x=Var(init=0.8))
```

Note, with the `merge` keyword in simulate!, init/start values are directly
given as a value (`x = 0.8`) and are not defined with `Var(..)`.

Function `simulate!` returns `sol` which is the value that is returned by function
[DifferentialEquations.solve](https://diffeq.sciml.ai/stable/features/ensemble/#Solving-the-Problem).
Functions of `DifferentialEquations` that operate on this return argument can therefore also be
used on the return argument `sol` of `simulate!`.


## 3.4 Plotting

A short overview of the most important plot commands is given in section
section [Results and Plotting](@ref)


## 3.5 State selection (DAEs)

Modia has a sophisticated symbolic engine to transform high index DAEs
(Differential Algebraic Equations) automatically to ODEs (Ordinary Differential Equations in
state space form). During the transformation, equations might be analytically
differentiated and code might be generated to solve linear equation systems
numerically during simulation. The current engine **cannot** transform a DAE to ODE form, if the
**DAE contains nonlinear algebraic equations**. There is an (internal) prototype available to transform
nearly any DAE system to a special index 1 DAE system that can be solved with standard DAE integrators.
After a clean-up phase, this engine will be made publicly available at some time in the future.
Some of the algorithms used in Modia are
described in [Otter and Elmqvist (2017)](https://modelica.org/events/modelica2017/proceedings/html/submissions/ecp17132565_OtterElmqvist.pdf). Some algorithms are not yet published.

Usually, the symbolic engine is only visible to the modeler, when the model has errors,
or when the number of ODE states is less than the number of DAE states. The latter case
is discussed in this section.

The following object diagram shows two rotational inertias that are connected
by an ideal gear. One inertia is actuated with a sinusoidal torque:

![TwoInertiasAndIdealGear](../../resources/images/TwoInertiasAndIdealGear.png)

In order to most easily understand the issues, this model
is provided in a compact, "flattened" form:

```julia
TwoInertiasAndIdealGearTooManyInits = Model(
    J1    = 50,
    J2    = 100,
    ratio = 2,
    f     = 3, # Hz

    phi1 = Var(init = 0.0), # Absolute angle of inertia1
    w1   = Var(init = 0.0), # Absolute angular velocity of inertia1
    phi2 = Var(init = 0.0), # Absolute angle of inertia2
    w2   = Var(init = 0.0), # Absolute angular velocity of inertia2

    equations = :[
        tau = 2.0*sin(2*3.14*f*time/u"s")

        # inertia1
        w1 = der(phi1)
        J1*der(w1) = tau - tau1

        # ideal gear
        phi1 = ratio*phi2
        ratio*tau1 = tau2

        # inertia2
        w2 = der(phi2)
        J2*der(w2) = tau2
    ]
)

drive1 = @instantiateModel(TwoInertiasAndIdealGearTooManyInits)
simulate!(drive1, Tsit5(), stopTime = 1.0, logStates=true)
plot(drive1, [("phi1", "phi2"), ("w1", "w2")])
```

The option `logStates=true` results in the following output:

```
... Simulate model TwoInertiasAndIdealGearTooManyInits

│ # │ state  │ init │ unit │ nominal │
├───┼────────┼──────┼──────┼─────────┤
│ 1 │ phi2   │ 0.0  │      │ NaN     │
│ 2 │ w2     │ 0.0  │      │ NaN     │
```

This model translates and simulates without problems.

Changing the init-value of `w2` to `1.0` and re-simulating:

```julia
simulate!(drive1, Tsit5(), stopTime = 1.0, logStates=true, merge = Map(w2=1.0))
```

results in the following error:

```
... Simulate model TwoInertiasAndIdealGearTooManyInits

│ # │ state │ init │ unit │ nominal │
├───┼───────┼──────┼──────┼─────────┤
│ 1 │ phi2  │ 0.0  │      │ NaN     │
│ 2 │ w2    │ 1.0  │      │ NaN     │


Error from simulate!:
The following variables are explicitly solved for, have init-values defined
and after initialization the init-values are not respected
(remove the init-values in the model or change them to start-values):

│ # │ name │ beforeInit │ afterInit │
├───┼──────┼────────────┼───────────┤
│ 1 │ w1   │ 0.0        │ 2.0       │
```

The issue is the following:

Every variable that is used in the `der` operator is a **potential ODE state**.
When an `init` value is defined for such a variable, then Modia either utilizes
this initial condition (so the variable has this value after initialization), or an
error is triggered, as in the example above.

The model contains the equation:

```julia
phi1 = ratio*phi2
```

So the potential ODE states `phi1` and `phi2` are constrained, and only one of them
can be selected as ODE state, and the other variable is computed from this equation.
Since parameter `ratio` can be changed before simulation is started, it can be changed
also to a value of `ratio = 0`. Therefore, only when `phi2` is selected as ODE state,
`phi1` can be uniquely computed from this equation. If `phi1` would be selected as ODE state,
then a division by zero would occur, if `ratio = 0`, since `phi2 = phi1/ratio`. For this
reason, Modia selects `phi2` as ODE state. This means the **`init` value of `phi1`
has no effect**. This is uncritical, as long as initialization computes this init value
from the constraint equation above, as done in the example above.

When differentiating the equation above:

```julia
der(phi1) = ratio*der(phi2)  # differentiated constraint equation
       w1 = der(phi1)
       w2 = der(phi2)
```

it becomes obvious, that there is also a hidden constraint equation for 
the potential ODE states `w1, w2`:

```julia
w1 = ratio*w2  # hidden constraint equation
```

Again, Modia selects `w2` as ODE state, and ignores the `init` value of `w1`.
In the second simulation, the `init` value of `w1` (= 0.0) is no longer consistent to the
init value of `w2` (= 1.0). Therefore, an error occurs.

The remedy is to remove the `init` values of `phi1, w1` from the model:

```julia
drive2 = @instantiateModel(TwoInertiasAndIdealGearTooManyInits |
                                Map(phi1 = Var(init=nothing),
                                    w1   = Var(init=nothing)) )
simulate!(drive2, Tsit5(), stopTime = 1.0, logStates=true, merge = Map(w2=1.0))
```

and simulation is successful!

Modia tries to respect `init` values during symbolic transformation.
In cases as above, this is not possible and the reported issue occurs.
In some cases, it might not be obvious, why Modia selects a particular
variable as an ODE state. You can get more information by setting
`logStateSelection=true`:

```julia
drive1 = @instantiateModel(TwoInertiasAndIdealGearTooManyInits, logStateSelection=true)
```

