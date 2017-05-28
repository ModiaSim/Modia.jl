# Modia
Modia is a domain specific extension of [Julia](http://julialang.org/ "Julia") for **modeling and simulation of physical systems**. The first version of Modia is planned to be uploaded in Summer 2017. Papers about Modia:

- [Overview about Modia](http://link.springer.com/chapter/10.1007%2F978-3-319-47169-3_15) ([ISoLA conference Oct. 2016](http://www.isola-conference.org/isola2016/)).
- [Overview and new features in Modia](https://www.modelica.org/events/modelica2017/proceedings/html/submissions/ecp17132693_ElmqvistHenningssonOtter.pdf) ([12th International Modelica Conference, Mai 2017](https://www.modelica.org/events/modelica2017/proceedings/html/index.html)).
- [New algorithms in Modia](https://www.modelica.org/events/modelica2017/proceedings/html/submissions/ecp17132565_OtterElmqvist.pdf) ([12th International Modelica Conference, Mai 2017](https://www.modelica.org/events/modelica2017/proceedings/html/index.html); slides for this paper in [pptx](https://modiasim.github.io/slides/Modelica2017-DAAE-Transformation.pptx) and [pdf](https://modiasim.github.io/slides/Modelica2017-DAAE-Transformation.pdf) format).

Modia is designed to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by "expression = expression" equations. The defined model is symbolically processed (for example, equations might be analytically differentiated), JIT compiled and simulated with [Sundials IDA solver](http://computation.llnl.gov/projects/sundials/ida) with the KLU sparse matrix package. By this approach it's possible and convenient to build models with hundred thousands of equations describing the dynamics of a car, an airplane, a power plant, etc. and simulate them. The authors used previous experience from the design of the modeling language [Modelica](https://www.modelica.org/). Modia will also be used to design and evaluate features for future Modelica versions. 

Component models are defined by @model macros. Such models contain definition of variables with various attributes such as start values, min, max, SI unit, etc. An @equations macro is used to define the equations of such a component. Coupling between components is expressed using a connect statement involving groups of variables. The semantics is either to constrain connected variables to be equal or to constrain a set of variables to sum to zero, for example to model Kirchhoff's current law. 

"Hello world" Modia model:

    # T*dx/dt + x = u(t)
    #
    using Modia
    @model FirstOrder begin
       x = Variable(start=1) 
       T = Parameter(0.5, "Time constant")
       u = 2.0 # Same as Parameter(2.0)
    @equations begin
       T*der(x) + x = u 
       end
    end
    result = simulateModel(FirstOrder, linspace(0,2,500))

