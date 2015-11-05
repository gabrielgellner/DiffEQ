## Test Suite
As this package is doing reimplementation of ODE solvers in pure Julia it is
essential that it has a good test suite. Thing to look at:
* The DETEST paper and its examples
* The Julia package `IVPTestSuite.jl`
I will want to have this level of functionality at a minimum. Some other ideas
is to have a general sweet of linear problems where we can compare the analytic
to the approximate solutions. As well as a set of problems that we can compare
to solvers from things like `Matlab` and `Mathematica` which are known good
solvers.
