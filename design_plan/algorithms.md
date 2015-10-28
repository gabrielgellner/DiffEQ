## Planned Algorithms to support

### Initial Value Problem (IVP) Solvers
As my work is so heavily based on the code from `ODE.jl` I have an excellent
framework for including most of the commonly Runge-Kutta algorithms as they have
them nicely setup. This work is heavily based on the treatise by Hairer and
Wanner (1993).

Currently I am just using the:
* Dormand Prince (4, 5) method

I plan to add support for
* Dormand Prince (8, 5, 1) approach

I need to think if I want to add all the other RK methods that `ODE.jl` supports
like the Fehlberg set and the non adaptive sets. To start this is a very low
priority. The only thing that might be interesting about this is that I could
have it so that the user might have a way of specifying this combinations
instead of just calling a specific method symbol (something like that would
implement the idea of "Use the RK with order (2, 3) and adaptive stepper X").
This feels similar in some ways to how Mathematica and the GSL library work.
That being said I am not sure how important this is and my default urge is to
just follow the Matlab way of having a set of good solvers instead of an
exhaustive selection.

### Stiff Solvers
What I currently lack is a solver that implements the multistep methods of Adams
or Gear. These methods seem like they are much more complex, but essential in my
mind for having a good solver. My current feeling is that I will do a port from
the Sundials `CVODE` code. I think this should be doable as a lot of the
complexity of that code comes from having to move along all the LAPACK and
sparse matrix commands. We shall see.

## Future plans
I really will want to ultimately also support:
* BVP solvers
* DAE solvers
* DDE solvers
