## Next Steps

### Callbacks
So currently we have implemented a type system that is similar to how the
fortran and C codes pass around memory for the solvers, currently we call this
an `ODESystem`.

Looking at https://github.com/JuliaLang/ODE.jl/pull/33 and
https://github.com/JuliaLang/ODE.jl/issues/30 a related issue is how to think
about how the callback function is used internally. They key is to have the
ability to change the function to work inplace instead of having:
```jl
dydt = f(t, y)
```
to have instead:
```jl
f!(t, y, dydt)
```
That way we don't incur the memory allocation on each call. The way that the
above suggests to to have a special type for the function and then make a high
level wrapper that takes a regular ```dydt = f(t, y)``` and place it into a
dummy version that does the call inplace. This way the user has the option to
make their own efficient version if they choose. This is similar to how I think
the `Optim.jl` package uses `DifferentiableFunction` and
`TwiceDifferentiableFunction`. Looking at the way the `Optim.jl` uses this to
also allow for efficient calculation of both the function and the jacobian in
one function call will be good to think about for when we try to implement
stiff solvers.

So in summary it is worth looking into having as part of the `ODESystem` type
some kind of `ODEFunction` or some such type to wrap the callback. Then have
all internal calls use the inplace version, dummy or not.

In the issue 30 discussion there is a more subtle issue of how to model the
form of the ODE in general citing the package `PETSc` which uses the model
```
F(t, u, u') = G(t, u)
```
for the general case. This seems interesting, but I have so little experience
with this kind of issue I am not sure how to deal with it.

### Runge-Kutta Methods
Currently we closely follow the structure of the `ODE.jl` packages structure
which is extremely generic, allowing for inputs of any Butcher tableau with
a generic stepper, one for fixed step and one for adaptive. We have removed
the fixed steppers as it is of very little utility outside of very specific
use cases. In generally it increasingly feels like there is little practical
benefit for this generality. From the treatise by Hairer and Wanner (1993) it
is clear that `Dopri5/Dop853` are a clear winners among the explicit Runge-Kutta
solvers, crushing `RKFehlberg45` codes. The only exception might be the more
recent `ode23` like codes from Shampine. Whatever the case we will move towards
only supporting these 2-3 solvers for the explicit Runge-Kutta family. To this
end some nice simplifications can be made:

* Simplify the Butcher tableau type
    * The FSAL property should be a boolean field that is set at construction
    * The order of the tableau should be a field set at construction and not a
      type parameter

Once these changes are implemented the question comes on how to call the `Dopri`
codes. Currently we have specific types for the specific type in `Dopri54`.
No since we want to support both `Dopri54` and `Dopri853` the question becomes
whether it makes sense to use the current method or instead have a general
type `Dopri` that then has an argument for the order `Dopri(func, order = 54)`.
And then for the regular case we just default to the best default solver
(`Dopri853` I believe). The only issue with this is how to do dispatch in this
case as the solver method might need the order information for things like
the Hermite interpolation and how the step embedding works. The only solution
for this would be to either have manual dispatch with `if` statements. Or have
the order be type instances that can be delegated on, like
`Dopri(func, order = dp54)` where `dp45 <: RKOrder` or some such.

### Missing Feature from the RK Codes
The code from `ODE.jl` is clearly based on the Hairer and Wanner monograph, but
it is not clear at times how it has been changed. Things that need to be
investigated:

* Some of the code has some magical constants like `timeout_const = 5` which
  feel unmotivated in the current form.
* It is not clear that the stiffness check has been ported over.

Features that are missing that are known include:
* the Hermite interpolation is of the incorrect order
    * Dopri5 should have 4th order
    * Dop853 should have 7th order
  and really there are special forms derived by Shampine for the forms of these
  Hermite coefficients that take advantage of the ode formulation. These
  should simply be ported over from the original codes.
* Event detection.

## Benchmarks
It would be great to have a set of benchmarks that output figures like in
Hairer and Wanner showing the tolerance vs the number of function calls, also
we could add tolerance vs cpu time.

It would also be good to see this for increasing the size of the problem using
a general code for linear systems. That way we could see how the solvers
perform for growing problem size.
