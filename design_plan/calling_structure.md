## How the Code Should be Structured and function
Something that I have been thinking about is how I want this code to be
structured as far as call structures.

For example in Matlab you have a variety of call styles for the `tspan`
argument.

    1) giving it a 2-element array `[tstart tend]`
    2) giving it explicit points in array `trange`

For the first case it is not clear to me is what points are returned if the
2-elemnt version are given. Does it return the underlying steps taken? Upon
further thinking version 1) is silly to implement. What does it by you vs the
explicit version `linspace(tstart, tend, nsteps)` or similarly `tstart:stepsize:tend`
as both of these are `AbstractArray` types in Julia I have no need to support
the earlier Matlab syntax.

For a different look using `NDSolve` in Mathematica returns an interpolating
function instead of an array of values. You are then responsible for calling
this function to get the points you need. I have found that using this default
can be very nice, but is more times than not a bit of a hassle. That being said
I do like the idea of a version that will return such an object for use.

In truth if the `[tstart, tend]` version is used I should really only return
the Mathematica like object as it is kind of crazy otherwise (why would I
possibly want and array of the natural adaptive steps?). Really want I want
from this form is the dense output.

Also from my days as a Fortran user I really miss the ability to just call the
the solver for a single step inside a loop. I wonder if Julia is efficient
enough to do something like this. (This would be the iterator approach)

## New API ideas
Currently we follow the `ODE.jl` versions which is really just a simplified
version of Matlab's api. One thing that might be worth doing is making a
settings interface. Matlab uses `odeset` to do this. For Julia we would want
to do this with a custom type.

### Settings
Things that I need:
* `reltol = 1.0e-5`,
* `abstol = 1.0e-8`,
* `minstep = abs(tspan[end] - tspan[1])/1e18`,
* `maxstep = abs(tspan[end] - tspan[1])/2.5`,
* `initstep = 0.0` or `initial_step = 0.0`?

### Universial calling function
I was thinking of using the name `desolve` like R's similar package. I also
thought of names like `dsolve` `ndsolve`. I am not sure if `desolve` is the best
as it kind of reads like we are "unsolving" something ... Modern SciPy uses
`ode` for there driver and object/method access to updating the parameters like
`ode(func, tspan).set_method("dopri5")` etc. I guess the nice thing about this
name is that it is similar to the Matlab `odeXX` with the integer codes removed
as this can call any number of them. That being said I would like to be able to
solve delay and DAE problems with the same driver, not just ODE's. In truth I
think if I am going to go for the central function with method keyword I will
use `dsolve` there is no reason to worry about thinking about it re Mathematica
since thinking that `dsolve` should be symbolic by default is simply not true
in Julia.

So what would the master function look like:

```jl
dsolve(model, y0, tspan) # but what would the default solver be, dopri5 I think
dsolve(model, y0, tspan; method = :dopri5) # symbol version which seems bad, see how `Optim.jl` is moving over to type dispatch
dsolve(model, y0, tspan; method = Dopri(5, 4))
dsolve(model, y0, tspan; method = ERK(5, 4))
dsolve(model, y0, tspan; method = ExplicitRungeKutta(5, 4))
```

vs

```jl
ode45(model, y0, tspan)
ode45dp(model, y0, tspan)
explicit_rungekutta{5, 4}(model, y0, tspan)
variable_bdf(model, y0, tspan)
```

etc

It is really hard for me to think about what is better. Clearly the specialty
functions will in general have much shorter names/calling. That being said if
I don't use the Matlab like names it is not clear what I would name everything,
as truly descriptive names (like the bdf) can feel very vague/non-descriptive.

So I think I will go with the `dsolve` function and using type dispatch. This
also has the nice behavior of making it more similar to `Optim.jl`. Now an
issue is how to deal with the different kinds of solvers, largely from the
RungeKutta family.

Can I do type dispatch on something like `method = RungeKutta(tableau)`? I think
I would need to have this be a parameterized type like `method = RungeKutta{tableau}`
but I am not sure this is possible.

## Desirable Behavior
Now I would like it that I could get the following different behavior. If I want
a specific set of values I get back to solution table `sol` type like I have
currently implemented. But ideally I would like to be able to do:
* A step at a time iterator version, where I can change the settings to the
  integrator, like how I would use a function like this in Fortran.
* I give a range `[tstart, tstop]` and I get back an object that has the solver
  stored points, but I use like the `InterpolatingFunction` from Mathematica.
  That is I have a Hermite interpolation object that I can get any value between
  `[tstart, tend]`.
How to split up my code so that this works well, and as efficiently as possible
will be an interesting challenge, even before I do any work looking at multistep
methods.

## API from `ODE.jl` Discussions
It seems that the Julia community has been having a lot of back and forth
trying to solve the interface issue in a general way with little traction. The
issue seems to be the best way to express:

```jl
f(t, y)
f(t, y, y')
M(t)f(t, y)
M(t, y)f(t, y)
```

As problem definitions for basic IVP problems to complex DAE problems. A key
part of this discussion is how to do this efficiently, but also with a
convenient syntax. One discussion was to use a `OdeProblem` type that could
allocate some of the memory and temporary steps used by the solvers for
efficiency. I need to look this over and see how much of this I care about.
Also there turns out to be a pure Julia `DASSL` implementation which also has
iterator support using a coroutine. If I use this I will want to change the
naming from camelCase to the current Julian way.

### OdeProblem Object
I wonder if it might work to have a type that can be called with a function
argument that would build a problem object then you could do delegation on
the type

```jl
dsolve(RungeKutta(func, y0), tout)
dsolve(dopri5(func), y0, tout)
```

I am not sure what information should be in the object, `y0`, `tout`? Also
how I would want to name them.

Also the point of this for efficiency is that the function and jacobian
would work inplace, so you would need a way to specify this

```jl
dsolve(dopri5!(func!), y0, tout)
```

my guess would be to have this form for each method type so that I would
then have the different call signature of

```jl
func(t, y, ywork)
```

With this framework (`OdeProblem`) I could allocate the temporary arrays for
the dopri5 solver so that repeated calls to the solver would not have to
reallocate. Though for small-med problems I am not sure how much that will
matter. It also makes the slightly strange issue that `dsolve` would actually
change things inplace. Also it would seem natural to store the solver
diagnostics in the `OdeProblem` object, and also in the returned `OdeSolution`
object which feels odd.

Where this might truly help is if I made the itertor version as I would not have
the high overhead for repeated, small calls.

Really if I use the `OdeProblem` as a kind of store of the memory needed for
the solver to run then what I am doing with the solver function is acting like
a sort of settings object. That is

```jl
dsolve(dopri5(func), y0, tout; abstol = 1e-5)
```

etc is really just using the solver call as a way to set the options (hence why
it might make sense for `yo` and `tout` to be part of the `OdeProblem` though
I am still not sure).

Also if I go this route I am thinking I will want separate functions for if I
want an array of points back, a dense interpolating function like object, or
and iterator. Though I am not sure what names for this would be.
* aode (array ode)
* dode (dense ode)
* iode (itertor ode)
these names seem dangerously terse. But full spelling it out is clearly too
long. That being said I kind of like them.

```jl
aode(dopri5(func), y0, linspace(0, 100, 50); abstol = 1e-5)
aode(cvode(func), y0 tout; max_order = 10)
```

and then the code for `aode` would do dispatch on the type of the first
argument

```jl
aode(prob::dopri5, y0::Vector, tout::Vector; kwargs{for dopri5}...) = rk_adapt(prob, y0, tout, btab::RKExplicitTableau = bt_dopri5)
```

This is actually looking really good. But this makes me notice that really I
will want the problem types to be CapsCase. `Dopri5`, `VODE` (don't really want
the C as it is not the C code ... as I will port this to pure Julia).

Now if I use the `xode` name then I really won't want to support DDE as these
are not ODE's as far as I understand, whereas DAE's are.

### OdeProblem -> OdeSystem
I have decided to call this an `OdeSystem` as this gives the nice short variable
name `sys` or `osys` vs `prob` which I don't really.

Now currently all the fields are top level but I am toying with the idea of
adding the `y0` initial value to the system, thereby solving the `ndim` issue as
well as being a bit closer to containing all the problem/system information.
With this in mind I see that doing things like `sys.y0 = [blah]` could be common
what I am scared of is if changing the system fields like this is common that
the user might much around with the work arrays. Though this should largely
be harmless it seems ugly to me. maybe I want a subtype `Workspace` that
contains the memory stuff so the user would need to do `sys.workspace.ks` etc
which would make it clear that this is not a top level.
