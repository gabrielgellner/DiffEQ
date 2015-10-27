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

Also from my days as a Fortran user I really miss the ability to just call the
the solver for a single step inside a loop. I wonder if Julia is efficient
enough to do something like this.

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

dsolve(model, y0, tspan) # but what would the default solver be, dopri5 I think
dsolve(model, y0, tspan; method = :dopri5) # symbol version which seems bad, see how `Optim.jl` is moving over to type dispatch
dsolve(model, y0, tspan; method = Dopri(5, 4))
dsolve(model, y0, tspan; method = ERK(5, 4))
dsolve(model, y0, tspan; method = ExplicitRungeKutta(5, 4))

vs

ode45(model, y0, tspan)
ode45dp(model, y0, tspan)
explicit_rungekutta{5, 4}(model, y0, tspan)
variable_bdf(model, y0, tspan)

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
I would need to have this be a parametrized type like `method = RungeKutta{tableau}`
but I am not sure this is possible.
