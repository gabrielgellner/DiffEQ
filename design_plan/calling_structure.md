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
* `initstep = 0.0`
