import Base.call

## This file implements a basic dense output type using the Hermite
## interpolation that is related to the underlying solver interface of things
## like the DopriX methods.
##
## This work will try as best to follow convenctions set out in the packages
## `Interpolation.jl` and its successor `Grid.jl`. Hopefully with time these
## packages can be used inplace, or as the underlying technology of these types.
type DenseODESolution
    ndim::Int
    order::Int # Interpolation order
    xvals::Array{Float64, 1}
    yvals::Array{Float64, 2}
    fvals::Array{Any, 1} # hermite coefficients
end

function DenseODESolution(xvals::Array{Float64, 1}, yvals::Array{Float64, 2}, fvals::Array{Any, 1})
    issorted(xvals) || throw(ArgumentError("x points must be in ascending order"))
    # I am only using 3-order interpolation at the moment, might want to be able to set this in the future
    DenseODESolution(size(yvals, 2), 3, xvals, yvals, fvals)
end

function (sol::DenseODESolution)(xval::Float64)
    ((xval < sol.xvals[1]) | (xval > sol.xvals[end])) && throw(ArgumentError("Outside of Range"))
    xr = searchsorted(sol.xvals, xval)
    i1 = xr.stop
    i2 = xr.start
    # if both are the same we are on a knot point, just return it
    i2 - i1 == 0 && return sol.yvals[:, i1]
    # otherwise carry out the interpolation
    dt = sol.xvals[i2] - sol.xvals[i1]
    yout = hermite_shampine_interp(xval, sol.xvals[i1], dt, sol.fvals[i1])
    return yout
end

(sol::DenseODESolution)(xvals::AbstractArray{Float64}) = vcat([sol(xval) for xval in xvals]...) # return a table
