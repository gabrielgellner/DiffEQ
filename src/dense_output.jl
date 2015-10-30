## This file implements a basic dense output type using the Hermite
## interpolation that is related to the underlying solver interface of things
## like the DopriX methods.
##
## This work will try as best to follow convenctions set out in the packages
## `Interpolation.jl` and its successor `Grid.jl`. Hopefully with time these
## packages can be used inplace, or as the underlying technology of these types.
type DenseOdeSolution
    order::Int # Interpolation order
    xvals::Array{Float64, 1}
    fvals::Array{Float64, 2}
    mesh::Array{Float64, 2}
end

getindex(sol::DenseOdeSolution, xval::Float64) = hermite_interp(xval, )
