# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# T is the type of the coefficients
abstract Tableau{Name, S, T <: Real}
Base.eltype{N, S, T}(b::Tableau{N, S, T}) = T
order(b::Tableau) = b.order
# Subtypes need to define a convert method to convert to a different
# eltype with signature:
Base.convert{Tnew <: Real}(::Type{Tnew}, tab::Tableau) = error("Define convert method for concrete Tableau types")

## Solver Output Types
abstract AbstractOdeSolution

#TODO: add solver information. Look at what is going on in `Optim.jl`
type RKOdeSolution <: AbstractOdeSolution
    # I currently hard code the types as this is the most common case. I can
    # look into generalizing this, but I will be wary of making the code overly
    # complex to do so
    x::Array{Float64, 1}
    y::Array{Float64, 2}
end

############################
# OdeProblem Types
############################
## Design
# the idea behind this setup is to have a bundle of memory and the functions
# needed for the ode solvers. This will also be used as the type dispatch
# for the generic ode solver interfaces (aode, dode, iode, etc).
#
abstract AbstractOdeSystem
abstract RungeKuttaSystem <: AbstractOdeSystem

##TODO: is it worth having a RungeKuttaOdeProblem <: AbstractOdeProblem abstract
## type to bundle the RK types?
type Dopri5 <: RungeKuttaSystem
    ndim::Int
    func::Function
    y0::Array{Float64, 1}
    # maybe this should include the initial_conditions y0? That is really a part of the problem, and it could be changed by sys.y0 = <blah> if needed
    ks::Array{Float64, 2}
    ywork::Array{Float64, 1} # this maybe should be called ystart/yinit
    ytrial::Array{Float64, 1}
    yerr::Array{Float64, 1}
    ytmp::Array{Float64, 1}
end

function Dopri5(func::Function, y0::Array{Float64, 1})
    #I have hard coded the stages into this `7` I think this makes the most
    #sense as each RK type will need to have its own constructor like this,
    #so a parametric type isn't needed.
    ndim = length(y0)
    Dopri5(
        ndim, # ndim
        func, # dydt
        y0, # y0
        Array(Float64, 7, ndim), #ks
        Array(Float64, ndim), #ywork
        Array(Float64, ndim), #ytrial
        Array(Float64, ndim), #yerr
        Array(Float64, ndim) #ytmp
    )
end
