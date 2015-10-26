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
abstract AbstractODESolution

#TODO: add solver information. Look at what is going on in `Optim.jl`
type RKODESolution <: AbstractODESolution
    # I currently hard code the types as this is the most common case. I can
    # look into generalizing this, but I will be wary of making the code overly
    # complex to do so
    x::Array{Float64, 1}
    y::Array{Float64, 2}
end

#################################################
# Settings/Options Type
#################################################
type IVPOptions
    reltol::Float64 #TODO this can also be an array
    abstol::Float64 #TODO this can also be an array
    minstep::Float64
    maxstep::Float64
    initstep::Float64 #TODO or initial_step?
    method::Symbol
    #norm_func #Do I want this?
    #max_num_steps
end

##TODO: SO I have an issue in that min/max{step} requires access to the time
## range of the call to generate the actual Float64 numbers.
ivpoptions(;
    abstol = 1e-3,
    reltol = 1e-6,
    minstep = :auto,
    maxstep = :auto
    initstep = 0.0,
    method = :dopri45
) = IVPOptions(reltol, abstol, minstep, maxstep, initstep, method)
