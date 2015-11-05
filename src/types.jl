################################################################################
## Solver Output Types
################################################################################
abstract AbstractODESolution

#TODO: add solver information. Look at what is going on in `Optim.jl`
type RKODESolution <: AbstractODESolution
    # I currently hard code the types as this is the most common case. I can
    # look into generalizing this, but I will be wary of making the code overly
    # complex to do so
    x::Array{Float64, 1}
    y::Array{Float64, 2}
end

################################################################################
## ODESystem Types
################################################################################
## Design
# the idea behind this setup is to have a bundle of memory and the functions
# needed for the ode solvers. This will also be used as the type dispatch
# for the generic ode solver interfaces (aode, dode, iode, etc).
#
abstract AbstractODESystem
abstract RungeKuttaSystem <: AbstractODESystem

type RKWorkspace
    ndim::Int
    ks::Array{Float64, 2}
    yinit::Array{Float64, 1}
    ytrial::Array{Float64, 1}
    yerr::Array{Float64, 1}
    ytmp::Array{Float64, 1}
    out_i::Int # used for fixed size output ##TODO think of a better way
end

##TODO: I am not entirely sold on this name. Clearly it is a more explicit version of
## Dopri5 which is reminicint of the fortran codes, and will be nice for when I
## implement Dopri853. But would something like RKDP54 be better? I don't like
## how many capitals and jargony that feels. I coudl do something like
## ODE54 ... but that has capitals and is vague. Need to think about this
type Dopri54 <: RungeKuttaSystem
    ##TODO: I shouln't let ndim change
    func::Function
    y0::Array{Float64, 1}
    work::RKWorkspace
end

function Dopri54(func::Function, y0::Array{Float64, 1})
    #I have hard coded the stages into this `7` I think this makes the most
    #sense as each RK type will need to have its own constructor like this,
    #so a parametric type isn't needed.
    ndim = length(y0)
    Dopri54(
        func, # dydt
        y0, # y0
        RKWorkspace(
            ndim, #ndim
            Array(Float64, ndim, 7), #ks
            Array(Float64, ndim), #ywork
            Array(Float64, ndim), #ytrial
            Array(Float64, ndim), #yerr
            Array(Float64, ndim), #ytmp
            0 # out_i
        )
    )
end
