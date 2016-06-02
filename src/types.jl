################################################################################
## Solver Output Types
################################################################################
abstract AbstractODESolution
abstract AbstractODESolutionStatistics

type ODESolutionStatistics <: AbstractODESolutionStatistics
    nfevals::Int
    nsteps::Int
    nfailed::Int
end

#TODO: add solver information. Look at what is going on in `Optim.jl`
type RKODESolution <: AbstractODESolution
    # I currently hard code the types as this is the most common case. I can
    # look into generalizing this, but I will be wary of making the code overly
    # complex to do so
    x::Array{Float64, 1}
    y::Array{Float64, 2}
    stats::ODESolutionStatistics
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
    ydim::Int
    ks::Array{Float64, 2}
    yinit::Array{Float64, 1}
    ytrial::Array{Float64, 1}
    yerr::Array{Float64, 1}
    ytmp::Array{Float64, 1}
    ycont::Array{Float64, 2} # used for dense/continous output
    tstart::Float64
    tend::Float64
    dt::Float64
    tdir::Int # direction of the time steps (should be 1 or -1)
    out_i::Int # used for fixed size output ##TODO think of a better way
    order::Int # order of the runge-kutta method
    laststep::Bool
end

type RKOptions
    reltol::Float64
    abstol::Float64
    maxstep::Float64
    minstep::Float64
    initstep::Float64
end

function RKOptions(;
    reltol = 1.0e-5,
    abstol = 1.0e-8,
    maxstep = 5.0, # dopri5.f suggests between 1.5 - 5.0
    minstep = 1.0/maxstep,
    initstep = 0.0
    )
    RKOptions(reltol, abstol, maxstep, minstep, initstep)
end

##TODO: I am not entirely sold on this name. Clearly it is a more explicit version of
## Dopri5 which is reminicint of the fortran codes, and will be nice for when I
## implement Dopri853. But would something like RKDP54 be better? I don't like
## how many capitals and jargony that feels. I coudl do something like
## ODE54 ... but that has capitals and is vague. Need to think about this
type Dopri54 <: RungeKuttaSystem
    ##TODO: I shouln't let ndim change
    func::Function
    y0::Array{Float64, 1} ##TODO: allow any AbstractFloat
    work::RKWorkspace
end

# the dopri5.f code seems to use the maximum not the minimum (for the order) -- whereas
#`ODE.jl` uses the formulas from the book which use the minimum. This needs to be resolved.
# It seems to be an issue when using embedded methods whether to use the larger or smaller
# order method for extrapolation, though the dormand prince pairs where specifically
# designed for the larger pairs being used for extrapolation to be less problamatic
# (as described in the Butcher 2008 book)
#TODO: allow any AbstractFloat
function Dopri54(func::Function, y0::Array{Float64, 1})
    #I have hard coded the stages into this `7` I think this makes the most
    #sense as each RK type will need to have its own constructor like this,
    #so a parametric type isn't needed.
    ydim = length(y0)
    Dopri54(
        func, # dydt
        y0, # y0
        RKWorkspace(
            ydim, #ndim
            Array{Float64}(ydim, 7), #ks
            Array{Float64}(ydim), #ywork
            Array{Float64}(ydim), #ytrial
            Array{Float64}(ydim), #yerr
            Array{Float64}(ydim), #ytmp
            Array{Float64}(ydim, 5), #ycont 5 comes from the fact that it is 5th order interpolant
            0.0, # tstart
            0.0, # tend
            0, # tdir
            0.0, # dt
            0, #out_i
            5, #order
            false, #laststep
        )
    )
end
