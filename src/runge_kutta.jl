# Explicit Runge-Kutta solvers
##############################
# (Hairer & Wanner 1992 p.134, p.165-169)

##############################
# Adaptive Runge-Kutta methods
##############################
#NOTE: naming convenction bt and btab are shorthand for Butcher Tableaus
##TODO: get rid of the kwargs... and be explicit
aode(sys::Dopri54, tspan; kwargs...) = rksolver_array(sys, tspan, bt_dopri54; kwargs...)
dode(sys::Dopri54, tspan; kwargs...) = rksolver_dense(sys, tspan, bt_dopri54; kwargs...)

function rksolver_array(sys::RungeKuttaSystem,
                        tspan::AbstractVector{Float64},
                        btab::TableauRKExplicit;
                        reltol = 1.0e-5,
                        abstol = 1.0e-8,
                        minstep = abs(tspan[end] - tspan[1])/1e18,
                        maxstep = abs(tspan[end] - tspan[1])/2.5,
                        initstep = 0.0
                        )
    # parameters
    order = minimum(btab.order) #TODO: why is the order the minimum?

    ## Initialization
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]

    # initialize work arrays
    sys.work.yinit = copy(sys.y0)

    # output ys
    nsteps_fixed = length(tspan)
    ##Note: it is more column major to think of an array of points joined along
    ## columns. When returned it must be transposed.
    ys = Array(Float64, sys.work.ndim, nsteps_fixed)
    ys[:, 1] = sys.y0

    # Time
    dt, tdir, sys.work.ks[:, 1] = hinit(sys, tstart, tend, order, reltol, abstol) # sets ks[:, 1] = f0
    if initstep != 0
        dt = sign(initstep) == tdir ? initstep : error("initstep has wrong sign.")
    end

    ## Integration loop
    dts, errs, steps = rk_stepper!(sys, t, dt, tdir, tend, tspan, ys, [], btab, order, abstol, reltol, minstep, maxstep, norm, rk_array_output!)

    return RKODESolution(tspan, ys')
end

function rksolver_dense(sys::RungeKuttaSystem,
                        tspan::AbstractVector{Float64},
                        btab::TableauRKExplicit;
                        reltol = 1.0e-5,
                        abstol = 1.0e-8,
                        minstep = abs(tspan[end] - tspan[1])/1e18,
                        maxstep = abs(tspan[end] - tspan[1])/2.5,
                        initstep = 0.0
                        )
    # parameters
    order = minimum(btab.order) # it might be worth adding this as a field to the btab

    ## Initialization
    if length(tspan) > 2
        error("Dense output requires tspan to be two points (tstart, tend).")
    end
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]
    tout = [copy(tstart)]

    # initialize work arrays
    sys.work.yinit = copy(sys.y0)

    # ys is an array of arrays so that it can grow as needed, this will be
    # converted to an array (table: nsteps x ndim array) at output
    ys = Array{Float64, 1}[copy(sys.y0)]
    fs = Array{Float64, 1}[sys.func(tstart, sys.y0)]

    # Time
    dt, tdir, sys.work.ks[:, 1] = hinit(sys, tstart, tend, order, reltol, abstol) # sets ks[:, 1] = func(y0)
    if initstep != 0
        dt = sign(initstep) == tdir ? initstep : error("initstep has wrong sign.")
    end

    # Integrate
    ##TODO: work on this argument list!
    dts, errs, steps = rk_stepper!(sys, t, dt, tdir, tend, tout, ys, fs, btab, order, abstol, reltol, minstep, maxstep, norm, rk_dense_output!)

    # Output solution
    return DenseODESolution(tout, hcat(ys...)', hcat(fs...)')
end

# estimator for initial step based on book
# "Solving Ordinary Differential Equations I" by Hairer et al., p.169
function hinit(sys::AbstractODESystem, tstart, tend, order, reltol, abstol)
    # Returns first step, direction of integration and F evaluated at t0
    tdir = sign(tend - tstart)
    tdir == 0 && error("Zero time span")
    tau = max(reltol*norm(sys.work.yinit, Inf), abstol)
    d0 = norm(sys.work.yinit, Inf)/tau
    f0 = sys.func(tstart, sys.work.yinit)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    # perform Euler step
    x1 = sys.work.yinit + tdir*h0*f0
    f1 = sys.func(tstart + tdir*h0, x1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 1e-15
        h1 = max(1e-6, 1e-3*h0)
    else
        pow = -(2.0 + log10(max(d1, d2)))/(order + 1.0)
        h1 = 10.0^pow
    end
    return tdir*min(100*h0, h1, tdir*(tend - tstart)), tdir, f0
end

##TODO: it might make more sense just to pass the {N, S} instead of using type generic for it
function rk_stepper!(sys::RungeKuttaSystem, t, dt, tdir, tend, tout, ys, fs, btab::TableauRKExplicit, order, abstol, reltol, minstep, maxstep, norm, output_func!::Function)
    # after step reduction do not increase step for `timeout_const` steps
    const timeout_const = 5

    # Diagnostics
    dts = Float64[]
    errs = Float64[]
    steps = [0, 0]  # [accepted, rejected]

    ## Integration loop
    nout = size(ys, 2) # needed for array_output
    islaststep = abs(t + dt - tend) <= eps(tend) ? true : false
    timeout = 0 # for step-control
    sys.work.out_i = 2 # the index into tout and ys
    while true
        # do one step (assumes ks[1, :] == f0)
        rk_embedded_step!(sys, t, dt, btab)
        # Check error and find a new step size:
        err, newdt, timeout = stepsize_hw92!(sys, dt, tdir, order, timeout, abstol, reltol, maxstep, norm)

        if err <= 1.0 # accept step
            # diagnostics
            steps[1] += 1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            ##NOTE in `ODE.jl` as they are using array of arrays which means the line
            ## f0 = ks[1] is a view not a copy
            f0 = sub(sys.work.ks, :, 1)
            # FSAL -> First Same As Last
            f1 = isfsal(btab) ? sub(sys.work.ks, :, btab.nstages) : fn(t + dt, sys.work.ytrial)
            # interpolate onto given output points

            # output/save step
            output_func!(sys, ys, fs, tout, t, dt, f0, f1, nout, islaststep)

            ####################################################################
            # prepare for next iteration
            ####################################################################
            sys.work.ks[:, 1] = f1 # load ks[:, 1] == f0 for next step

            # Break (end integration loop) if this was the last step:
            islaststep && break

            # Swap bindings of y and ytrial, avoids one copy
            sys.work.yinit, sys.work.ytrial = sys.work.ytrial, sys.work.yinit

            # Update t to the time at the end of current step:
            t += dt
            dt = newdt

            # Hit end point exactly if next step within 1% of end:
            if tdir*(t + dt*1.01) >= tdir*tend
                dt = tend - t
                islaststep = true # next step is the last, if it succeeds
            end
        elseif abs(newdt) < minstep  # minimum step size reached, break
            ##TODO; make this a real error, this is likely the "stiffness" check from the original code. Compare and give a meaningful error name if so
            println("Warning: dt < minstep.  Stopping.")
            break
        else # redo step with smaller dt
            islaststep = false
            steps[2] += 1
            dt = newdt
            timeout = timeout_const
        end
    end
    return dts, errs, steps
end

##TODO: I have broken the tdir variable, so I need to add checks for
## integrating in reverse
function rk_array_output!(sys, ys, fs, tspan, t, dt, f0, f1, nout, islaststep)
    # interpolate onto requested times in (t, t + dt)
    ##TODO: I am missing the last timestep
    ## we need out_i - 1 < nout so that we don't have infinite loop at laststep
    while sys.work.out_i - 1 < nout && (tspan[sys.work.out_i] < t + dt || islaststep)
        # TODO: 3rd order only!
        hermite_interp!(sub(ys, :, sys.work.out_i), tspan[sys.work.out_i], t, dt, sys.work.yinit, sys.work.ytrial, f0, f1)
        sys.work.out_i += 1
    end
end

##TODO: currently I have to pass in nout and islaststep because of the array_output
## find out how to avoid this
function rk_dense_output!(sys, ys, fs, tout, t, dt, f0, f1, nout, islaststep)
    push!(ys, deepcopy(sys.work.ytrial))
    push!(tout, t + dt)
    # Also save the derivatives for the hermite interpolation
    push!(fs, f1) # ys is the right point so take f1
end

function rk_embedded_step!(sys, t, dt, btab::TableauRKExplicit)
    # Does one embedded R-K step updating ytrial, yerr and ks.
    #
    # Assumes that work.ks[1, :] is already calculated!
    #
    # Modifies work.ytrial, work.yerr, work.ks, and work.ytmp
    ##NOTE: currently hard coded to be Float64
    fill!(sys.work.ytrial, 0.0)
    fill!(sys.work.yerr, 0.0)
    for d = 1:sys.work.ndim
        sys.work.ytrial[d] += btab.b[1, 1]*sys.work.ks[d, 1]
        sys.work.yerr[d] += btab.b[2 ,1]*sys.work.ks[d, 1]
    end
    for s = 2:btab.nstages
        calc_next_k!(sys, s, t, dt, btab)
        for d = 1:sys.work.ndim
            sys.work.ytrial[d] += btab.b[1, s]*sys.work.ks[d, s]
            sys.work.yerr[d] += btab.b[2, s]*sys.work.ks[d, s]
        end
    end
    for d = 1:sys.work.ndim
        sys.work.yerr[d] = dt*(sys.work.ytrial[d] - sys.work.yerr[d])
        sys.work.ytrial[d] = sys.work.yinit[d] + dt*sys.work.ytrial[d]
    end
end

function stepsize_hw92!(sys, dt, tdir, order, timeout, abstol, reltol, maxstep, norm)
    # Estimates the error and a new step size following Hairer &
    # Wanner 1992, p167 (with some modifications)
    #
    # If timeout > 0 no step size increase is allowed, timeout is
    # decremented in here.
    #
    # Returns the error, newdt and the number of timeout-steps
    #
    # TODO:
    # - allow component-wise reltol and abstol?
    timout_after_nan = 5
    ##TODO: this is a very complicated way to calculate 0.8 ;)
    fac = [0.8, 0.9, 0.25^(1/(order + 1)), 0.38^(1/(order + 1))][1]
    facmax = 5.0 # maximal step size increase. 1.5-5
    facmin = 1.0/facmax  # maximal step size decrease. ?

    # in-place calculate xerr./tol
    for d = 1:sys.work.ndim
        # if outside of domain (usually NaN) then make step size smaller by maximum
        isoutofdomain(sys.work.ytrial[d]) && return 10.0, dt*facmin, timout_after_nan
        sys.work.yerr[d] = sys.work.yerr[d]/(abstol + max(norm(sys.work.yinit[d]), norm(sys.work.ytrial[d]))*reltol) # Eq 4.10
    end
    err = norm(sys.work.yerr, 2) # Eq. 4.11
    newdt = min(maxstep, tdir*dt*max(facmin, fac*(1/err)^(1/(order + 1)))) # Eq 4.13 modified
    if timeout > 0
        newdt = min(newdt, dt)
        timeout -= 1
    end
    return err, tdir*newdt, timeout
end

function calc_next_k!(sys, s, t, dt, btab)
    # Calculates the next ks and puts it into ks[s, :]
    # - ks and ytmp are modified inside this function.
    sys.work.ytmp[:] = sys.work.yinit
    for ss = 1:(s - 1), d = 1:sys.work.ndim
        sys.work.ytmp[d] += dt*sys.work.ks[d, ss]*btab.a[s, ss]
    end
    sys.work.ks[:, s] = sys.func(t + btab.c[s]*dt, sys.work.ytmp)
    nothing
end
