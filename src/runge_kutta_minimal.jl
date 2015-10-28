# Explicit Runge-Kutta solvers
##############################
# (Hairer & Wanner 1992 p.134, p.165-169)

##############################
# Adaptive Runge-Kutta methods
##############################
#NOTE: naming convenction bt and btab are shorthand for Butcher Tableaus
##TODO: get rid of the kwargs... and be explicit
aode(sys::Dopri5, tspan; kwargs...) = oderk_adapt(sys, tspan, bt_dopri5; kwargs...)

function oderk_adapt{N, S}(
                     sys::RungeKuttaSystem,
                     tspan::AbstractVector{Float64},
                     btab::TableauRKExplicit{N, S};
                     reltol = 1.0e-5,
                     abstol = 1.0e-8,
                     minstep = abs(tspan[end] - tspan[1])/1e18,
                     maxstep = abs(tspan[end] - tspan[1])/2.5,
                     initstep = 0.0
                     )
    # parameters
    order = minimum(btab.order)
    const timeout_const = 5 # after step reduction do not increase step for
                            # timeout_const steps

    ## Initialization
    #ndim = length(y0)
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]

    # initialize work arrays
    sys.ywork = copy(sys.y0)

    # output ys
    nsteps_fixed = length(tspan)
    ys = Array(Float64, nsteps_fixed, sys.ndim)
    ys[1, :] = sys.y0

    # Time
    dt, tdir, sys.ks[1, :] = hinit(sys.func, sys.ywork, tstart, tend, order, reltol, abstol) # sets ks[1, :] = f0
    if initstep != 0
        dt = sign(initstep) == tdir ? initstep : error("initstep has wrong sign.")
    end
    # Diagnostics
    dts = Float64[]
    errs = Float64[]
    steps = [0, 0]  # [accepted, rejected]

    ## Integration loop
    islaststep = abs(t + dt - tend) <= eps(tend) ? true : false
    timeout = 0 # for step-control
    iter = 2 # the index into tspan and ys
    while true
        # do one step (assumes ks[1, :] == f0)
        rk_embedded_step!(sys.ytrial, sys.yerr, sys.ks, sys.ytmp, sys.ywork, sys.func, t, dt, sys.ndim, btab)
        # Check error and find a new step size:
        ##TODO: we call the function with y, ytrial, yerr names, but the function
        ## uses names like x0, xtrial, xerr which is confusing
        err, newdt, timeout = stepsize_hw92!(dt, tdir, sys.ywork, sys.ytrial, sys.yerr, order, timeout,
                                             sys.ndim, abstol, reltol, maxstep, norm)

        if err <= 1.0 # accept step
            # diagnostics
            steps[1] += 1
            ##TODO: as I have removed variable output size I likely don't need
            ## to grow these
            push!(dts, dt)
            push!(errs, err)

            # Output:
            ##NOTE in `ODE.jl` as they are using array of arrays which means the line
            ## f0 = ks[1] is a view not a copy
            f0 = sub(sys.ks, 1, :)
            ##FSAL -> First Same As Last, this code seems like it could be set in a btab field instead
            ## of always being recalculated
            f1 = isFSAL(btab) ? sub(sys.ks, S, :) : fn(t + dt, sys.ytrial)
            # interpolate onto given output points
            while iter - 1 < nsteps_fixed && (tdir*tspan[iter] < tdir*(t + dt) || islaststep) # output at all new times which are < t+dt
                # TODO: 3rd order only!
                hermite_interp!(sub(ys, iter, :), tspan[iter], t, dt, sys.ywork, sys.ytrial, f0, f1)
                iter += 1
            end
            sys.ks[1, :] = f1 # load ks[1, :] == f0 for next step

            # Break if this was the last step:
            islaststep && break

            # Swap bindings of y and ytrial, avoids one copy
            sys.ywork, sys.ytrial = sys.ytrial, sys.ywork

            # Update t to the time at the end of current step:
            t += dt
            dt = newdt

            # Hit end point exactly if next step within 1% of end:
            if tdir*(t + dt*1.01) >= tdir*tend
                dt = tend - t
                islaststep = true # next step is the last, if it succeeds
            end
        elseif abs(newdt) < minstep  # minimum step size reached, break
            println("Warning: dt < minstep.  Stopping.")
            break
        else # redo step with smaller dt
            islaststep = false
            steps[2] += 1
            dt = newdt
            timeout = timeout_const
        end
    end
    return RKOdeSolution(tspan, ys)
end

# estimator for initial step based on book
# "Solving Ordinary Differential Equations I" by Hairer et al., p.169
function hinit(fn, y0, tstart, tend, order, reltol, abstol)
    # Returns first step, direction of integration and F evaluated at t0
    tdir = sign(tend - tstart)
    tdir == 0 && error("Zero time span")
    tau = max(reltol*norm(y0, Inf), abstol)
    d0 = norm(y0, Inf)/tau
    f0 = fn(tstart, y0)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    # perform Euler step
    x1 = y0 + tdir*h0*f0
    f1 = fn(tstart + tdir*h0, x1)
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

function rk_embedded_step!{N, S}(ytrial, yerr, ks, ytmp, y, fn, t, dt, ndim, btab::TableauRKExplicit{N, S})
    # Does one embedded R-K step updating ytrial, yerr and ks.
    #
    # Assumes that ks[:, 1] is already calculated!
    #
    # Modifies ytrial, yerr, ks, and ytmp
    ##NOTE: currently hard coded to be Float64
    fill!(ytrial, 0.0)
    fill!(yerr, 0.0)
    for d = 1:ndim
        ytrial[d] += btab.b[1, 1]*ks[1, d]
        yerr[d] += btab.b[2 ,1]*ks[1, d]
    end
    for s = 2:S
        calc_next_k!(ks, ytmp, y, s, fn, t, dt, ndim, btab)
        for d = 1:ndim
            ytrial[d] += btab.b[1, s]*ks[s, d]
            yerr[d] += btab.b[2, s]*ks[s, d]
        end
    end
    for d = 1:ndim
        yerr[d] = dt*(ytrial[d] - yerr[d])
        ytrial[d] = y[d] + dt*ytrial[d]
    end
end

function stepsize_hw92!(dt, tdir, x0, xtrial, xerr, order,
                       timeout, ndim, abstol, reltol, maxstep, norm)
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
    for d = 1:ndim
        # if outside of domain (usually NaN) then make step size smaller by maximum
        isoutofdomain(xtrial[d]) && return 10.0, dt*facmin, timout_after_nan
        xerr[d] = xerr[d]/(abstol + max(norm(x0[d]), norm(xtrial[d]))*reltol) # Eq 4.10
    end
    err = norm(xerr, 2) # Eq. 4.11
    newdt = min(maxstep, tdir*dt*max(facmin, fac*(1/err)^(1/(order + 1)))) # Eq 4.13 modified
    if timeout > 0
        newdt = min(newdt, dt)
        timeout -= 1
    end
    return err, tdir*newdt, timeout
end

function calc_next_k!(ks::Matrix{Float64}, ytmp::Vector{Float64}, y, s, fn, t, dt, ndim, btab)
    # Calculates the next ks and puts it into ks[s, :]
    # - ks and ytmp are modified inside this function.
    ytmp[:] = y
    for ss = 1:(s - 1), d = 1:ndim
        ytmp[d] += dt*ks[ss, d]*btab.a[s, ss]
    end
    ks[s, :] = fn(t + btab.c[s]*dt, ytmp)
    nothing
end
