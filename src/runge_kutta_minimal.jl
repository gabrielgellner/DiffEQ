# Explicit Runge-Kutta solvers
##############################
# (Hairer & Wanner 1992 p.134, p.165-169)

##############################
# Adaptive Runge-Kutta methods
##############################
#NOTE: naming convenction bt and btab are shorthand for Butcher Tableaus

## High lever interface, this will change
#ode21(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_rk21; kwargs...)
#ode23(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_rk23; kwargs...)
#ode45_fe(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_rk45; kwargs...)
ode45_dp(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_dopri5; kwargs...)
# Use Dormand-Prince version of ode45 by default
const ode45 = ode45_dp
#ode78(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_feh78; kwargs...)

# This function is the meat of the adaptive solvers.
function oderk_adapt{N, S}(fn, y0::Vector{Float64}, tspan::Vector{Float64},
                     btab::TableauRKExplicit{N, S};
                     reltol = 1.0e-5,
                     abstol = 1.0e-8,
                     minstep = abs(tspan[end] - tspan[1])/1e18,
                     maxstep = abs(tspan[end] - tspan[1])/2.5,
                     initstep = 0.0
                     )

    ##TODO: think about how to get rid of this check. Why not have two types of
    # tableaus? So I can just have the dispatch deal with this?
    !isadaptive(btab) && error("Can only use this solver with an adaptive RK Butcher table")

    #TODO: get rid of this call, I want to just have the tableaus be concrete
    # types of Float64
    #btab = convert(Float64, btab_)

    # parameters
    order = minimum(btab.order)
    timeout_const = 5 # after step reduction do not increase step for
                      # timeout_const steps

    ## Initialization
    dof = length(y0)
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]

    # work arrays:
    y = copy(y0)      # y at time t
    ytrial = Array(Float64, dof) # trial solution at time t+dt
    yerr = Array(Float64, dof) # error of trial solution
    ks = Array(Float64, S, dof)
    # allocate!(ks, y0, dof) # no need to allocate as fn is not in-place
    ytmp = Array(Float64, dof)

    # output ys
    nsteps_fixed = length(tspan) # these are always output
    ys = Array(Float64, nsteps_fixed, dof)
    #allocate!(ys, y0, dof)
    ys[1, :] = y0

    # Time
    dt, tdir, ks[1, :] = hinit(fn, y, tstart, tend, order, reltol, abstol) # sets ks[1] = f0
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
        # do one step (assumes ks[1] == f0)
        rk_embedded_step!(ytrial, yerr, ks, ytmp, y, fn, t, dt, dof, btab)
        # Check error and find a new step size:
        err, newdt, timeout = stepsize_hw92!(dt, tdir, y, ytrial, yerr, order, timeout,
                                             dof, abstol, reltol, maxstep, norm)

        if err <= 1.0 # accept step
            # diagnostics
            steps[1] += 1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            f0 = ks[1, :]
            f1 = isFSAL(btab) ? ks[S, :] : fn(t + dt, ytrial)
            # interpolate onto given output points
            while iter - 1 < nsteps_fixed && (tdir*tspan[iter] < tdir*(t + dt) || islaststep) # output at all new times which are < t+dt
                ys[iter, :] = hermite_interp(tspan[iter], t, dt, y, ytrial, f0, f1) # TODO: 3rd order only!
                iter += 1
            end
            ks[1, :] = f1 # load ks[1] == f0 for next step

            # Break if this was the last step:
            islaststep && break

            # Swap bindings of y and ytrial, avoids one copy
            y, ytrial = ytrial, y

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
    return RKODESolution(tspan, ys)
end

function rk_embedded_step!{N, S}(ytrial, yerr, ks, ytmp, y, fn, t, dt, dof, btab::TableauRKExplicit{N, S})
    # Does one embedded R-K step updating ytrial, yerr and ks.
    #
    # Assumes that ks[:, 1] is already calculated!
    #
    # Modifies ytrial, yerr, ks, and ytmp
    fill!(ytrial, zero(eltype(ytrial)))
    fill!(yerr, zero(eltype(ytrial)))
    for d = 1:dof
        ytrial[d] += btab.b[1, 1]*ks[1, d]
        yerr[d] += btab.b[2 ,1]*ks[1, d]
    end
    for s = 2:S
        calc_next_k!(ks, ytmp, y, s, fn, t, dt, dof, btab)
        for d = 1:dof
            ytrial[d] += btab.b[1, s]*ks[s, d]
            yerr[d] += btab.b[2, s]*ks[s, d]
        end
    end
    for d = 1:dof
        yerr[d] = dt*(ytrial[d] - yerr[d])
        ytrial[d] = y[d] + dt * ytrial[d]
    end
end

function stepsize_hw92!(dt, tdir, x0, xtrial, xerr, order,
                       timeout, dof, abstol, reltol, maxstep, norm)
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
    fac = [0.8, 0.9, 0.25^(1/(order + 1)), 0.38^(1/(order + 1))][1]
    facmax = 5.0 # maximal step size increase. 1.5-5
    facmin = 1.0/facmax  # maximal step size decrease. ?

    # in-place calculate xerr./tol
    for d = 1:dof
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

function calc_next_k!(ks::Matrix, ytmp::Vector, y, s, fn, t, dt, dof, btab)
    # Calculates the next ks and puts it into ks[s]
    # - ks and ytmp are modified inside this function.

    ytmp[:] = y
    for ss = 1:(s - 1), d = 1:dof
        ytmp[d] += dt * ks[ss, :][d]*btab.a[s, ss]
    end
    ks[s, :] = fn(t + btab.c[s]*dt, ytmp)
    nothing
end
