################################################################################
# Helper functions:
################################################################################

# I have factored this out so that I can use it in a clean namespaced import
# so that it is clear that things like `vcat_nosplat` are not builtins.
# Currently I don't know how to do the import logic to get this. I have issues
# with getting the types from the `runge_kutta_tableaus.jl`

function hermite_interp!(y, tquery, t, dt, y0, y1, f0, f1)
    # For dense output see Hairer & Wanner p.190 using Hermite
    # interpolation. Updates y in-place.
    #
    # f_0 = f(x_0 , y_0) , f_1 = f(x_0 + h, y_1 )
    # this is O(3). TODO for higher order.
    theta = (tquery - t)/dt
    for i = 1:length(y0)
        y[i] = ((1 - theta)*y0[i] + theta*y1[i] + theta*(theta - 1) *
               ((1 - 2*theta)*(y1[i] - y0[i]) + (theta - 1)*dt*f0[i] + theta*dt*f1[i]))
    end
    nothing
end

function hermite_interp(tquery, t, dt, y0, y1, f0, f1)
    # Returns the y instead of in-place
    y = similar(y0)
    hermite_interp!(y, tquery, t, dt, y0, y1, f0, f1)
    return y
end

# estimator for initial step based on book
# "Solving Ordinary Differential Equations I" by Hairer et al., p.169
function hinit(fn, x0, t0, tend, p, reltol, abstol)
    # Returns first step, direction of integration and F evaluated at t0
    tdir = sign(tend - t0)
    tdir == 0 && error("Zero time span")
    tau = max(reltol*norm(x0, Inf), abstol)
    d0 = norm(x0, Inf)/tau
    f0 = fn(t0, x0)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    # perform Euler step
    x1 = x0 + tdir*h0*f0
    f1 = fn(t0 + tdir*h0, x1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 1e-15
        h1 = max(1e-6, 1e-3*h0)
    else
        pow = -(2.0 + log10(max(d1, d2)))/(p + 1.0)
        h1 = 10.0^pow
    end
    return tdir*min(100*h0, h1, tdir*(tend - t0)), tdir, f0
end

# isoutofdomain takes the state and returns true if state is outside
# of the allowed domain.  Used in adaptive step-control.
isoutofdomain(x) = isnan(x)
