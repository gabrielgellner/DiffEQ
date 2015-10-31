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

# isoutofdomain takes the state and returns true if state is outside
# of the allowed domain.  Used in adaptive step-control.
isoutofdomain(x) = isnan(x)
