################################################################################
# Helper functions:
################################################################################

# in Hairer & Wanner we have the following explainting of the dense output:
# given that for a single step t -> t + dt we have the following
# f(t, y(t)), f(t + dt, y(t + dt)) and y(t) and y(t + dt) then clearly the
# 3rd degree Hermite polynomial is available for interpolating between (t, t + dt).
# now to bootstrap this interpolant without additional funciton calls we can evaluate
#
const hermite_shampine_coef = [
    -12715105075/11282082432,
    0.0, #TODO: this will have to be 0 if I want to make these codes work for any AbstractFloat -- and I will need to explcitly convert
    87487479700/32700410799,
    -10690763975/1880347072,
    701980252875/199316789632,
    -1453857185/822651844,
    69997945/29380423
]

function setup_hermite!(sys, dt)
    ##TODO: I like the names ytrial -> y1, and yinit -> y0
    ydiff = sys.work.ytrial - sys.work.yinit
    bspl = dt*sub(sys.work.ks, :, 1) - ydiff
    sys.work.ycont[:, 1] = sys.work.yinit
    sys.work.ycont[:, 2] = ydiff
    sys.work.ycont[:, 3] = bspl
    sys.work.ycont[:, 4] = -dt*sub(sys.work.ks, :, 7) + ydiff - bspl
    sys.work.ycont[:, 5] = dt*sys.work.ks*hermite_shampine_coef
    return nothing
end

function hermite_shampine_interp!(y, tquery, t, dt, ycont)
    ##TODO: ycont could actually have a better name ... lke hcoef or something
    theta = (tquery - t)/dt
    theta1 = 1.0 - theta
    for i = 1:length(y)
        y[i] = ycont[i, 1] + theta*(ycont[i, 2] + theta1*(
            ycont[i, 3] + theta*(ycont[i, 4] + theta1*ycont[i, 5])
        ))
    end
    return nothing
end

function hermite_shampine_interp(tquery, t, dt, ycont)
    # Returns the y instead of in-place
    y = Array{Float64}(size(ycont, 1))
    hermite_shampine_interp!(y, tquery, t, dt, ycont)
    return y
end

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
