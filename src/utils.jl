################################################################################
# Helper functions:
################################################################################

# in Hairer & Wanner describe how we can get an Hermite interpolation of order 5
# without additional function calls in pgs. 191-194. The coefficient used are
# based on extensions to the original dormand prince versions by shampine.
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

function setup_hermite!(sys)
    ydiff = sys.work.ytrial - sys.work.yinit
    bspl = sys.work.dt*sub(sys.work.ks, :, 1) - ydiff
    sys.work.ycont[:, 1] = sys.work.yinit
    sys.work.ycont[:, 2] = ydiff
    sys.work.ycont[:, 3] = bspl
    sys.work.ycont[:, 4] = -sys.work.dt*sub(sys.work.ks, :, 7) + ydiff - bspl
    sys.work.ycont[:, 5] = sys.work.dt*sys.work.ks*hermite_shampine_coef
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
islaststep(sys::AbstractODESystem) = sys.work.laststep

function hairer_norm(sys, abstol, reltol)
    # In harier and wanner pg. 168 the following error norm is described, which they mention
    # others have used max norms instead like we currently do. But the fortran codes use
    # this formulation
    # Mathematically the formula is:
    # err = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(\frac{y_{1i} - \yhat{y_{1i}}}{sc_i}^2)}
    # with sc_i = atol_i + \max(\abs(y_{0i}), \abs(y_{1i}))*rtol_i
    total = 0.0
    for i in 1:sys.work.ydim
        sc = abstol + reltol*max(abs(sys.work.yinit[i]), abs(sys.work.ytrial[i]))
        total += sqrt(((sys.work.ks[:, 1] - sys.work.ks[:, 7])/sys.work.ydim)^2)
    end
    return total
end
