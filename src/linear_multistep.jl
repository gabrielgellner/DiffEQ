##TODO: this code will not be used in the final versions. It is just to help get
## a sense of the structure of the code. I will need a real adaptive solver
function adams_fixedstep_solver(sys::AdamsSystem, tout;
                                h = 1e-2 # fixed stepsize
                                order = 2 # order of the Adams step between 1-12
                                )
    # setup system memory

end

"""linmult_initial_step

Setup the initial approximation history used by the multistep methods. Will use
a fixed step euler method
"""
function linmult_initial_step(sys, h0)
    ##TODO: this code is just copied and pasted from `hinit` in `runge_jutta.jl`
    tdir = sign(tend - tstart)
    tdir == 0 && error("Zero time span")
    tau = max(reltol*norm(sys.work.yinit, Inf), abstol)
    d0 = norm(sys.work.yinit, Inf)/tau
    f0 = sys.func(tstart, sys.work.yinit)
    d1 = norm(f0, Inf)/tau
    #if d0 < 1e-5 || d1 < 1e-5
    #    h0 = 1e-6
    #else
    #    h0 = 0.01*(d0/d1)
    #end
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
