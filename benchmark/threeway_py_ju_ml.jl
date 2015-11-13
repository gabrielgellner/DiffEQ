# This code tests the same simple problem using three implementations using the same
# dopri5 methods
# 1) scipy's ode using the dopri5 method as this should be a straight wrapper of the fotran code
# 2) ODE.jl's version
# 3) DiffEQ.jl's modified ODE.jl's version
# All of these codes ill be compared to a matlab's ode45 version at an extremely high
# tolerance (atol = rtol = 1e-13)
#
# The codes will look at how going from low to high tolerances how they compare

using ODE
using DiffEQ
using PyCall
using PyPlot

# scipy imports
@pyimport scipy.integrate as ivp

function rigid(t, y)
    ydot = zeros(3)
    ydot[1] = y[2]*y[3]
    ydot[2] = -y[1]*y[3]
    ydot[3] = -0.51*y[1]*y[2]
    ydot
end

tout = linspace(0.0, 12.0, 100)
y0 = [0.0, 1.0, 1.0]

# the rigid system from matlab using atol = rtol = 1e-13
reft = 12.0
refsol = [-0.705397809522307, -0.708811632467256, 0.863846690370256]

#tols = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
tols = [10.0^-i for i in 2:10]
data_de = Array{Float64}(3, length(tols))
data_o = Array{Float64}(3, length(tols))
data_sp = Array{Float64}(3, length(tols))
for (i, tol) in enumerate(tols)
    sol = aode(Dopri54(rigid, y0), tout, abstol = tol, reltol = tol)
    data_de[:, i] = abs(refsol - sol.y[end, :]')

    t, y = ode45(rigid, y0, tout, abstol = tol, reltol = tol, points = :specified)
    data_o[:, i] = abs(refsol - y[end])

    # scipy's interface is a bit verbose
    r = ivp.ode(rigid)[:set_integrator]("dopri5", atol = tol, rtol = tol, safety = 0.8)
    r[:set_initial_value](y0, 0.0)
    ylast = r[:integrate](tout[end])
    data_sp[:, i] = abs(refsol - ylast)
end

# read in the same run in matlab 2015a
data_ml = readdlm("matlab_rigid.csv", ',')
println(size(data_ml))

subplot(131)
loglog(tols, data_de'[:, 1], label = "aode")
loglog(tols, data_o'[:, 1], "--k", label = "ode45")
loglog(tols, data_sp'[:, 1], "-.r", label = "sp")
loglog(data_ml[:, 1], data_ml[:, 2], "m", label = "matlab")
# it is usefull to see a tol:tol[-1] line to see how the tolerance vs accuracy is working
loglog([10.0^-i for i in 2:10], [10.0^(1 - i) for i in 2:10] , ":g")
xlabel("tol")
ylabel("error in y1 abs(ref - est)")
legend(loc = 2)

subplot(132)
loglog(tols, data_de'[:, 2], label = "aode")
loglog(tols, data_o'[:, 2], "--k", label = "ode45")
loglog(tols, data_sp'[:, 2], "-.r", label = "sp")
loglog(data_ml[:, 1], data_ml[:, 3], "m", label = "matlab")
# it is usefull to see a tol:tol[-1] line to see how the tolerance vs accuracy is working
loglog([10.0^-i for i in 2:10], [10.0^(1 - i) for i in 2:10] , ":g")
xlabel("tol")
ylabel("error in y2 abs(ref - est)")
legend(loc = 2)

subplot(133)
loglog(tols, data_de'[:, 3], label = "aode")
loglog(tols, data_o'[:, 3], "--k", label = "ode45")
loglog(tols, data_sp'[:, 3], "-.r", label = "sp")
loglog(data_ml[:, 1], data_ml[:, 4], "m", label = "matlab")
# it is usefull to see a tol:tol[-1] line to see how the tolerance vs accuracy is working
loglog([10.0^-i for i in 2:10], [10.0^(1 - i) for i in 2:10] , ":g")
xlabel("tol")
ylabel("error in y3 abs(ref - est)")
legend(loc = 2)

show()
