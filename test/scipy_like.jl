# here we take as many of the tests from the scipy.odeint as we can
# The license for scipy is:
#=
Copyright (c) 2001, 2002 Enthought, Inc.
2 All rights reserved.
3

4 Copyright (c) 2003-2012 SciPy Developers.
5 All rights reserved.
6

7 Redistribution and use in source and binary forms, with or without
8 modification, are permitted provided that the following conditions are met:
9

10   a. Redistributions of source code must retain the above copyright notice,
11      this list of conditions and the following disclaimer.
12   b. Redistributions in binary form must reproduce the above copyright
13      notice, this list of conditions and the following disclaimer in the
14      documentation and/or other materials provided with the distribution.
15   c. Neither the name of Enthought nor the names of the SciPy Developers
16      may be used to endorse or promote products derived from this software
17      without specific prior written permission.
18

19

20 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
21 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
22 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
23 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
24 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
25 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
26 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
27 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
28 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
29 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
30 THE POSSIBILITY OF SUCH DAMAGE.
=#

# with the specific byline for the scipy.odeint tests:
# Authors: Nils Wagner, Ed Schofield, Pauli Virtanen, John Travers

############################
# End Legal
############################
using DiffEQ
using Base.Test

## Problems
"""Free vibration of a simple oscillator::
    m \ddot{u} + k u = 0, u(0) = u_0 \dot{u}(0) \dot{u}_0
Solution::
    u(t) = u_0*cos(sqrt(k/m)*t)+\dot{u}_0*sin(sqrt(k/m)*t)/sqrt(k/m)
"""
function simple_oscillator(t, y)
    k = 4.0
    m = 1.0
    A= zeros(2, 2)
    A[1, 2] = 1.0
    A[2, 1] = k/m
    return A*y
end

function simple_oscillator_analytic(t)
    k = 4.0
    m = 1.0 ##TODO I should put this in the top scope
    z0 = [1.0, 0.1]
    omega = sqrt(k/m)
    z0[1]*cos(omega*t) + z0[2]*sin(omega*t)/omega
end

tout = linspace(0, 100, 100)
sol = ode45(simple_oscillator, [0.5, 0.5], tout)
@test isapprox(sol.y[:, 1],  simple_oscillator_analytic(tout))
