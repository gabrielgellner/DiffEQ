format long

y0 = [0.0 1.0 1.0];
tspan = linspace(0.0, 12.0, 1000);

% this is from a matlab run with tol = 1e-13, the maximum base 10 accuracy 
refsol = [-0.705397809522307 -0.708811632467256 0.863846690370256];

toli = 2:10;
ntol = length(toli);
tols = zeros(ntol, 1);
for i = 1:ntol
    tols(i) = 10.0^-toli(i);
end

data = zeros(ntol, 3);
nsteps = zeros(ntol, 1);
for i = 1:ntol
    options = odeset('RelTol', tols(i), 'AbsTol', tols(i));
    sol = ode45(@rigid, tspan, y0, options);
    nsteps(i) = sol.stats.nsteps;
    data(i, :) = abs(refsol - deval(sol, tspan(end))');
end

if 1
    loglog(tols, data)
    % add a line showing the minus one exponent of the tolerance to see how
    % the solver is doing for precision
    hold on
    loglog(tols(2:end), tols(1:end - 1), 'k--') 
    hold off
end

dlmwrite('../matlab_rigid.csv', [tols data nsteps], 'precision', 16)