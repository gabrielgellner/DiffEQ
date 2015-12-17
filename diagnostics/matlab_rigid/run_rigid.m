options = odeset('RelTol',1e-13,'AbsTol',1e-13);
tout = linspace(0.0, 12.0, 10000);
tic;
[T,Y] = ode45(@rigid, tout, [0.0 1.0 1.0], options);
toc
format long
disp(Y(end, :))