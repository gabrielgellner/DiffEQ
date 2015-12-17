function dydt = rigid(t, y)
dydt = zeros(3, 1);
dydt(1) = y(2)*y(3);
dydt(2) = -y(1)*y(3);
dydt(3) = -0.51*y(1)*y(2);
end

