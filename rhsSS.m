function dxdt = rhsSS(t, x, k, dilution)
    dx1dt = k(4) - k(5)*x(1) + k(1)*(x(1)^2*k(2)^(-2) + x(2)^2*k(3)^(-2)) * ...
        (1 + x(1)^2*k(2)^(-2) + x(2)^2*k(3)^(-2))^(-1);
    dx2dt = k(8) - k(9)*x(2) + k(6) * (1 + x(1)^2*k(7)^(-2))^(-1);
    dxdt = [dx1dt - dilution(t, k); dx2dt - dilution(t, k)];
end
    