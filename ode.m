function dxdt=ode(~, x, f, k)
%Note: This is copied from STAR method (Scholes et al.)
%
%Simple function that provides necessary inputs for running the ode solvers
%ode15s, ode23 and ode45.
%The function and parameter values need to be provided by the user (f, k)
dxdt = f(x, k);
end