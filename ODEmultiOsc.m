function [saver] = ODEmultiOsc(c_ini,n,k,f_ode,tspan) 
%Note: This was directly copied from STAR Code (from Scholes et al.)
%
%This is a function to run multiple ODE equations for systems that show 
%damped oscillations, the time span for simulations is expanded to 2*tspan 
%to allow for sufficient convergence. The steady state is estimated by 
%calculating the median of the ODE trajectories after a burn in correction 
%of 0.5*tspan.
    
%Input(s)
% - c_ini: 
% - n: node number
% - k: parameters
% - f_ode: ode function
% - tspan: time span for ODE simulation

%Output(s)
% - saver: set of suggested steady states (redundant entries present). Each
% row defines one steady state. The columns represent the nodes
% x(1)...x(n)

    %Initialise saver
    saver = zeros(length(c_ini(:,1)),n);
    %Simulate for all possible initial conditions (c_ini)
    for p1 = 1:length(c_ini(:,1))
       y_0 = c_ini(p1,:);
       try [t,xout] = ode15s(@(t, x)ode(t,x,f_ode,k),tspan*2,y_0);
       catch
          [t, xout] = ode23s(@(t, x)ode(t,x,f_ode,k),tspan*2,y_0);
       end
       %Estimate steady state by calculating the median of trajectory after
       % burn in time 0.5*t
       saver(p1,:)= BurnInCorrection2(t,xout);
    end
end