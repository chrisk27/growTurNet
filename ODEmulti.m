function [saver] = ODEmulti(c_ini,n,k,f_ode,tspan)
%Note: This was directly copied from STAR Code (from Scholes et al.)
%
%This is a function to run multiple ODE simualtions for different initial
%conditions defined in c_ini. Each trajectorie's endpoint is defined as
%possible steady state and this matrix is outputted from the function.

%Input(s)
% - c_ini: 
% - n: node number
% - k: parameters
% - f_ode: ode function
% - tspan: time span for ODE simulation

%Output(s)
% - saver: set of steady state suggestions. Each row represents one steady
% state estimation, columns represent nodes x(1)...x(n)

    %Intialise saver matrix
    saver = zeros(length(c_ini(:,1)),n);
    %Simulate ODE's for all possible initial conditions (c_ini)
    for p1 = 1:length(c_ini(:,1))
       y_0 = c_ini(p1,:);
       try [~,xout] = ode15s(@(t, x)ode(t,x,f_ode,k),tspan,y_0);
       catch
          [~, xout] = ode23s(@(t, x)ode(t,x,f_ode,k),tspan,y_0);
       end
       saver(p1,:)=xout(length(xout(:,1)),:); %save endpoints of trajectory.
    end
end