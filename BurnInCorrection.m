function [m] = BurnInCorrection(t,xout)
%Note: This was directly copied from STAR Code (from Scholes et al.)
%This is a function to calculate a threshold for estimating local
%maxima using function findpeaks from a given ODE trajectory.

%Input(s): 
% - t: vector define time stamp values of simulation
% - xout: nxt dimensional matrix that defines the state variables for each
% node and timepoint.

%Output(s):
% - m: 1% of maximum value of node trajectory from time t > 500 to t_final

    for i = 1:length(t)
        if t(i) > 500
            j = i;
            break
        end
    end
    m = max(xout(j:length(t),:))*0.01;
end