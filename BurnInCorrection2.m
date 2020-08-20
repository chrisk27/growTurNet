function [m] = BurnInCorrection2(t,xout)
%Note: This was directly copied from STAR Code (from Scholes et al.)
%
%This is a function to estimate a steady state from a damped oscillator.

%Input(s): 
% - t: vector define time stamp values of simulation
% - xout: nxt dimensional matrix that defines the state variables for each
% node and timepoint.

%Output(s):
% - m: vector of median value of simulation for each node from t>500 to
% t-final

    for i = 1:length(t)
        if t(i) > 500
            j = i;
            break
        end
    end
    m = median(xout(j:length(t),:));
end