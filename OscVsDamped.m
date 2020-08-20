function state = OscVsDamped(Peak1,Peak2) 
%Note: This was directly copied from STAR Code (from Scholes et al.)
%
%This function takes two peaks from an ODE trajctory and checks whether a 
%damped oscillation or oscillation is expected given a threshold definition of: 
%(Peak1-Peak2)/Peak1 < 0.001. 

%Input(s)
% - Peak1: third to last local maxima found in trajectory
% - Peak2: second to last local maxima found in trajectory

%Output(s)
% - state: 2 defines oscillator, 3 defines damped oscillator

    if (Peak1-Peak2)/Peak1 < 0.001
        state = 2;
    else
        state = 3;
    end
end