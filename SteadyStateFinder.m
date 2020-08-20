%Run SystemDefinition_example first

%__________________________________________________________________________
%ODE optimisation parameters
options = optimoptions('fsolve','FunctionTolerance', 1.0e-12,'Display','off'); %ODE solver setting see, Matlab documentation for further details
t_final = 1000; %Time threshold for ODE simulation
%__________________________________________________________________________
%Initialise some matrices (Note: Some of these will only be used when doing
%Turint analysis, not just steady state analysis. Will sort later.
y0 = rand.*ones(n, 1); %Initial conditions for first simulation     
total = length(k_grid(1,:)); %Calculates number of iterations needed to sample all parameter combinations
z = k_length+2; 
System_save = zeros(total,z); %System save is a matrix to save all parameters as well as information to the amount of Steady States found for each parameter set and the classfication of the system according to 0,1,2,3 (Stable,unstable,oscillating,damped oscillating). The behavior is estimated from the simulation trajectories.
SteadyState_save = zeros(2*total,n+6); %This is a matrix that saves all steady states indexed to each parameter set.
Turing_save = zeros(2*total,n+3); %A matrix to save which parameter combination and corresponding steady state form Turing Instabilities. In addition, the diffusion values, and what type of Turing Instability is present, are defined.
tspan = [0 t_final]; %Vector defining time for solving ODE equation using the Matlab ODE suite
count_i = 1;
len = length(c_ini(:,1));
Saver=zeros(len^n,n);
counter = 2+n+k_length; 
final = counter;
%warning('off', 'optim:fsolve:NonSquareSystem'); %Turns off a repeating warning that occurs when fsolve switches method to solve. May correct later.
%__________________________________________________________________________
%Workflow Part 1: Find steady states numerically

tic;

for i = 1:length(k_grid(1,:)) %Loops over all sets of parameters
    k = k_grid(:,i)'; %Choose the set of parameters I want
    
    %To find the steady state, first need to find h(inf) as the dilution
    %term. Then, will write the function for the steady-state
    syms t_sym
    h_inf = double(limit(dilution_func(t_sym, k), t_sym, inf)); %Dilution function at infinite time
    ss_func = @(x, k) [rxn_func{1}(x, k) - h_inf * x(1);...
        rxn_func{2}(x, k) - h_inf *x(2)];%Steady-state function, that should solve to zero
    
    %Run first ODE simulation
    try
        [t1, xout1] = ode15s(@(t, x)ode(t, x, ss_func, k), tspan, y0);
    catch
        [t1, xout1] = ode23s(@(t, x)ode(t, x, ss_func, k), tspan, y0);
    end
    
    ss2 = xout1(length(xout1(:,1)),:); %Suggest a steady state as last entry in the time series
    max_thresh = BurnInCorrection(t1,xout1); %Calculate a peak threshold as 1% of max value of each node
    state = 0;
    
    %Check to see if system is oscillating (normal or damped)
    min_max_count = 0;
    for p1 = 1:n
        local_max = (findpeaks(xout1(:,p1),'MinPeakProminence',max_thresh(p1)));
        if length(local_max) > 1
            min_max_count = min_max_count + 1;
            local_max_end(p1) = local_max(length(local_max)-1);
            local_max_save = local_max;
        end
    end
    %Distinguish between damped oscillations from normal oscillations
    if min_max_count > 0
        state = OscVsDamped(local_max_save(length(local_max_save)-1), local_max_save(length(local_max_save)));
        ss = median(xout1(:,:));
    end
    
    %Now: Evaluate if the steady state is even a possibilty to look into
    %solutions (elimate non-steady states or oscillations, keep constant
    %values and damped oscillations)

    %Exclude negative and imaginary solutions
    ss2(ss2<=0) = 0;
    if all(ss2) == 0 || isreal(ss2) < 1
        state = 1;
    end
    
    if state == 0 % For systems that don't oscillate and that seem stable
        %Will check for multiple steady states by simulating the trajectory
        %with multiple initial conditions (defined by c_ini)
        Saver = ODEmulti(c_ini, n, k, ss_func, tspan); %Performs ODE simulations and outputs endpoints
        %Perform cluster analysis to check for redundant entries
        C = ClusterAnalysis(Saver, n, k, ss_func, tspan);
        if isempty(C) == 1
            state = 1;
        end
        %Save all non-redundant steady states
        for p1 = 1:length(C(:,1))
            System_save(i, length(k)+2) = length(C(:,1)); %Save number of steady states
            SteadyState_save(count_i,1) = i;
            SteadyState_save(count_i, 2:n+1) = C(p1,:); %save steady sate (steady state values for each node)
            count_i = count_i + 1;
        end
    end
    
    if state == 2 %System is found to be an oscillator and will not be further analyzed
        SteadyState_save(count_i,1) = i;
        SteadyState_save(count_i, 2:n+1) = ss;
        SteadyState_save(count_i, n+2) = state;
        count_i = count_i + 1;
    end
    
    if state == 3 %If it's a damped oscillator
        %Similar to non-oscillating system, will scan for multiple initial
        %conditions
        Saver = ODEmultiOsc(c_ini, n, k, ss_func, tspan);
        %Perform Cluster Analysis: using f_solve, solutions are optimized
        %and then k_means clustering is used to distinguish steady states
        C = ClusterAnalysisOsc(Saver, n, k, ss_func, options);
        System_save(i, length(k)+2) = length(C(:,1)); %save number of steady states
        for p1 = 1:length(C(:,1))
            SteadyState_save(count_i,1) = i;
            SteadyState_save(count_i,2:n+1) = C(p1,:); %Save steady state values for each node
            count_i = count_i + 1;
        end
    end
    
    if state == 1 %if system doesn't converge sufficiently within the time given in tspan, define system as unstable
        System_save(i,length(k)+2) = 0;
    end
    
    System_save(i, 1:length(k)) = k;
    System_save(i, length(k) + 1) = state;
end
    
%Eliminate all empty entries
for i = 1:length(SteadyState_save(:,1))
    if SteadyState_save(i, 1) == 0
        j = i- 1;
        break
    end
end

SteadyState_save = SteadyState_save(1:j,:);

toc %Finish timing this part of the process
    
    
    
    
        
    
        
    
    
    
