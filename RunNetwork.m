%Run SystemDefinition_example first. Also, please use your own email
%settings (listed below) to report completion.

%Set up the file system that you will use. I put in my default settings,
%but you should change this if you want to save your own data
savePath = uigetdir('\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\Network Results\',...
    strcat('Select directory to save ', ID));
savePath = strcat(savePath, '\');

emailAddress = 'ckonow@brandeis.edu';
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
SteadyState_save = zeros(2*total,n+4); %This is a matrix that saves all steady states indexed to each parameter set.
Turing_save = zeros(2*total,4); %A matrix to save which parameter combination and corresponding steady state form Turing Instabilities. In addition, the diffusion values, and what type of Turing Instability is present, are defined.
tspan = [0 t_final]; %Vector defining time for solving ODE equation using the Matlab ODE suite
count_i = 1;
len = length(c_ini(:,1));
Saver=zeros(len^n,n);
counter = 2+n+k_length; 
final = counter;

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
        rxn_func{2}(x, k) - h_inf * x(2)];%Steady-state function, that should solve to zero
    
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
disp('Steady State search finished in: ');

toc %Finish timing this part of the process

%Now, I am writing another script to determine whether the steady state is
%stable or not in the absence of diffusion. This is much more complicated
%(and deviates much more from Scholes et al.'s STAR Code, since I am going
%to try and use the first two stability conditions from Madsvamuse et al.
%instead of the "standard" linear stability analysis. However, parts of
%this will indeed look like linear stability analysis, since I need to
%perturb around the steady state calculated in SteadyStateFinder.

disp('Analyzing Steady States.');
tic;

%The first portion is building up the functions that are required to 
X = sym('x', [1 n]); %Creates symbolic nodes
K = sym('k', [1 k_length]); %Creates symbolic parameter values
rxn_sym = [rxn_func{1}(X, K); rxn_func{2}(X, K)]; %Creates a symbolic function of the reaction component
rxn_jac = jacobian(rxn_sym, X); %Calculates Jacobian Matrix
rxn_jac_func = matlabFunction(rxn_jac); %Converts jacobian to symbolic function
rxn_jac_in = symvar(rxn_jac); % read out input parameters of jacobian Z
Xoverlap = ismember(X,rxn_jac_in); %find nodes that have to be substituted in Jacobian
Koverlap = ismember(K, rxn_jac_in); %find parameters that have to be substituted in Jacobian

count = 0; %Counts how many systems meet conditions for Turing patterns

for i = 1:length(SteadyState_save(:, 1))
    if System_save(SteadyState_save(i,1),length(k)+1) ~= 1 && SteadyState_save(i,n+2) ~= 2 %Exclude oscillators and unstable solutions
        
        %Calculate appropriate timescales. The growth timescale, T_gr, is the
        %infinmum time that it takes for the domain size to increase by a
        %factor of e. The dynamics timescale, T_dyn, is the maximum of the 
        %diffusion timescale and the kinetics timescale. Since we are only
        %considering aspymtotically slow growth, we will use the parameter
        %epsilon = T_dyn / T_gr << 1. This should be (by Madzvamuse et al.)
        %between 5 x 10^(-4) and 0.12
        %
        %To start, we will calculate the dynamic timescale at the initial
        %condition. This is done by finding the largest real part of the
        %eigenvalue of the jacobian at the intial condition.
        kk = System_save(SteadyState_save(i,1),1:k_length);
        x_init = SteadyState_save(i,2:n+1);
        init_input = num2cell([kk(Koverlap), x_init(Xoverlap)]);
        eig_init = eig(rxn_jac_func(init_input{:})); %Calculate the eigenvalues
        max_eig_init = max(abs(real(eig_init))); %Find max real part of the eigenvalue
        T_dyn = 1 / (gamma*max_eig_init); %Dynamic Timescale
        syms ts %Define symbolic time to solve for growth timescale
        T_gr = eval(solve(exp(1) == grow_func(ts, kk), ts)); %Growth timescale
        if isempty(T_gr) %If there isn't a growth timescale (i.e. no growth)
            T_gr = inf;
        end
        
        %Next, we must calculate the midpoint time, t_mid. This is where we
        %will evaluate the concentration profiles of the two nodes, and
        %then check the stability.
        t_star = T_dyn /2; %Uses the midpoint approximation
        %Lay out the system of ODEs we need to solve
        nonDiffEq = @(x, t, k, m, gamma) [gamma * rxn_func{1}(x, k) ...
            - m * dilution_func(t, k) * x(1); gamma * rxn_func{2}(x, k) ...
            - m * dilution_func(t, k) * x(2)];
        grow_tspan = [0 t_star];
        try
            [t_nonDiff, conc_nonDiff] = ode15s(@(t, x)nonDiffEq(x, t, kk, m, gamma),...
                grow_tspan, x_init);
        catch
            [t_nonDiff, conc_nonDiff] = ode23s(@(t, x)nonDiffEq(x, t, kk, m, gamma),...
                grow_tspan, x_init);
        end
        x_tstar = conc_nonDiff(end, :);
        
        %Now, we will check the first two stability conditions
        h_tstar = dilution_func(t_star, kk);
        stability_jac_input = num2cell([kk(Koverlap), x_tstar(Xoverlap)]);
        stability_jac = rxn_jac_func(stability_jac_input{:});
        tr_jac = trace(stability_jac);
        det_jac = det(stability_jac);
        SteadyState_save(i, n+3) = tr_jac;
        SteadyState_save(i, n+4) = det_jac;
        
        cond_1 = (-1) * gamma * tr_jac + 2*h_tstar;
        cond_2 = (-1) * h_tstar * gamma * tr_jac +...
            gamma^2 * det_jac;
        
        %Fill in the conditioned values
        SteadyState_save(i, n+5) = cond_1;
        SteadyState_save(i, n+6) = cond_2;
        if cond_1 <= 0 || cond_2 <= 0 % If homogeneous steady state is not stable
            SteadyState_save(i,n+2) = 1;
            System_save(SteadyState_save(i,1),k_length+2) = System_save(SteadyState_save(i,1),k_length+2)-1;
            %SteadyState_save(i, n+5) = NaN;
            %SteadyState_save(i, n+6) = NaN;
        else %If it seems stable, continue on to test for Turing Instability
            %Note: I'm ignoring most of what the STAR code does for
            %creating different diffusion ratios, since I'm only looking at
            %two nodes, so it doesn't need to generalize to 3 or more. So,
            %a lot of the stuff that was generated in the definition step
            %is useless.
            %
            %Will always assume the diffusion constant for the first node
            %is one, and then loop through the other possibilities listed
            %in d_range 
            for j = 1:length(d_range)
                d_val = 1 / d_range(j);
                cond_3 = (-1) * gamma * (d_val*stability_jac(1, 1) + stability_jac(2, 2)) ...
                    + h_tstar * (1 + d_val);
                cond_4 = (h_tstar*(1+d_val) - gamma * (d_val*stability_jac(1, 1) + stability_jac(2, 2)))^2 ...
                    - 4*d_val*(gamma^2 * det_jac - gamma*h_tstar*tr_jac);
                %SteadyState_save(i, n+5) = cond_3;
                %SteadyState_save(i, n+6) = cond_4;
                
                %Now, check to see if conditions are met
                if cond_3 < 0 && cond_4 > 0 %If it does meet conditions
                    count = count + 1;
                    Turing_save(count, 1) = i; %Saves index for steady state conditions
                    Turing_save(count, 2) = 1;
                    Turing_save(count, 3) = d_range(j); %Saves diffusion of node 1 (above) and node 2
                    Turing_save(count, 4:7) = [cond_1, cond_2, cond_3, cond_4];
                end
            end
        end
    end
end

disp('Finished Turing analysis of steady states in: ')
toc

% Delete all empty entries in Turing_save
%Delete all zero entries
for i = 1:length(Turing_save(:,1))
    if Turing_save(i,1) ==0 
        j = i-1;
        break
    end
end

Turing_save = Turing_save(1:j,:); %Get rid of unused entries except for 1 

%Save outputs to a file named ID (see Systemdefinition_example.m)
save(strcat(savePath,ID),'System_save','SteadyState_save','Turing_save', 'savePath', 'rxn_sym')

try
    setpref('Internet', 'E_mail', 'matlabck27@gmail.com');  % sender "from" address, typically same as username, e.g. 'xyz@gmail.com'
    setpref('Internet', 'SMTP_Username', 'matlabck27');
    setpref('Internet', 'SMTP_Password', 'turing1234');
    setpref('Internet', 'SMTP_Server',   'smtp.gmail.com');
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth',                'true');  % Note: 'true' as a string, not a logical value!
    props.setProperty('mail.smtp.starttls.enable',     'true');  % Note: 'true' as a string, not a logical value!
    props.setProperty('mail.smtp.socketFactory.port',  '465');   % Note: '465'  as a string, not a numeric value!
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    sendmail(emailAddress, 'Turing Network Simulations',...
    'Please go look at results.')
catch
    disp('Error in email settings.')
end