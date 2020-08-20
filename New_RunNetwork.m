%% Set Up Save Path - Make sure to run New_DefineSystem.m before
savePath = uigetdir('\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\Network Results\',...
    strcat('Select directory to save ', ID));
savePath = strcat(savePath, '\');

emailAddress = 'ckonow@brandeis.edu';

%% Optimization Parameters
options = optimoptions('fsolve','MaxFunctionEvaluations',2000,'FunctionTolerance', 1.0e-8,'Display','off','FiniteDifferenceType','central'); %ODE solver setting see, Matlab documentation for further details
t_final = 1000; %Time threshold for initial ODE simulation

%% Initialize Matrices and other Values
SS_Array = cell(length(h_range) + 2, 2); %This is the cell array where all the steady state data is stored
    %It is structured as follows:
    %Each row is a different h value, which is stored in the first column
    %The second column is a matrix of steady state values, where the first
    %0-3 columns are the values of the variables being changed, and the
    %last n columns are the steady-state values for the nodes.
    %The second-to-last row contains a symbolic reaction function (col 1)
    %and the symbolic jacobian (col 2).
    %The last row contains a cell array for the variables that are changing
    %(col 1) and the k_input array (col 2)
for i = 1:length(h_range)
    SS_Array{i, 1} = h_range(i);
end

%Determine scanning variables
numMats = 0;
scanVars = {};
idxChange = [];
for i = 1:length(k_input) %Determines number of variables that we're spanning
    if isscalar(k_input{i}) == 0
        numMats = numMats + 1;
        scanVars{numMats} = strcat('k' + string(i)); %saves the variables that we're changing (in order)
        idxChange(numMats) = i;
    end
end
SS_Array{end, 1} = scanVars;
SS_Array{end, 2} = k_input;

%Jacobian Stuff
X = sym('x', [1 n]); %Creates symbolic nodes
K = sym('k', [1 k_length]); %Creates symbolic parameter values
rxn_sym = [rxn_func{1}(X, K); rxn_func{2}(X, K)]; %Creates a symbolic function of the reaction component
rxn_jac = jacobian(rxn_sym, X); %Calculates Jacobian Matrix
rxn_jac_func = matlabFunction(rxn_jac); %Converts jacobian to symbolic function
rxn_jac_in = symvar(rxn_jac); % read out input parameters of jacobian Z
Xoverlap = ismember(X,rxn_jac_in); %find nodes that have to be substituted in Jacobian
Koverlap = ismember(K, rxn_jac_in); %find parameters that have to be substituted in Jacobian
SS_Array{length(h_range) + 1, 1} = rxn_sym;
SS_Array{length(h_range) + 1, 2} = rxn_jac;

%Define the basic matrix which each growth value will test from
baseMat = zeros(length(k_grid(1, :)), length(idxChange) + n);
for i = 1:length(k_grid(1, :))
    for j = 1:length(idxChange)
        baseMat(i, j) = k_grid(idxChange(j), i);
    end
end
startSS = length(idxChange) + 1; %Index to mark start of where I put SS

%% Numerically Integrate for h = 0 case
tic;
tspan_init = [0 1000];
h_0_mat = baseMat;
for i = 1:length(k_grid(1, :))
    %Use ODE solvers to try original steady state
    rng shuffle %chooses a new seed
    rand_init = rand.*ones(n, 1);
    kk = k_grid(:, i)';
    ss_func = @(x,k) [rxn_func{1}(x,k); rxn_func{2}(x,k)];
    try
        [t1, xout1] = ode15s(@(t, x)ode(t, x, ss_func, kk), tspan_init, rand_init);
    catch
        [t1, xout1] = ode23s(@(t, x)ode(t, x, ss_func, kk), tspan_init, rand_init);
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
        ss2 = median(xout1(:,:));
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
        Saver = ODEmulti(c_ini, n, kk, ss_func, tspan_init); %Performs ODE simulations and outputs endpoints
        %Perform cluster analysis to check for redundant entries
        C = ClusterAnalysis(Saver, n, kk, ss_func, tspan_init);
        if isempty(C) == 1
            state = 1;
        elseif length(C(:,1)) ~= 1
            disp(strcat('Multiple Steady States at index ' + string(i)))
        end
    end
    
    if state == 3 %If it's a damped oscillator
        %Similar to non-oscillating system, will scan for multiple initial
        %conditions
        Saver = ODEmultiOsc(c_ini, n, kk, ss_func, tspan_init);
        %Perform Cluster Analysis: using f_solve, solutions are optimized
        %and then k_means clustering is used to distinguish steady states
        C = ClusterAnalysisOsc(Saver, n, kk, ss_func, options);
        if isempty(C) == 1
            state = 1;
        elseif length(C(:, 1)) ~= 1
            disp(strcat('Multiple Osc Steady States at index ' + string(i)))
        end
    end
    
    if state == 1 || state == 2
        %disp(strcat('Could not find steady state at index ' + string(i)))
    elseif length(C(:,1)) ~= 1
        diffC = zeros(length(C(:, 1)), 1);
        for j = 1:length(C(:, 1))
            diffC(j) = abs(var(C(j, :)));
        end
        for j = 1:length(diffC)
            if diffC(j) == min(diffC)
                steadstat = C(j, :);
            end
        end
        disp(strcat('Had to choose minimal steady state at index ' +string(i)))
        h_0_mat(i, startSS:end) = steadstat;
    else
        h_0_mat(i, startSS:end) = C(:,:);
    end
end
SS_Array{1, 2} = h_0_mat;
toc

%% Use fsolve for all h ~= 0 cases. If error, use ODE solver
tic;
for gr = 2:length(h_range)
    h_gr_new = baseMat;
    grow_val = h_range(gr);
    ss_gr_func = @(x,k) [rxn_func{1}(x,k) - grow_val*x(1); rxn_func{2}(x,k) - grow_val*x(2)];
    lb = [eps; eps]; %Lower bounds on any solve problem
    ub = [Inf; Inf]; %Upper bounds
    for i = 1:length(k_grid(1, :))
        kk = k_grid(:, i)';
        prevSS = SS_Array{gr - 1, 2}(i, startSS:end); %Steady state with previous h value
        if all(prevSS == 0)
            x_guess = rand(2, 1);
        else
            x_guess = prevSS';
        end
        
        try
            %x_new_SS = lsqnonlin(@(x)ss_gr_func(x, kk), x_guess, lb, ub, options);
            x_new_SS = fsolve(@(x)ss_gr_func(x, kk), x_guess, options);
            if any(x_new_SS <= 0)
                ErX = MException('ConcValue:NonRealistic', 'Unreasonable concentration value obtained');
                %disp('Fsolve gave concentration less than or equal to zero')
                throw(ErX)
            end
        catch
            %disp(strcat('Could not use fsolve on index ', string(i)))
            rando_init = rand(2, 1);
            try
                [t2, xout2] = ode15s(@(t, x)ode(t, x, ss_gr_func, kk), tspan_init, rando_init);
            catch
                [t2, xout2] = ode23s(@(t, x)ode(t, x, ss_gr_func, kk), tspan_init, rando_init);
            end
            
            ss3 = xout2(length(xout2(:,1)),:); %Suggest a steady state as last entry in the time series
            state = 0;

            %Now: Evaluate if the steady state is even a possibilty to look into
            %solutions (elimate non-steady states or oscillations, keep constant
            %values and damped oscillations)

            %Exclude negative and imaginary solutions
            ss3(ss3<=0) = 0;
            if all(ss3) == 0 || isreal(ss3) < 1
                state = 1;
            end

            if state == 0 % For systems that don't oscillate and that seem stable
                %Will check for multiple steady states by simulating the trajectory
                %with multiple initial conditions (defined by c_ini)
                Saver = ODEmulti(c_ini, n, kk, ss_gr_func, tspan_init); %Performs ODE simulations and outputs endpoints
                %Perform cluster analysis to check for redundant entries
                C = ClusterAnalysis(Saver, n, kk, ss_gr_func, tspan_init);
                if isempty(C) == 1
                    state = 1;
                elseif length(C(:,1)) ~= 1
                    disp(strcat('Multiple Steady States at index ' + string(i)))
                end
            end
            
            if state == 1 || state == 2
                disp(strcat('Could not find steady state at index ', string(i), ' for h = ', string(h_range(gr))))
            elseif length(C(:,1)) ~= 1
                diffC = zeros(length(C(:, 1)));
                for j = length(C(:, 1))
                    diffC(j) = abs(var(C(j, :)));
                end
                for j = length(diffC)
                    if diffC(j) == min(diffC)
                        steadstat = C(j, :);
                    end
                end
                disp(strcat('Had to choose minimal steady state at index ', string(i), ' for h = ', string(h_range(gr))))
                x_new_SS = steadstat';
            else
                x_new_SS = C(:,:)';
            end
        end
        
        h_gr_new(i, startSS:end) = x_new_SS';
    end
    if all(x_new_SS) >= 0
        SS_Array{gr, 2} = h_gr_new;
    end
end
toc

%% Generate Storage for Calculated Condition Values
%First, will set up array for condition values. The first column will store
%the h values, and the last 4 columns will store the 4 conditions (in order
%as presented by Madzvamuse et al.). The last row will contain the same as
%SS_Array, but will also have the diffusion range in the third column.
Cond_Array = cell(length(h_range)+1, 5);
for i = 1:length(h_range)+1
    Cond_Array{i, 1} = SS_Array{i, 1};
end
Cond_Array{end, 2} = SS_Array{length(h_range) + 1, 2};
Cond_Array{end, 3} = d_range;

%Now, I will set up the different condition matricies we need to store in
%the Cond_Array. The first two conditions (which do not deal with
%diffusion) only have one entry in addition to the variable parameters.
%However, the last two, which do deal with diffusion, will have
%length(d_range) additional columns.
condMat = zeros(length(k_grid(1,:)), length(idxChange) + 1);
condMat_long = zeros(length(k_grid(1,:)), length(idxChange) + length(d_range));
condMat(:, 1:length(idxChange)) = h_0_mat(:, 1:length(idxChange));
condMat_long(:, 1:length(idxChange)) = h_0_mat(:, 1:length(idxChange));

%Now, loop through each of the growth values and calculate the conditions.
%Store all of the results in a condMat, which will then be placed into the
%Cond_Array.
for gr = 1:length(h_range)
    grow_val = h_range(gr);
    ss_vals = SS_Array{gr, 2}(:, startSS:end); %Matrix containing just the steady state values.
    cM1 = condMat;
    cM2 = condMat;
    cM3 = condMat_long;
    cM4 = condMat_long;
    
    for i = 1:length(k_grid(1, :))
        kk = k_grid(:, i)';
        xVal = ss_vals(i, :);
        stability_jac_input = num2cell([kk(Koverlap), xVal(Xoverlap)]);
        stability_jac = rxn_jac_func(stability_jac_input{:});
        tr_jac = trace(stability_jac);
        det_jac = det(stability_jac);
        
        %Plug in Condition Values for first two conditions
        cM1(i, end) = (-1) * gamma * tr_jac + 2*grow_val;
        cM2(i, end) = (-1) * grow_val * gamma * tr_jac + gamma^2 * det_jac;
        
        %Fill in Condition Values for last two conditions for each
        %diffusion value
        for d = 1:length(d_range)
            d_val = 1 / d_range(d); %Diffusion constant ratio
            cM3(i, d + length(xVal)) = (-1) * gamma * (d_val*stability_jac(1, 1) + stability_jac(2, 2)) ...
                    + grow_val * (1 + d_val);
            cM4(i, d + length(xVal)) = (grow_val*(1+d_val) - gamma * (d_val*stability_jac(1, 1) + stability_jac(2, 2)))^2 ...
                - 4*d_val*(gamma^2 * det_jac - gamma*grow_val*tr_jac);
        end
    end
    
    %Put condition matrices back into the array
    Cond_Array{gr, 2} = cM1;
    Cond_Array{gr, 3} = cM2;
    Cond_Array{gr, 4} = cM3;
    Cond_Array{gr, 5} = cM4;
end

%% Plot or Find Some Way to Visualize
%Determine axis labels
xlb = min(k_input{idxChange(1)});
xub = max(k_input{idxChange(1)});
ylb = min(k_input{idxChange(2)});
yub = max(k_input{idxChange(2)});

XtickVals = logspace(log10(xlb), log10(xub), log10(xub) - log10(xlb) + 1);
YtickVals = logspace(log10(ylb), log10(yub), log10(yub) - log10(ylb) + 1);

for gr = 1:length(h_range)
    %For each condition, using the sineage (> or < 0) in Klika et al. 2017
    cond1 = Cond_Array{gr, 2};
    cond2 = Cond_Array{gr, 3};
    cond3 = Cond_Array{gr, 4};
    cond4 = Cond_Array{gr, 5};
 
    
    fig = figure('Units', 'inches', 'Position', [4.5 1 4.5 7], 'Name', ...
        strcat('h = ', string(h_range(gr))));
    
    %Condition 1 Plot
    ax1 = subplot(4, 2, 1);
    CondPlot2D(cond1(:, 1), cond1(:, 2), cond1(:, 3:end));
    ax1.FontSize = 8;
    ax1.XScale = 'log';
    ax1.YScale = 'log';
    ax1.XLimMode = 'manual';
    ax1.YLimMode = 'manual';
    ax1.XLim = [xlb, xub];
    ax1.YLim = [ylb, yub];
    ax1.XTick = XtickVals;
    ax1.YTick = YtickVals;
    ax1.Title.String = 'Condition 1';
    
    
    %Condition 2 Plot
    ax2 = subplot(4, 2, 2);
    CondPlot2D(cond2(:, 1), cond2(:, 2), cond2(:, 3:end));
    ax2.FontSize = 8;
    ax2.XScale = 'log';
    ax2.YScale = 'log';
    ax2.XLimMode = 'manual';
    ax2.YLimMode = 'manual';
    ax2.XLim = [xlb, xub];
    ax2.YLim = [ylb, yub];
    ax2.XTick = XtickVals;
    ax2.YTick = YtickVals;
    ax2.Title.String = 'Condition 2';
    
    %Condition 3 Plot
    ax3 = subplot(4, 2, 3);
    CondPlot2D(cond3(:, 1), cond3(:, 2), cond3(:, 3:end), 'ConditionGreaterThan', false);
    ax3.FontSize = 8;
    ax3.XScale = 'log';
    ax3.YScale = 'log';
    ax3.XLimMode = 'manual';
    ax3.YLimMode = 'manual';
    ax3.XLim = [xlb, xub];
    ax3.YLim = [ylb, yub];
    ax3.XTick = XtickVals;
    ax3.YTick = YtickVals;
    ax3.Title.String = 'Condition 3';
    
    %Condition 4 Plot
    ax4 = subplot(4, 2, 4);
    CondPlot2D(cond4(:, 1), cond4(:, 2), cond4(:, 3:end), 'ConditionGreaterThan', true);
    ax4.FontSize = 8;
    ax4.XScale = 'log';
    ax4.YScale = 'log';
    ax4.XLimMode = 'manual';
    ax4.YLimMode = 'manual';
    ax4.XLim = [xlb, xub];
    ax4.YLim = [ylb, yub];
    ax4.XTick = XtickVals;
    ax4.YTick = YtickVals;
    ax4.Title.String = 'Condition 4';
    
    %Full Turing Space Plot
    ax5 = subplot(4, 2, 5:8);
    
    %Plot Turing Space
    successes = checkTuringCond(cond1(:, 3:end), cond2(:, 3:end), cond3(:, 3:end), cond4(:, 3:end));
    coloring = zeros(length(cond1(:,1)), 3);
    for i = 1:length(successes)
        if successes(i)
            coloring(i, :) = [1 0 0];
        end
    end
    xVal = cond1(:, 1);
    yVal = cond2(:, 2);
    scatter(xVal, yVal, 12, coloring, 'MarkerFaceColor', 'flat');
    ax5.FontSize = 12;
    ax5.XScale = 'log';
    ax5.YScale = 'log';
    ax5.XLimMode = 'manual';
    ax5.YLimMode = 'manual';
    ax5.XLim = [xlb, xub];
    ax5.YLim = [ylb, yub];
    ax5.XTick = XtickVals;
    ax5.YTick = YtickVals;
    ax5.Title.String = 'Meets Turing Conditions';
    
    %Save figure to Current folder
    figureName = 'img' + string(gr) + '_h_' + string(h_range(gr)) + '_TuringSpace.png';
    saveas(gcf, char(savePath+figureName));    
    close;
end
    

%% Save everything to some file (as defined by folder before)
filename = strcat(savePath, ID, '_newAnalysis.mat');
save(filename, 'SS_Array', 'Cond_Array')

disp('Done with Simulations')

%Email myself to let me know I'm finished
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

