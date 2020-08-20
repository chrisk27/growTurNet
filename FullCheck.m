%% Load in the "New" analysis file. Will use this to compare.
baseFolder = '\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\Network Results\';
if exist('SS_Array', 'var') == 1 && exist('Cond_Array', 'var') == 1
    if exist('savePath', 'var') ~= 1
        oldPath = uigetdir(baseFolder, 'Select Parent Directory');
    else
        disp('Using existing variables and path in workspace')
        oldPath = savePath;
    end
else
    while exist('Cond_Array', 'var') ~= 1 || exist('SS_Array', 'var') ~= 1
        [filename, oldPath] = uigetfile(baseFolder, 'Select proper .mat file');
        load(strcat(oldPath, filename));
    end
end

mkdir oldPath FullAnalysis

saveDir = strcat(oldPath, '\', 'FullAnalysis\');
%% Pull out all of the conditions and reaction equations
numGrowth = length(Cond_Array) - 1;
numEntries = length(Cond_Array{1, 2});

rxn1 = SS_Array{numGrowth+1, 1}(1); %symbolic reaction function of 1
rxn2 = SS_Array{numGrowth+1, 1}(2); %symbolic reaction function of 2

k_input = SS_Array{numGrowth+2, 2}; %reaction parameters
k_length = length(k_input);
k_grid = combvec(k_input{:});

h_range = zeros(numGrowth, 1); %Vector with every growth value
for i = 1:numGrowth
    h_range(i) = SS_Array{i, 1};
end

d_range = Cond_Array{numGrowth+1, 3}; %potential diffusion coefficients

%% Other Misc Parameters (usually won't change)
gamma = 1; %Scales influence of reaction on PDE
n = 2; %Number of nodes in Network
m = 1; %Number of dimensions (usually 1)
x_max = 20001;
int = 10000;
c_ini = permn(1:int:x_max,n); %Grid from which algorithm will sample. Here 3 initial condition
options = optimoptions('fsolve','MaxFunctionEvaluations',2000,'FunctionTolerance', 1.0e-8,'Display','off','FiniteDifferenceType','central'); %ODE solver setting see, Matlab documentation for further details
t_final = 1000; %Time threshold for initial ODE simulation

%% Recreate actual reaction functions
syms x1 x2
if k_length == 9
    syms k1 k2 k3 k4 k5 k6 k7 k8 k9
    rxn1func = matlabFunction(rxn1, 'Vars', {[x1; x2], [k1, k2, k3, k4, k5, k6, k7, k8, k9]});
    rxn2func = matlabFunction(rxn2, 'Vars', {[x1; x2], [k1, k2, k3, k4, k5, k6, k7, k8, k9]});
elseif k_length == 10
    syms k1 k2 k3 k4 k5 k6 k7 k8 k9 k10
    rxn1func = matlabFunction(rxn1, 'Vars', {[x1; x2], [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]});
    rxn2func = matlabFunction(rxn2, 'Vars', {[x1; x2], [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]});
end

%% Set Up Output Files and Other Containers
%First I need to figure out what I want to output, will do this later
New_State_Out = cell(length(h_range) + 2, 3);
for i = 1:length(h_range)
    New_State_Out{i, 1} = h_range(i);
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
New_State_Out{end, 1} = scanVars;
New_State_Out{end, 2} = k_input;

%Define the basic matrix which each growth value will test from
baseMat = zeros(length(k_grid(1, :)), length(idxChange) + n);
for i = 1:length(k_grid(1, :))
    for j = 1:length(idxChange)
        baseMat(i, j) = k_grid(idxChange(j), i);
    end
end
startSS = length(idxChange) + 1; %Index to mark start of where I put SS

%% WorkFlow Part 1: Find Steady States Numerically
tic;
tspan_init = [0 t_final];
for hSpot = 1:length(h_range)
    h_val = h_range(hSpot);
    hMat = baseMat;
    for i = 1:length(k_grid(1, :))
        %Use ODE solvers to find original steady state
        rng shuffle %chooses a new seed
        rand_init = rand.*ones(n, 1);
        kk = k_grid(:, i)';
        ss_func = @(x,k)[rxn1func(x,k) - h_val*x(1); rxn2func(x,k) - h_val*x(2)];
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
            hMat(i, startSS:end) = steadstat;
        else
            hMat(i, startSS:end) = C(:,:);
        end
    end
    New_State_Out{hSpot, 2} = hMat;
    hMessage = strcat('Done with h = ', string(h_val));
    disp(hMessage)
end
toc
disp('Finished Steady State Analysis')

%% Calculate Jacobian and all that jazz
X = [x1; x2];
if k_length == 9
    K = [k1, k2, k3, k4, k5, k6, k7, k8, k9];
elseif k_length == 10
    K = [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10];
end

rxnJac = jacobian([rxn1func(X, K), rxn2func(X, K)],X);
rxnJacFunc = matlabFunction(rxnJac);
rxnJacVar = symvar(rxnJac);
Xoverlap = ismember(X,rxnJacVar); %find nodes that have to be substituted in Jacobian
Koverlap = ismember(K, rxnJacVar); %find parameters that have to be substituted in Jacobian

%% Run through and do the last bit of integration for each of the growth values
tic;
for hSpot = 1:length(h_range)
    h_val = h_range(hSpot);
    hMatGrow = baseMat;
    for i = 1:length(k_grid(1, :))
        kk = k_grid(:, i)';
        x_steady = New_State_Out{hSpot, 2}(i, startSS:end)';
        for value = 1:length(x_steady)
            if x_steady(value) == 0
                x_steady(value) = eps^2;
            end
        end
        init_input = num2cell([kk(Koverlap), x_steady(Xoverlap)']);
        eig_init = eig(rxnJacFunc(init_input{:})); %Calculate the eigenvalues
        max_eig_init = max(abs(real(eig_init))); %Find max real part of the eigenvalue
        T_dyn = 1 / (gamma*max_eig_init); %Dynamic Timescale
        
        %Note: I'm only doing this for exponential type growth now, would
        %need to vary it to make linear (or any other type of growth) a
        %nonconstant function. In fact, this just extends it for a few
        %seconds longer than the initial steady state. So, probably won't
        %get too much interesting behavior. I think, up until I either add
        %in the diffusion to the integration (which I don't think I have to
        %based on Madsvamuse et al) or use linear growth, this won't do
        %much in the end.
        
        t_star = T_dyn /2; %Uses the midpoint approximation
        %Lay out the system of ODEs we need to solve
        nonDiffEq = @(x, t, k, m, gamma) [gamma * rxn1func(x, k) ...
            - m * h_val * x(1); gamma * rxn2func(x, k) ...
            - m * h_val * x(2)];
        grow_tspan = [0 t_star];
        
        try
            [t_nonDiff, conc_nonDiff] = ode15s(@(t, x)nonDiffEq(x, t, kk, m, gamma),...
                grow_tspan, x_steady);
        catch
            [t_nonDiff, conc_nonDiff] = ode23s(@(t, x)nonDiffEq(x, t, kk, m, gamma),...
                grow_tspan, x_steady);
        end
        x_tstar = conc_nonDiff(end, :);
        hMatGrow(i, startSS:end) = x_tstar;
    end
    New_State_Out{hSpot, 3} = hMatGrow;
    hMessage = strcat('Done with h = ', string(h_val));
    disp(hMessage)
end
toc
disp('End of time-sensitive "perturbation"')

%% Generate Storage for Calculated Condition Values
%First, will set up array for condition values. The first column will store
%the h values, and the last 4 columns will store the 4 conditions (in order
%as presented by Madzvamuse et al.). The last row will contain the same as
%SS_Array, but will also have the diffusion range in the third column.
New_Cond_Array = cell(length(h_range)+1, 5);
for i = 1:length(h_range)+1
    New_Cond_Array{i, 1} = New_State_Out{i, 1};
end
New_Cond_Array{end, 2} = New_State_Out{length(h_range) + 1, 2};
New_Cond_Array{end, 3} = d_range;

%Now, I will set up the different condition matricies we need to store in
%the Cond_Array. The first two conditions (which do not deal with
%diffusion) only have one entry in addition to the variable parameters.
%However, the last two, which do deal with diffusion, will have
%length(d_range) additional columns.
condMat = zeros(length(k_grid(1,:)), length(idxChange) + 1);
condMat_long = zeros(length(k_grid(1,:)), length(idxChange) + length(d_range));
condMat(:, 1:length(idxChange)) = New_State_Out{1, 2}(:, 1:length(idxChange));
condMat_long(:, 1:length(idxChange)) = New_State_Out{1, 2}(:, 1:length(idxChange));

%Now, loop through each of the growth values and calculate the conditions.
%Store all of the results in a condMat, which will then be placed into the
%Cond_Array.
for gr = 1:length(h_range)
    grow_val = h_range(gr);
    int_vals = New_State_Out{gr, 3}(:, startSS:end); %Matrix containing just the steady state values.
    cM1 = condMat;
    cM2 = condMat;
    cM3 = condMat_long;
    cM4 = condMat_long;
    
    for i = 1:length(k_grid(1, :))
        kk = k_grid(:, i)';
        xVal = int_vals(i, :)';
        stability_jac_input = num2cell([kk(Koverlap), xVal(Xoverlap)']);
        stability_jac = rxnJacFunc(stability_jac_input{:});
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
    New_Cond_Array{gr, 2} = cM1;
    New_Cond_Array{gr, 3} = cM2;
    New_Cond_Array{gr, 4} = cM3;
    New_Cond_Array{gr, 5} = cM4;
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
    cond1 = New_Cond_Array{gr, 2};
    cond2 = New_Cond_Array{gr, 3};
    cond3 = New_Cond_Array{gr, 4};
    cond4 = New_Cond_Array{gr, 5};
 
    
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
    saveas(gcf, char(saveDir+figureName));    
    close;
end
    

%% Save everything to some file (as defined by folder before)
filename = strcat(saveDir, filename, '_timeDiff.mat');
save(filename, 'New_State_Out', 'New_Cond_Array')

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
        



