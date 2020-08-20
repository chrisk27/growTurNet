%% Description
%This script will take the range of values that meet all Turing conditions,
%and re-simulate them with actual time values and constants, etc. The
%objective is to make sure that the true ODE-simulated parameter space
%matches the one that the quick analysis (New_RunNetwork) shows. I do
%expect there to be some of the "messy" phase space cut out, but the major
%part of it should hopefully stay relatively constant.
%
%This script should run independently of New_DefineSystem and
%New_RunNetwork, as long as there is a SS_Array and a Cond_Array available,
%which we can load in. In addition, you need to type in a new rxn_function,
%which is the first step below. For now, I'll just assume that we're doing
%exponential growth so we don't have to worry about linear, but I'll code
%that in later as well.

%% Reaction Function
rxn_func = {@(x,k)(k(4) - k(5)*x(1) + k(1)*((x(1)/k(2))^2 + (x(2)/k(3))^2)...
    * (1 + (x(1)/k(2))^2 + (x(2)/k(3))^2)^(-1));
     @(x,k)(k(8) - k(9)*x(2) + k(6) / (1 + (x(1)/k(7))^2))};

%% Check and/or Load in Arrays
baseFile = '\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\Network Results\';
if ~exist('Cond_Array', 'var') || ~exist('SS_Array', 'var')
    currentFolder = cd(baseFile);
    [fileName, pathName] = uigetfile('*.mat');
    load(strcat(pathName, fileName))
    cd(currentFolder)
    savePath = strcat(pathName, '\');
else
    savePath = uigetdir(baseFile, 'Select directory to save');
end

%% Pull Basic Information from Arrays
tic;
%h values
numGR = size(Cond_Array, 1) - 1; %Number of growth rates tested
h_range = zeros(1, numGR);
for i = 1:numGR
    h_range(i) = Cond_Array{i, 1};
end

%k grid (not sure if I'll need the whole grid yet)
k_input = SS_Array{end, 2};
k_grid = combvec(k_input{:});

%d values
d_range = Cond_Array{end, 3};

%Use symbolic arrays to see if the reaction function is correct
X = sym('x', [1 2]); %Creates symbolic nodes
K = sym('k', [1 length(k_input)]); %Creates symbolic parameter values
rxn_sym = [rxn_func{1}(X, K); rxn_func{2}(X, K)]; %Creates a symbolic function of the reaction component

if rxn_sym(1) ~= Cond_Array{end, 1}(1) || rxn_sym(2) ~= Cond_Array{end, 1}(2)
    ME = MException('RxnFunc:NotAMatch', 'Reaction functions do not match');
    throw(ME)
end

rxn_jac = jacobian(rxn_sym, X); %Calculates Jacobian Matrix
rxn_jac_func = matlabFunction(rxn_jac); %Converts jacobian to symbolic function
rxn_jac_in = symvar(rxn_jac); % read out input parameters of jacobian Z
Xoverlap = ismember(X,rxn_jac_in); %find nodes that have to be substituted in Jacobian
Koverlap = ismember(K, rxn_jac_in); %find parameters that have to be substituted in Jacobian

if rxn_jac ~= Cond_Array{end, 2}
    ME1 = MException('RxnJac:NotAMatch', 'Reaction functions do not match');
    throw(ME1)
end

%% Set up other parameters and export containers
%Other parameters
gamma = 1; %Scales influence of reaction on PDE
n = 2; %Number of nodes in Network
x_max = 20001;
int = 10000;
m =1; %#of spatial dimensions
c_ini = permn(1:int:x_max,n); %Grid from which algorithm will sample. Here 3 initial conditions

%Set up variable parameters
arrayLoc = [];
for j = 1:size(SS_Array{end, 1}, 2) %Find the right array locations
    str = SS_Array{end, 1}{j}{1};
    arrayLoc(j) = str2double(str(2:end));
end
        
%Export containers - Will have an equivalent Cond_Array (newCondVals) and a
%new state array that will hold the x values (concentrations) like the
%SS_Array
newStateArray = SS_Array(1:size(SS_Array, 1)-2, :);
for i = 1:size(newStateArray, 1)
    newMat = zeros(size(SS_Array{i, 2}));
    newMat(:, 1:2) = SS_Array{i, 2}(:, 1:2);
    newStateArray{i, 2} = newMat;
end

newCondArray = Cond_Array(1:size(Cond_Array, 1) - 1, :);
for i = 1:size(newCondArray, 1)  % Delete existing condition value matrices
    for j = 2:size(newCondArray, 2)
        newMat = zeros(size(Cond_Array{i, j}));
        newMat(:, 1:2) = Cond_Array{i, j}(:, 1:2);
        newCondArray{i, j} = newMat;
    end
end


%% Run Simulations and Record Results
for gr = 1:length(h_range)
    grow_val = h_range(gr);
    SS_List = SS_Array{gr, 2};
    
    stateMat = newStateArray{gr, 2};
    newCond1 = newCondArray{gr, 2};
    newCond2 = newCondArray{gr, 3};
    newCond3 = newCondArray{gr, 4};
    newCond4 = newCondArray{gr, 5};
    
    for i = 1:length(k_grid(1, :))
        kk = k_grid(:, i)';
        
        %First, check to make sure we're looking at the right set of
        %parameters
        for j = 1:length(arrayLoc)
            if kk(arrayLoc(j)) ~= SS_List(i, j)
                ME2 = MException('RxnParams:NotAMatch', 'Variable Parameters do not match.');
                throw(ME2)
            end
        end
        
        %Calculate the appropriate timescales using the eigenvalues
        x_init = SS_List(i, 3:end);
        if any(isnan(x_init)) ||any(x_init == 0)
            continue
        end
        init_input = num2cell([kk(Koverlap), x_init(Xoverlap)]);
        eig_init = eig(rxn_jac_func(init_input{:})); %Calculate the eigenvalues
        max_eig_init = max(abs(real(eig_init))); %Find max real part of the eigenvalue
        T_dyn = 1 / (gamma*max_eig_init); %Dynamic Timescale
        t_star = T_dyn /2; %Uses the midpoint approximation
        
        %Lay out the system of ODEs we need to solve numerically
        nonDiffEq = @(x, t, k, gr_val, m, gamma) [gamma * rxn_func{1}(x, k) ...
            - m * gr_val * x(1); gamma * rxn_func{2}(x, k) ...
            - m * gr_val * x(2)];
        grow_tspan = [1000, 1000+t_star];
        try
            [t_nonDiff, conc_nonDiff] = ode15s(@(t, x)nonDiffEq(x, t, kk, grow_val, m, gamma),...
                grow_tspan, x_init);
        catch
            [t_nonDiff, conc_nonDiff] = ode23s(@(t, x)nonDiffEq(x, t, kk, grow_val, m, gamma),...
                grow_tspan, x_init);
        end
        x_tstar = conc_nonDiff(end, :);
        stateMat(i, 3:end) = x_tstar; %inserts value back into the output array
        
        %Calculate new values using the jacobian
        stability_jac_input = num2cell([kk(Koverlap), x_tstar(Xoverlap)]);
        stability_jac = rxn_jac_func(stability_jac_input{:});
        tr_jac = trace(stability_jac);
        det_jac = det(stability_jac);
        
        %Calculate condition values
        %Plug in Condition Values for first two conditions
        newCond1(i, end) = (-1) * gamma * tr_jac + 2*grow_val;
        newCond2(i, end) = (-1) * grow_val * gamma * tr_jac + gamma^2 * det_jac;
        
        for d = 1:length(d_range)
            d_val = 1 / d_range(d); %Diffusion constant ratio
            newCond3(i, d + length(x_init)) = (-1) * gamma * ...
                (d_val*stability_jac(1, 1) + stability_jac(2, 2)) ...
                    + grow_val * (1 + d_val);
            newCond4(i, d + length(x_init)) = (grow_val*(1+d_val) - gamma...
                * (d_val*stability_jac(1, 1) + stability_jac(2, 2)))^2 ...
                - 4*d_val*(gamma^2 * det_jac - gamma*grow_val*tr_jac);
        end
    end
    
    %Replace all condition and state matrices back into their arrays
    newStateArray{gr, 2} = stateMat;
    newCondArray{gr, 2} = newCond1;
    newCondArray{gr, 3} = newCond2;
    newCondArray{gr, 4} = newCond3;
    newCondArray{gr, 5} = newCond4;
end

%% Plot and compare to Previous Results
%Left side will be old plots, right side will be new plots

%Determine axis labels
xlb = min(k_input{arrayLoc(1)});
xub = max(k_input{arrayLoc(1)});
ylb = min(k_input{arrayLoc(2)});
yub = max(k_input{arrayLoc(2)});

XtickVals = logspace(log10(xlb), log10(xub), log10(xub) - log10(xlb) + 1);
YtickVals = logspace(log10(ylb), log10(yub), log10(yub) - log10(ylb) + 1);

for gr = 1:length(h_range)
    oldCond1 = Cond_Array{gr, 2};
    oldCond2 = Cond_Array{gr, 3};
    oldCond3 = Cond_Array{gr, 4};
    oldCond4 = Cond_Array{gr, 5};
    newCond1 = newCondArray{gr, 2};
    newCond2 = newCondArray{gr, 3};
    newCond3 = newCondArray{gr, 4};
    newCond4 = newCondArray{gr, 5};
    
    fig = figure('Units', 'inches', 'Position', [2.5 1 9 7], 'Name', ...
        strcat('h = ', string(h_range(gr))));
    
    %Condition 1 Plots
    ax1o = subplot(4, 4, 1);
    CondPlot2D(oldCond1(:, 1), oldCond1(:, 2), oldCond1(:, 3:end));
    ax1o.FontSize = 8;
    ax1o.XScale = 'log';
    ax1o.YScale = 'log';
    ax1o.XLimMode = 'manual';
    ax1o.YLimMode = 'manual';
    ax1o.XLim = [xlb, xub];
    ax1o.YLim = [ylb, yub];
    ax1o.XTick = XtickVals;
    ax1o.YTick = YtickVals;
    ax1o.Title.String = 'Old Condition 1';
    
    ax1n = subplot(4, 4, 3);
    CondPlot2D(newCond1(:, 1), newCond1(:, 2), newCond1(:, 3:end));
    ax1n.FontSize = 8;
    ax1n.XScale = 'log';
    ax1n.YScale = 'log';
    ax1n.XLimMode = 'manual';
    ax1n.YLimMode = 'manual';
    ax1n.XLim = [xlb, xub];
    ax1n.YLim = [ylb, yub];
    ax1n.XTick = XtickVals;
    ax1n.YTick = YtickVals;
    ax1n.Title.String = 'New Condition 1';
    
    %Condition 2
    ax2o = subplot(4, 4, 2);
    CondPlot2D(oldCond2(:, 1), oldCond2(:, 2), oldCond2(:, 3:end));
    ax2o.FontSize = 8;
    ax2o.XScale = 'log';
    ax2o.YScale = 'log';
    ax2o.XLimMode = 'manual';
    ax2o.YLimMode = 'manual';
    ax2o.XLim = [xlb, xub];
    ax2o.YLim = [ylb, yub];
    ax2o.XTick = XtickVals;
    ax2o.YTick = YtickVals;
    ax2o.Title.String = 'Old Condition 2';
    
    ax2n = subplot(4, 4, 4);
    CondPlot2D(newCond2(:, 1), newCond2(:, 2), newCond2(:, 3:end));
    ax2n.FontSize = 8;
    ax2n.XScale = 'log';
    ax2n.YScale = 'log';
    ax2n.XLimMode = 'manual';
    ax2n.YLimMode = 'manual';
    ax2n.XLim = [xlb, xub];
    ax2n.YLim = [ylb, yub];
    ax2n.XTick = XtickVals;
    ax2n.YTick = YtickVals;
    ax2n.Title.String = 'New Condition 2';
    
    %Condition 3
    ax3o = subplot(4, 4, 5);
    CondPlot2D(oldCond3(:, 1), oldCond3(:, 2), oldCond3(:, 3:end), 'ConditionGreaterThan', false);
    ax3o.FontSize = 8;
    ax3o.XScale = 'log';
    ax3o.YScale = 'log';
    ax3o.XLimMode = 'manual';
    ax3o.YLimMode = 'manual';
    ax3o.XLim = [xlb, xub];
    ax3o.YLim = [ylb, yub];
    ax3o.XTick = XtickVals;
    ax3o.YTick = YtickVals;
    ax3o.Title.String = 'Old Condition 3';
    
    ax3n = subplot(4, 4, 7);
    CondPlot2D(newCond3(:, 1), newCond3(:, 2), newCond3(:, 3:end), 'ConditionGreaterThan', false);
    ax3n.FontSize = 8;
    ax3n.XScale = 'log';
    ax3n.YScale = 'log';
    ax3n.XLimMode = 'manual';
    ax3n.YLimMode = 'manual';
    ax3n.XLim = [xlb, xub];
    ax3n.YLim = [ylb, yub];
    ax3n.XTick = XtickVals;
    ax3n.YTick = YtickVals;
    ax3n.Title.String = 'New Condition 3';
    
    %Condition 4
    ax4o = subplot(4, 4, 6);
    CondPlot2D(oldCond4(:, 1), oldCond4(:, 2), oldCond4(:, 3:end), 'ConditionGreaterThan', true);
    ax4o.FontSize = 8;
    ax4o.XScale = 'log';
    ax4o.YScale = 'log';
    ax4o.XLimMode = 'manual';
    ax4o.YLimMode = 'manual';
    ax4o.XLim = [xlb, xub];
    ax4o.YLim = [ylb, yub];
    ax4o.XTick = XtickVals;
    ax4o.YTick = YtickVals;
    ax4o.Title.String = 'Old Condition 4';
    
    ax4n = subplot(4, 4, 8);
    CondPlot2D(newCond4(:, 1), newCond4(:, 2), newCond4(:, 3:end), 'ConditionGreaterThan', true);
    ax4n.FontSize = 8;
    ax4n.XScale = 'log';
    ax4n.YScale = 'log';
    ax4n.XLimMode = 'manual';
    ax4n.YLimMode = 'manual';
    ax4n.XLim = [xlb, xub];
    ax4n.YLim = [ylb, yub];
    ax4n.XTick = XtickVals;
    ax4n.YTick = YtickVals;
    ax4n.Title.String = 'New Condition 4';
    
    %Full Turing Space
    ax5o = subplot(4, 4, [9, 10, 13, 14]);
    successes = checkTuringCond(oldCond1(:, 3:end), oldCond2(:, 3:end), oldCond3(:, 3:end), oldCond4(:, 3:end));
    coloring = zeros(length(oldCond1(:,1)), 3);
    for i = 1:length(successes)
        if successes(i)
            coloring(i, :) = [1 0 0];
        end
    end
    xVal = oldCond1(:, 1);
    yVal = oldCond2(:, 2);
    scatter(xVal, yVal, 12, coloring, 'MarkerFaceColor', 'flat');
    ax5o.FontSize = 12;
    ax5o.XScale = 'log';
    ax5o.YScale = 'log';
    ax5o.XLimMode = 'manual';
    ax5o.YLimMode = 'manual';
    ax5o.XLim = [xlb, xub];
    ax5o.YLim = [ylb, yub];
    ax5o.XTick = XtickVals;
    ax5o.YTick = YtickVals;
    ax5o.Title.String = 'Old Turing Space';
    
    ax5n = subplot(4, 4, [11, 12, 15, 16]);
    successes = checkTuringCond(newCond1(:, 3:end), newCond2(:, 3:end), newCond3(:, 3:end), newCond4(:, 3:end));
    coloring = zeros(length(newCond1(:,1)), 3);
    for i = 1:length(successes)
        if successes(i)
            coloring(i, :) = [1 0 0];
        end
    end
    xVal = newCond1(:, 1);
    yVal = newCond2(:, 2);
    scatter(xVal, yVal, 12, coloring, 'MarkerFaceColor', 'flat');
    ax5n.FontSize = 12;
    ax5n.XScale = 'log';
    ax5n.YScale = 'log';
    ax5n.XLimMode = 'manual';
    ax5n.YLimMode = 'manual';
    ax5n.XLim = [xlb, xub];
    ax5n.YLim = [ylb, yub];
    ax5n.XTick = XtickVals;
    ax5n.YTick = YtickVals;
    ax5n.Title.String = 'New Turing Space';
    
    figureName = 'compareImg' + string(gr) + '_h_' + string(h_range(gr)) + '_TuringSpace.png';
    saveas(gcf, char(savePath+figureName));    
end
toc
    