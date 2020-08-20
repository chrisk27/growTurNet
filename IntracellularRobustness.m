%This script is used for analyzing the "intracellular robustness" measure
%defined by Scholes et al. Basically, it will loop through the different
%values of the reaction parameters and then determine the fraction that
%yielded Turing patterns. I'll run a couple different sorts and store the
%results together, and then we'll be able to compare results more.
%
%To run, you will have needed both to load in the arrays from the Network
%Analysis script (those saved by RunNetwork) as well as compiled the
%results using TuringCompiler (or load in the array/column numbers from the
%output of that script.)
%
%Output of this script will be (all saved as Intra_Robustness.mat):
%The basicRobustness variable, which is just the number of intracellular
%parameter sets that yield Turing patterns divided by the total number of
%parameter sets tested AND
%The growthRobustness_Intra array, which stores the robustness for each of
%the "growth" parameters

%Load in (if not already done)
while exist('CompilerOutput', 'var') == 0
    [selectedCompiledFile, sCFPath] = uigetfile('\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\*.mat',...
        'Choose compiled file to analyze');
    load(strcat(sCFPath,selectedCompiledFile))
end
while exist('Turing_save', 'var') == 0
    [selectedFile, sFPath] = uigetfile('\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\*.mat',...
        'Choose original output file');
    load(strcat(sFPath,selectedFile))
end

%The first thing that we will do is just a basic robustness calculation,
%which will show how many parameter combinations actually yield Turing
%patterns.
last_SysSave_col = length(System_save(1,:)) - 3; %Subtract off three since last 2 cols detail stability/type of SS, and 3rd to last is growth param
unique_params_all = unique(System_save(:, 1:last_SysSave_col), 'rows'); %A list of all unique intracellular parameter sets
num_cond_tot = length(unique_params_all(:, 1));

unique_params_turing = unique(CompilerOutput(:, param_cols), 'rows'); %Unique parameters that yielded Turing condition
num_cond_turing = length(unique_params_turing(:, 1));

basicRobustness = num_cond_turing / num_cond_tot; %Number of conditions yielding TP / number of conditions tested

%Now, I want to sort it by growth rate and then determine robustness
SysSave_grow_col = length(System_save(1,:)) - 2;
all_growth_params = unique(System_save(:, SysSave_grow_col));
growthRobustness_Intra = zeros(length(all_growth_params), 2);
growthRobustness_Intra(:, 1) = all_growth_params';

for i = 1:length(all_growth_params)
    %Filters the domains to just the growth rate in question
    gp = all_growth_params(i);
    System_filter = System_save(System_save(:, SysSave_grow_col) == gp, :);
    Compiled_filter = CompilerOutput(CompilerOutput(:, growth_param) == gp, :);
    
    %Calculate robustness in much the same way as above
    cond_turing_grow = unique(Compiled_filter(:, param_cols), 'rows');
    num_turing_grow = length(cond_turing_grow(:, 1));
    cond_all_grow = unique(System_filter(:, 1:last_SysSave_col), 'rows');
    num_all_grow = length(cond_all_grow(:, 1));
    
    growthRobustness_Intra(i, 2) = num_turing_grow / num_all_grow;
end

if exist('savePath', 'var')
    save(strcat(savePath, 'Intra_Robustness'), 'basicRobustness', 'growthRobustness_Intra');
elseif exist('sFPath', 'var')
    save(strcat(sFPath, 'Intra_Robustness'), 'basicRobustness', 'growthRobustness_Intra');
else
    disp('Must save manually')
end


