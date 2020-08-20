%This script determines the "extracellular robustness" of the network of
%its given Turing conditions. In essence, it determines the proportion of
%diffusion parameters studied actually meet the criteria for Turing
%patterns.
%
%This code will end up looking very similar to IntracellularRobustness.m.
%Similarly, you will need to have loaded both the three matrix outputs from
%RunNetwork, as well as the outputs from TuringCompiler.m
%
%The output of this script will be (saved as Diff_Robustness):
%total_diffusion_robustness, which is the total number of diffusion
%conditions that yield the Turing conditions divided by the total number of
%diffusion conditions tests AND
%growthRobustness_diffusion, an array that lists the diffusion robustness
%per each "growth" parameter

%Load in the compiled Turing conditions
while exist('CompilerOutput', 'var') == 0
    [selectedCompiledFile, sCFPath] = uigetfile('\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\*.mat',...
        'Choose compiled file to analyze');
    load(strcat(sCFPath,selectedCompiledFile))
end

%Note: Since we don't use diffusion until the very end (until it's done
%checking for the first two stability conditions to be met), unlike the
%IntracellularRoubustness there isn't a good list of what conditions are
%there. Instead, I'm going to write this with the assumption that we're
%always using the below code to check for (all) diffusion coefficients
%analyzed, as was used both in my original simulations and Scholes et al.
d_range = logspace(-3,3,7); %This is the diffusion range that will be analysed. (From original code)
num_diffusion = length(d_range); %Number of different conditions that were tested

%First, I want to just see how many total different combos we see
different_diffusion = unique(CompilerOutput(:, diff_cols), 'rows');
num_different_diffusion = length(different_diffusion(:, 1));

total_diffusion_robustness = num_different_diffusion / num_diffusion;

%Then, I want to sort through and determine each per growth amount
SysSave_grow_col = length(System_save(1,:)) - 2;
all_growth_params = unique(System_save(:, SysSave_grow_col));
growthRobustness_Diffusion = zeros(length(all_growth_params), 2);
growthRobustness_Diffusion(:, 1) = all_growth_params';

for i = 1:length(all_growth_params)
    %Filters the domains to just the growth rate in question
    gp = all_growth_params(i);
    Compiled_filter = CompilerOutput(CompilerOutput(:, growth_param) == gp, :);
    
    grow_unique_diffusion = unique(Compiled_filter(:, diff_cols), 'rows');
    num_grow_diffusion = length(grow_unique_diffusion(:, 1));
    
    growthRobustness_Diffusion(i, 2) = num_grow_diffusion / num_diffusion;
end

if exist('savePath', 'var')
    save(strcat(savePath, 'Diff_Robustness'), 'total_diffusion_robustness', 'growthRobustness_Diffusion');
elseif exist('sFPath', 'var')
    save(strcat(sFPath, 'Diff_Robustness'), 'total_diffusion_robustness', 'growthRobustness_Diffusion');
else
    disp('Must save manually')
end

    
