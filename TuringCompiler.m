%This script will collect all of simulations that yield Turing patterns
%into one file. I will base the extracellular and intracellular robustness
%measures off of this (maybe - we'll see how it plays out). To work, you
%must have the Turing_save, System_save, and SteadyState_save files from a
%simulation loaded into the workplace (can be done by loading the .mat
%file). If you don't, it will open up a gui to go find it.
%
%Organization of CompilerOutput matrix:
%Each row is different set of conditions that form Turing patterns on
%growing domains.
%Cols 1 & 2 are the indices for the conditions in the System_save and
%SteadyState_save files, respectively.
%Cols in range param_cols are the reaction parameters. Note: will need to
%keep track of what each variable means some other way.
%Col growth_col is the growth parameter "r".
%The two columns after growth_col are the steady-state values.
%Cols in diff_cols contain the diffusion coefficients.
%The last 4 colums contain the values of the l.h.s. of the Turing
%conditions in Madvamuze.

while exist('Turing_save', 'var') == 0
    [selectedFile, sFPath] = uigetfile('\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\*.mat',...
        'Choose file to compile');
    load(strcat(sFPath,selectedFile))
end

ID_compiled = 'Compiled_Turing_conditions';

output = zeros(length(Turing_save(:, 1)), 8 + length(System_save(1, :)));

for i = 1:length(Turing_save(:, 1)) %Fill in output
    
    %First, input the indices to find where the rest of the data is
    %located. The first entry (col 1) is the System_save index, and the
    %second entry (col 2) is the SteadyState_save index.
    SS_save_idx = Turing_save(i, 1);
    Sys_save_idx = SteadyState_save(SS_save_idx, 1);
    output(i, 1) = Sys_save_idx;
    output(i, 2) = SS_save_idx;
    
    %Now, insert parameter and steady state values
    rxn_params = length(System_save(1, :)) - 3;
    growth_r_param = length(System_save(1, :)) - 2;
    output(i, 3:rxn_params+2) = System_save(Sys_save_idx, 1:rxn_params); %Adds in reaction parameters
    output(i, rxn_params+3) = System_save(Sys_save_idx, growth_r_param); %Adds in growth parameter
    
    SS_vals_idx_start = rxn_params + 4; %So I don't have to keep typing this myself
    output(i, SS_vals_idx_start:SS_vals_idx_start+1) = SteadyState_save(SS_save_idx, 2:3); %Puts in steady state conc vals
    
    %Now, Diffusion values
    Diff_idx = SS_vals_idx_start+2;
    output(i, Diff_idx:Diff_idx + 1) = Turing_save(i, 2:3);
    
    %Finally, the calculated condition values (I just moved these to the
    %Turing_save file)
    output(i, Diff_idx+2:Diff_idx+5) = Turing_save(i, 4:7);
    
end

param_cols = 3:rxn_params+2;
growth_param = rxn_params+3;
diff_cols = Diff_idx:Diff_idx+1;
CompilerOutput = output;

if exist('savePath', 'var')
    save(strcat(savePath, ID_compiled), 'CompilerOutput', 'param_cols', 'growth_param', 'diff_cols');
elseif exist('sFPath', 'var')
    save(strcat(sFPath, ID_compiled), 'CompilerOutput', 'param_cols', 'growth_param', 'diff_cols');
else
    disp('Must save manually')
end
        
        
    
    
    
    
    

