%% Load in the proper file containing the SS_Array and the Cond_Array
%If it's alreadly loaded, will choose output path
baseFolder = '\\files.brandeis.edu\epstein-lab\CKonow\Turing Networks\Network Results\';
if exist('SS_Array', 'var') == 1 && exist('Cond_Array', 'var') == 1
    if exist('savePath', 'var') ~= 1
        savePath = uigetdir(baseFolder, 'Select directory to save ');
    else
        disp('Using existing variables and path in workspace')
    end
else
    while exist('Cond_Array', 'var') ~= 1 || exist('SS_Array', 'var') ~= 1
        [filename, savePath] = uigetfile(baseFolder, 'Select proper .mat file');
        load(strcat(savePath, filename));
    end
end

%% Determine which conditions yield Turing patterns for each rate
numGrowth = length(Cond_Array) - 1;
numEntries = length(Cond_Array{1, 2});

storageMat = zeros(numEntries, numGrowth);
for i = 1:numGrowth %Compiles whether there are Turing patterns or not
    cond1 = Cond_Array{i, 2}(:,3:end);
    cond2 = Cond_Array{i, 3}(:,3:end);
    cond3 = Cond_Array{i, 4}(:,3:end);
    cond4 = Cond_Array{i, 5}(:,3:end);
    
    storageMat(:, i) = checkTuringCond(cond1, cond2, cond3, cond4);
end

%% Total Area Plot (may functionalize later)
growAxis = zeros(1, numGrowth);
totalSuccesses = zeros(1, numGrowth);

for i = 1: numGrowth
    growAxis(i) = Cond_Array{i, 1};
    totalSuccesses(i) = sum(storageMat(:, i));
end

plot(growAxis, totalSuccesses, 'Marker', 'o', 'XLabel', 'Growth Factor',...
    'YLabel', 'Number of Successful Conditions')

