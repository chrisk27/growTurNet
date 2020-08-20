function output = checkTuringCond(condition1, condition2, condition3, condition4, varargin)
%This function checks whether the Turing conditions are met, as per Klika
%et al.
%   To use this function, pass in the values or matrices that contain the 
%   left-hand-side condition values from Klika et al. 2017. Each vector or
%   matrix should have the same number of rows, though the columns can
%   vary. The optional entry will decide whether to accept or reject 
%   whether to accept matrix rows that have some values that meet the 
%   conditions and others that don't. The options are 'any' (as long 
%   as one value in the row meets condition), 'all' (all values must meet 
%   condition), and 'majority' (>= 50% must meet condition). Default option
%   is 'any'.
%
%   This function returns a vector (length = number of rows of initial)
%   of logical 1s if it meets all Turing conditions, or logical 0 if it
%   doesn't

    defaultMultipleHandling = 'matching';
    optionsMultipleHandling = {'any', 'all', 'majority', 'matching'};
    
    lengthVal = length(condition1(:, 1)); %Length that each condition must be
    lengthTest = @(x) (length(x(:, 1)) == lengthVal);
    
    p = inputParser;
    addRequired(p, 'condition1', lengthTest);
    addRequired(p, 'condition2', lengthTest);
    addRequired(p, 'condition3', lengthTest);
    addRequired(p, 'condition4', lengthTest);
    addOptional(p, 'MultipleHandling', defaultMultipleHandling, ...
        @(x) any(validatestring(x, optionsMultipleHandling)));
    parse(p, condition1, condition2, condition3, condition4, varargin{:});
    
    out = zeros(lengthVal, 1, 'logical');
    
    cond1 = p.Results.condition1;
    cond2 = p.Results.condition2;
    cond3 = p.Results.condition3;
    cond4 = p.Results.condition4;
    
    switch p.Results.MultipleHandling
        
        case 'any'
            for i = 1:lengthVal %Could vectorize this, maybe later
                if (any(cond1(i,:) > 0) && (any(cond2(i,:) > 0)) && ...
                        (any(cond3(i,:) < 0)) && (any(cond4(i,:) > 0)))
                    out(i) = 1;
                end
            end
            
        case 'all'
            for i = 1:lengthVal %Could vectorize this, maybe later
                if (all(cond1(i,:) > 0) && (all(cond2(i,:) > 0)) && ...
                        (all(cond3(i,:) < 0)) && (all(cond4(i,:) > 0)))
                    out(i) = 1;
                end
            end
            
        case 'majority'
            for i = 1:lengthVal %This is messy, may need to rethink
                prop1 = sum(cond1(i, :) > 0) / length(cond1(i, :));
                prop2 = sum(cond2(i, :) > 0) / length(cond2(i, :));
                prop3 = sum(cond3(i, :) < 0) / length(cond3(i, :));
                prop4 = sum(cond4(i, :) > 0) / length(cond4(i, :));
                
                if (prop1 >= 0.5) && (prop2 >= 0.5) && (prop3 >= 0.5) ...
                        && (prop4 >= 0.5)
                    out(i) = 1;
                end
            end
            
        case 'matching'
            for i = 1:lengthVal %I'm an idiot for not adding this before - definitely screwed up results
                if (cond1(i) > 0) && (cond2(i) > 0) %if steady state is stable
                    for d = 1:length(cond3(i, :))
                        if (cond3(i, d) < 0) && (cond4(i, d) > 0) %Gives system a shot regardless for each diffusion value
                            out(i) = 1;
                        end %If both conditions aren't matched for a single diffusion value, then it doesn't count
                    end
                end
            end
    end
    
    output = out;
        
end

