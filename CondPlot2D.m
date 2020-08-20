function s = CondPlot2D(x, y, data, varargin)
    %This function plots the parameter space, and color codes the system to
    %show where in parameter space the condition is met.
    %
    %The inputs x and y should be vectors containing the coordinates of
    %each scatter plot point. They should be the same length.
    %The data matrix (or vector) should have the same number of rows as x
    %and y, but can have any number of columns (will be able to choose any
    %or all values must meet condition).
    %
    %Optional Parameters:
    %
    %'AnyOrAll': string; determines whether all values for data point must 
    %   meet condition to qualify. Options are 'any', 'all', or 'majority';
    %   Default is 'any'.
    %   
    %'ConditionValue': scalar; determines value to compare data. Options 
    %   are any number; Default is 0.
    %
    %'ConditionGreaterThan': logical; determines whether to search for 
    %   data greater than (true) or less than (false). Default is true.
    %
    %'MetColor': triplet; sets color of data point if it meets condition.
    %   Options are 1x3 RGB vector (values from 0-1); Default is red ([1 0
    %   0]).
    %
    %'Title': string; sets title of plot. Default is 'Condition Plot'.
    %
    %'Xlabel': string; sets title of x axis. If 'none', does not include
    %   title. Default is 'none'.
    %
    %'Ylabel': string; sets title of y axis. If 'none', does not include
    %   title. Default is 'none'.
    
    defaultAnyOrAll = 'any';
    defaultConditionValue = 0;
    defaultConditionGreaterThan = true;
    defaultMetColor = [1 0 0];
    %defaultTitle = 'Condition Plot';
    %defaultXlabel = 'none';
    %defaultYlabel = 'none';
    possibleAnyOrAll = {'any', 'all', 'majority'};
    
    p = inputParser;
    validColor = @(x) isvector(x) && (length(x) == 3) && all(x <= 1) && all(x >= 0);
    addRequired(p, 'x', @isvector);
    addRequired(p, 'y', @isvector);
    addRequired(p, 'data', @isnumeric);
    addParameter(p, 'AnyOrAll', defaultAnyOrAll, @(x) any(validatestring(x, possibleAnyOrAll)));
    addParameter(p, 'ConditionValue', defaultConditionValue, @isscalar);
    addParameter(p, 'ConditionGreaterThan', defaultConditionGreaterThan, @islogical);
    addParameter(p, 'MetColor', defaultMetColor, validColor);
    %addParameter(p, 'Title', defaultTitle, @isstring);
    %addParameter(p, 'Xlabel', defaultXlabel, @isstring);
    %addParameter(p, 'Ylabel', defaultYlabel, @isstring);
    parse(p, x, y, data, varargin{:});
    
    xVal = p.Results.x;
    yVal = p.Results.y;
    dataMat = p.Results.data;
    
    if length(xVal) ~= length(yVal) || length(xVal) ~= length(dataMat(:, 1))
        ME = MException('ScatterPlot:nonequivalentSizes', ...
            'Nonequal number of entries');
        throw(ME)
    end
    
    %Process Data
    numCol = length(dataMat(1, :));
    condVal = p.Results.ConditionValue;
    colorOut = zeros(length(xVal), 3); %Coloring matrix for scatterplot
    switch p.Results.AnyOrAll
        case 'any'
            if p.Results.ConditionGreaterThan
                for i = 1:length(dataMat(:, 1))
                    if any(dataMat(i, :) > condVal)
                        colorOut(i, :) = p.Results.MetColor;
                    end
                end
            elseif ~p.Results.ConditionGreaterThan
                for i = 1:length(dataMat(:, 1))
                    if any(dataMat(i, :) < condVal)
                        colorOut(i, :) = p.Results.MetColor;
                    end
                end
            end
            
        case 'majority'
            if p.Results.ConditionGreaterThan
                for i = 1:length(dataMat(:, 1))
                    numMet = sum(dataMat(i, :) > condVal);
                    if (numMet/numCol) >= 0.5
                        colorOut(i, :) = p.Results.MetColor;
                    end
                end
            elseif ~p.Results.ConditionGreaterThan
                for i = 1:length(dataMat(:, 1))
                    numMet = sum(dataMat(i, :) < condVal);
                    if (numMet/numCol) >= 0.5
                        colorOut(i, :) = p.Results.MetColor;
                    end
                end
            end
            
        case 'all'
            if p.Results.ConditionGreaterThan
                for i = 1:length(dataMat(:, 1))
                    if all(dataMat(i, :) > condVal)
                        colorOut(i, :) = p.Results.MetColor;
                    end
                end
            elseif ~p.Results.ConditionGreaterThan
                for i = 1:length(dataMat(:, 1))
                    if all(dataMat(i, :) < condVal)
                        colorOut(i, :) = p.Results.MetColor;
                    end
                end
            end
    end
    
    %Now, create the scatterplot
    s = scatter(xVal, yVal, 8, colorOut, 'MarkerFaceColor', 'flat');
end
    
    
            
                    
    
    
    