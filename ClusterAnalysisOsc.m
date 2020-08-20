function [C] = ClusterAnalysisOsc(saver,n,k,f_ode,options)
%Note: This is copied from STAR method (Scholes et al.)
%
%The function analyses a set of supposed steady state values using the kmeans
%algorithm. First each node's steady state values are analysed seperately to
%find how many values are likely to be non-redundant entries. This is done by
%counting the number of local maxima in a histogram plot. The number of
%peaks for each node are then multiplied (=count) and the data clustered
%according to this count. fsolve is finally used to check if steady states
%are true steady states and redundant entries are hereafter removed from the 
%list of possible steady states. 

%Input(s):
% - saver: matrix with all end points of several ODE simulations
% - n: number of nodes
% - k: parameter values
% - f_ode: ODE function
% - options: fsolve solver options (see Matlab documentation for details)

%Output(s):
%- C: final suggested steady state(s) in form of a matrix. Rows define
%values for each node x(1)...x(n); each row represents a single steady
%state.

 max_system = max(saver(:)); %find maximum value within dataset
 %Make sure ymax >=1 10000
    if max_system < 10000
        ymax = 10000;
    else
        ymax = max_system;
    end  
    count = 1;
    
 %Calculate possible sets of steady states by multiplying peak counts for each species   
    for p1 = 1:n
        h = histcounts(saver(:,p1),'BinLimits',[0 ymax],'NumBins',100);
        [~,locs] = findpeaks(h);%,'MinPeakDistance',20,'MinPeakProminence',10); %minPeakdistance needs some adjustement...     
        if length(locs)>0 
            count = count*length(locs); %figure out number of clusters
        end
    end
    
    if count > length(saver(:,1))
        count = length(saver(:,1));
    end
    
 %Perform cluster analysis
    [~,C] = kmeans(saver,count);
    C = sortrows(C);

 %Using fsolve find "true" location of steady states (f_ode == 0)
    for p1 = 1:length(C(:,1))
        y_0 = C(p1,:)';
        C(p1,:)=fsolve(@(x)f_ode(x, k),y_0,options)'; %use fsolve for oscillators
    end
    
 %Remove redundant entries   
    q=1;m2 = 1;
    if count > 1
        true_cluster(q) = 1;
        for m2 = 1:count-1
           if abs(C(m2,:)-C(m2+1,:)) > C(m2,:)./100
              q = q+1;
              true_cluster(q) = m2+1;
           end
        end
        C = C(true_cluster,:);
    end
    q = 1;
    true_cluster2 = [];
    for i = 1:length(C(:,1))
        if all(C(i,:)> 0)
            true_cluster2(q) = i;
            q=q+1;
        end
    end 
    C = C(true_cluster2,:);
    
        true_cluster3 = [];
    q=1;
    for m3 = 1:length(C(:,1))
        if sum(abs(f_ode(C(m3,:)',k))) < 0.001
            true_cluster3(q) = m3;
            q = q+1;
        end
    end
    C = C(true_cluster3,:);
    
    if isempty(true_cluster3) == 1
        [~, xout] = ode15s(@(t, x)ode(t,x,f_ode,k),[0 2000],saver(1,:));
        if sum(abs(f_ode(xout(length(xout(:,1)),:),k))) < 0.001
            C = xout(length(xout(:,1)),:);
        end
   end    
end
