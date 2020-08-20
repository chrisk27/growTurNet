function [C] = ClusterAnalysis(saver,n,k,f_ode,tspan)
%Note: This is copied from STAR Code (Scholes et al.)
%
%This is a function that clusters a group of suggested steady states
%(endpoints of several ODE simulations) into a set of non-redundant steady
%states. First each nodes values are analysed for their ditribution.
%Multiple each number of these peaks for each node results in the number
%(= count) of possible steady states the system can exhibit. The data is
%clustered according to the value count. From each centroid of the
%clustering an ODE simulation is repeated and a non-redundant list of
%steady state estimations derived.

%Input(s)
% - saver: set of suggested steady states from multiple ODE simulations
% - n: node number
% - k: parameters
% - f_ode: ode function
% - tspan: time span for ODE simulation

%Output(s)
% - C: non redundant list of steady state estimations. Rows define steady
% state sets, columns represent each node x(1)...x(n)


    max_system = max(saver(:));
    if max_system < 10000
        ymax = 10000; %only use if any value above threshold.
    else
        ymax = max_system;
    end
    
    count = 1;
    for p1 = 1:n
        %h_val = imgaussfilt(saver(:,p1),10); %data smoothing
        h = histcounts(saver(:,p1),'BinLimits',[0 ymax],'NumBins',100);
        [~,locs] = findpeaks(h);%(h,'MinPeakDistance',20,'MinPeakProminence',10); %minPeakdistance needs some adjustement...     
        if ~isempty(locs) 
            count = count*length(locs); %figure out number of clusters
        end
    end
    %Perform cluster analysis
    if count > length(saver(:,1))
        count = length(saver(:,1));
    end
    [~,C] = kmeans(saver,count);
    
    C = sortrows(C);

    
    %Simulate to refine SS position
    
    for p1 = 1:length(C(:,1))
        y_0 = C(p1,:);
        [~, xout] = ode15s(@(t, x)ode(t,x,f_ode,k),tspan*2,y_0);
        C(p1,:) = xout(length(xout(:,1)),:);
    end
     
    q=1;%m2 = 1;
    if count > 1
        true_cluster(q) = 1;
        for m2 = 1:count-1
           if abs(C(m2,:)-C(m2+1,:)) > C(m2,:)/100
              q = q+1;
              true_cluster(q) = m2+1;
           end
        end
        C = C(true_cluster,:);
    end
    
    true_cluster2 = [];
    q=1;
    for m3 = 1:length(C(:,1))
        if sum(abs(f_ode(C(m3,:)',k))) < 0.001
            true_cluster2(q) = m3;
            q = q+1;
        end
    end
    C = C(true_cluster2,:);

end