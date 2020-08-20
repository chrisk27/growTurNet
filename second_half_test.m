%Now, I am writing another script to determine whether the steady state is
%stable or not in the absence of diffusion. This is much more complicated
%(and deviates much more from Scholes et al.'s STAR Code, since I am going
%to try and use the first two stability conditions from Madsvamuse et al.
%instead of the "standard" linear stability analysis. However, parts of
%this will indeed look like linear stability analysis, since I need to
%perturb around the steady state calculated in SteadyStateFinder.

disp('Analyzing Steady States.');
tic;

%The first portion is building up the functions that are required to 
X = sym('x', [1 n]); %Creates symbolic nodes
K = sym('k', [1 k_length]); %Creates symbolic parameter values
rxn_sym = [rxn_func{1}(X, K), rxn_func{2}(X, K)]; %Creates a symbolic function of the reaction component
rxn_jac = jacobian(rxn_sym, X); %Calculates Jacobian Matrix
rxn_jac_func = matlabFunction(rxn_jac); %Converts jacobian to symbolic function
rxn_jac_in = symvar(rxn_jac); % read out input parameters of jacobian Z
Xoverlap = ismember(X,rxn_jac_in); %find nodes that have to be substituted in Jacobian
Koverlap = ismember(K, rxn_jac_in); %find parameters that have to be substituted in Jacobian

count = 0; %Counts how many systems meet conditions for Turing patterns

for i = 1:length(SteadyState_save(:, 1))
    if System_save(SteadyState_save(i,1),k_length+1) ~= 1 && SteadyState_save(i,n+2) ~= 2 %Exclude oscillators and unstable solutions
        
        %Calculate appropriate timescales. The growth timescale, T_gr, is the
        %infinmum time that it takes for the domain size to increase by a
        %factor of e. The dynamics timescale, T_dyn, is the maximum of the 
        %diffusion timescale and the kinetics timescale. Since we are only
        %considering aspymtotically slow growth, we will use the parameter
        %epsilon = T_dyn / T_gr << 1. This should be (by Madzvamuse et al.)
        %between 5 x 10^(-4) and 0.12
        %
        %To start, we will calculate the dynamic timescale at the initial
        %condition. This is done by finding the largest real part of the
        %eigenvalue of the jacobian at the intial condition.
        kk = System_save(SteadyState_save(i,1),1:k_length);
        x_init = SteadyState_save(i,2:n+1);
        init_input = num2cell([kk(Koverlap), x_init(Xoverlap)]);
        eig_init = eig(rxn_jac_func(init_input{:})); %Calculate the eigenvalues
        max_eig_init = max(abs(real(eig_init))); %Find max real part of the eigenvalue
        T_dyn = 1 / (gamma*max_eig_init); %Dynamic Timescale
        syms ts %Define symbolic time to solve for growth timescale
        T_gr = eval(solve(exp(1) == grow_func(ts, kk), ts)); %Growth timescale
        if isempty(T_gr) %If there isn't a growth timescale (i.e. no growth)
            T_gr = inf;
        end
        
        %Next, we must calculate the midpoint time, t_mid. This is where we
        %will evaluate the concentration profiles of the two nodes, and
        %then check the stability.
        t_star = T_dyn /2; %Uses the midpoint approximation
        %Lay out the system of ODEs we need to solve
        nonDiffEq = @(x, t, k, m, gamma) [gamma * rxn_func{1}(x, k) ...
            - m * dilution_func(t, k); gamma * rxn_func{2}(x, k) ...
            - m * dilution_func(t, k)];
        grow_tspan = [0 t_star];
        try
            [t_nonDiff, conc_nonDiff] = ode15s(@(t, x)nonDiffEq(x, t, kk, m, gamma),...
                grow_tspan, x_init);
        catch
            [t_nonDiff, conc_nonDiff] = ode23s(@(t, x)nonDiffEq(x, t, kk, m, gamma),...
                grow_tspan, x_init);
        end
        x_tstar = conc_nonDiff(end, :);
        
        %Now, we will check the first two stability conditions
        h_tstar = dilution_func(t_star, kk);
        stability_jac_input = num2cell([x_tstar(Xoverlap), kk(Koverlap)]);
        stability_jac = rxn_jac_func(stability_jac_input{:});
        tr_jac = trace(stability_jac);
        det_jac = det(stability_jac);
        cond_1 = (-1) * gamma * tr_jac + 2*h_tstar;
        cond_2 = (-1) * h_tstar * gamma * tr_jac +...
            gamma^2 * det_jac;
        
        %Fill in the conditioned values
        SS_save_modified(i, n+3) = cond_1;
        SS_save_modified(i, n+4) = cond_2;
        if cond_1 <= 0 || cond_2 <= 0 % If homogeneous steady state is not stable
            SS_save_modified(i,n+2) = 1;
            Sys_save_modified(SteadyState_save(i,1),k_length+2) = System_save(SteadyState_save(i,1),k_length+2)-1;
            %SteadyState_save(i, n+5) = NaN;
            %SteadyState_save(i, n+6) = NaN;
        else %If it seems stable, continue on to test for Turing Instability
            %Note: I'm ignoring most of what the STAR code does for
            %creating different diffusion ratios, since I'm only looking at
            %two nodes, so it doesn't need to generalize to 3 or more. So,
            %a lot of the stuff that was generated in the definition step
            %is useless.
            %
            %Will always assume the diffusion constant for the first node
            %is one, and then loop through the other possibilities listed
            %in d_range 
            for j = 1:length(d_range)
                d_val = 1 / d_range(j);
                cond_3 = (-1) * gamma * (d_val*stability_jac(1, 1) + stability_jac(2, 2)) ...
                    + h_tstar * (1 + d_val);
                cond_4 = (h_tstar*(1+d_val) - gamma * (d_val*stability_jac(1, 1) + stability_jac(2, 2)))^2 ...
                    - 4*d_val*(gamma^2 * det_jac - gamma*h_tstar*tr_jac);
                %SteadyState_save(i, n+5) = cond_3;
                %SteadyState_save(i, n+6) = cond_4;
                
                %Now, check to see if conditions are met
                if cond_3 < 0 && cond_4 > 0 %If it does meet conditions
                    count = count + 1;
                    T_save_modified(count, 1) = i; %Saves index for steady state conditions
                    T_save_modified(count, 2) = 1;
                    T_save_modified(count, 3) = d_range(j); %Saves diffusion of node 1 (above) and node 2
                    T_save_modified(count, 4:7) = [cond_1, cond_2, cond_3, cond_4];
                end
            end
        end
    end
end

disp('Finished Turing analysis of steady states in: ')
toc

% Delete all empty entries in Turing_save
%Delete all zero entries
for i = 1:length(Turing_save(:,1))
    if Turing_save(i,1) ==0 
        j = i-1;
        break
    end
end

Turing_save = Turing_save(1:j,:); %Get rid of unused entries except for 1 

%Save outputs to a file named ID (see Systemdefinition_example.m)
save(strcat(savePath,ID),'Sys_save_modified','SS_save_modified','T_save_modified', 'savePath')