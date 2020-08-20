%First define the system of interest:
n=2; %n specifies the node number (Should always be 2 for now)
%Define name for saving data:       
%ID=strcat('example_',datestr(now,30)); %name of file that will be saved.
ID = 'Net16_Exp_NonComp';

%Define some original parameters that will not be scanned
gamma = 1; %Parameter for the time scale of the reaction, will almost always just be 1
m = 1; %Number of spatial dimensions, will almost always just use 1

%Specify your reaction functions:
rxn_func={@(x,k)(k(4) - k(5)*x(1) + k(1) * (1 + (k(2)/x(1))^2)^(-1) * ...
    (1 + (k(3)/x(2))^2)^(-1));
     @(x,k)(k(9) - k(10)*x(2) + k(6) * (1 + (x(2)/k(7))^2)^(-1) *...
     (1 + (x(1)/k(8))^2)^(-1))};

%Specify your growth function
grow_func = @(t,k)exp(k(11)*t); %Type in growth function (phi(t) in Klika et al.)
    %Typical functions include:
        %Linear Growth: phi(t) = r*t + 1
        %Exponential Growth: phi(t) = exp(r*t)
dilution_func=@(t,k)k(11); %Type in dilution function (h(t) in Klika et al.)
    %Should choose same type of growth as grow_func. Typical functions include:
        %Linear Growth: h(t) = (m*r) / (r*t + 1)
        %Exponential Growth: h(t) = m*r
%__________________________________________________________________________
%Define diffusion parameters: the algorithm requires four inputs for this
%m_min = 2; %m_min specifies the minimum number of nodes that can diffuse
%m_max = n; %m_max specifies the maximum number of nodes that can diffuse
%binary_diffusor = permn([0 1],n); %Define which nodes diffuse (1 specifying diffusing, 0 non-diffusing). 
%Here we take all possible permutations of nodes HOWEVER only systems that fulfil the m_min and m_max requirement will be sampled.
%If only specific nodes should diffuse set binary_diffurso to e.g. [1 1 0] for a three node system (A B C) where A and B diffuse. 
d_range = logspace(-3,3,7); %This is the diffusion range that will be analysed.

%Note: always adjust m_min and m_max according to the node number that
%should diffuse even if binary_diffuser consists of only one element. 

%__________________________________________________________________________
%Define parameter space for sampling as k_grid
%Example:
k_length = 11; %Number of parameters to be sampled
    %First are the parameters for the reaction functions
    %The last parameter is the growth constant r
ks = logspace(-1,2,3); %define the range and interval parameters should be sampled at
k_input = cell(1, k_length);
for i = 1:(k_length - 1) %Fills all reaction components
    k_input{i} = ks;
end
mu_params = [0.01 0.1 1];
k_input{5} = mu_params; %Node 1 decay
k_input{k_length - 1} = mu_params; %Node 2 decay
k_input{k_length} = [0.1 0.2 0.3 0.5]; %These are the r values for the growth/dilution
k_grid = combvec(k_input{:}); %Create a combinatorial matrix from the above defined ranges.

%Please note: k_grid should not exceed 500000 different combinations.

%__________________________________________________________________________
%Define initial conditions for ODE
x_max = 20001;
int = 10000;
c_ini = permn(1:int:x_max,n); %Grid from which algorithm will sample. Here 3 initial conditions per node are sampled.

