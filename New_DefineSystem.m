%% Define Functions and Simulation
ID = 'Net11_Comp_ReallyMake2';

%The reaction function is the "reaction" term of the PDE, including inflow
%and decay terms. It will be used to determine the steady state.
rxn_func = {@(x,k)(k(4) - k(5)*x(1) + k(1)*((x(1)/k(2))^2 + (x(2)/k(3))^2) *...
    (1 + (x(1)/k(2))^2 + (x(2)/k(3))^2)^(-1));
     @(x,k)(k(8) - k(9)*x(2) + k(6)*(x(1)/k(7))^2 / (1 + (x(1)/k(7))^2))};
 
%Specify your growth function
%Right now, I'm actually just going to be scanning values for h(t), so this
%won't be used for a moment
grow_func = @(t,r)exp(r*t); %Type in growth function (phi(t) in Klika et al.)
    %Typical functions include:
        %Linear Growth: phi(t) = r*t + 1
        %Exponential Growth: phi(t) = exp(r*t)
dilution_func=@(t,r)r; %Type in dilution function (h(t) in Klika et al.)
    %Should choose same type of growth as grow_func. Typical functions include:
        %Linear Growth: h(t) = (m*r) / (r*t + 1)
        %Exponential Growth: h(t) = m*r
        
%% Fill In Parameters
k_length = 9; %Number of reaction parameters.
    %Note: This does NOT include the growth "r" parameter, unlike previous
    %iterations of this code
k_input = cell(1, k_length); %Where we will store parameter values to loop over

%Below this is where we define reaction parameter values
k_input{1} = 100; %V1
k_input{2} = 10; %k11 (x1 autocatalysis)
k_input{3} = 0.1; %k12 (x2 inhibiting x1)
k_input{4} = 0.1; %b1 (natural production of x1)
k_input{5} = logspace(-2, 3, 100); %mu1 (decay rate of x1)
k_input{6} = 100; %V2
%k_input{7} = 10; %k22 (x2 inhibiting x2)
k_input{7} = 100; %k21 (x1 producing x2)
k_input{8} = 0.1; %b2
k_input{9} = logspace(-2, 3, 100); %mu2

k_grid = combvec(k_input{:});

%Below this is the diffusion constant value. Node 1 will always have a
%diffusion constant of 1, and node 2 will vary in this range.
d_range = logspace(-3, 3, 7);

%Below this is the range in which the h value will vary. For exponential
%growth, it will be between 0 and r (which should max out at about 0.3),
%and for linear growth it should be between 0 and 1 
h_range = linspace(0, 1, 21);

%% Other Misc Parameters (usually won't change)
gamma = 1; %Scales influence of reaction on PDE
n = 2; %Number of nodes in Network
x_max = 20001;
int = 10000;
c_ini = permn(1:int:x_max,n); %Grid from which algorithm will sample. Here 3 initial condition