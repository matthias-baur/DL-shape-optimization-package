%% Set parameters
format compact 
format long

% Define sum that is used as objective function J:
% J(lambda_1,...,lambda_n) = lambda_{n-sum_length+1} + ... + lambda_n
n = 2;                  % index of largest eigenvalue involved
sum_length = n;         % length of sum

interval_end = 7;       % upper limit for eigenfrequencies (choose slightly 
                        % larger than sqrt(lambda_n) of initial domain for gradient descent)

len_ab = 20;            % number of sine and cosine cofficients

% Set basic parameters for collocation points 
N_col  = 2^8;           % number of collocation points 
dt_col  = 2*pi/N_col;   % step width of angular paramter
t_col = 0:dt_col:(2*pi-dt_col); 
p_beta = 0.03;          % displacement distance for MFS eigenfunction source points
                        % along normals, corresponds to \delta_y in the report
                        % recommendations: 
                        % genetic algorithm     p_beta = 0.01;
                        % gradient descent      p_beta = 0.03;

% gradient descent parameters
eps = 10^(-4);          % termination condition: exit gradient descent loop 
                        % if J changes less than eps in one step
n_steps = 0;            % counter for steps in gradient descent loop
max_steps = 300;        % number of steps after which loop is terminated
beta_opt = 1;           % gradient step width parameter initialization
beta_rel = 0.01;        % gradient step relative accuracy
mult_tol = 0.1;         % multiplicity tolerance for n-th eigenvalue, 
                        % all eigenvalues in [lambda_n-mult_tol , lambda_n] 
                        % are interpreted as one eigenvalue with multiplicity > 1,
                        % so far only relevant when sum_length = 1

% create struct of all input parameters
param = parameter_struct(n,sum_length,interval_end,len_ab,N_col,dt_col,t_col,p_beta,beta_opt,beta_rel,mult_tol);
