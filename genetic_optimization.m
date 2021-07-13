%% Genetic optimization
% optimizes domains with selection and reproduction rules

% set parameters for genetic algorithm
N_best = 50;    % no. of fittest domains to keep in each step
                % also no. of newly added random domains in each step
N_pop  = 500;   % population size
N_iter = 10;    % number of iterations

% initialization of domain population
A_pop = zeros(N_pop,len_ab);
B_pop = zeros(N_pop,len_ab);
A_pop(:,1) = 1/sqrt(pi);

sum_pop = zeros(N_pop,1);

mult_tol = 0.1;
sigma = 0.1;    % base variance of coefficients
weight_vector = 1.1*exp(-0.3*(0:len_ab-2));     % modifies variance of coefficients
c_cutoff = 0.35;         % upper limit for coefficients
A_temp = bsxfun(@times,normrnd(0,sigma,N_pop,len_ab-1), weight_vector);     % normal distributed coefficients
B_temp = bsxfun(@times,normrnd(0,sigma,N_pop,len_ab-1), weight_vector);     % variance smaller for higher frequencies
A_rand_numbers = 0.1*rand(size(A_temp))-0.05;
B_rand_numbers = 0.1*rand(size(B_temp))-0.05;
A_temp(abs(A_pop(:,2:end))>c_cutoff) = A_rand_numbers(abs(A_pop(:,2:end))>c_cutoff);    % replace too large coefficients with
B_temp(abs(B_pop(:,2:end))>c_cutoff) = B_rand_numbers(abs(B_pop(:,2:end))>c_cutoff);    % smaller random number
A_pop(:,2:end) = A_temp;
B_pop(:,2:end) = B_temp;

A_pop_new = A_pop;
B_pop_new = B_pop;
sum_pop_new = sum_pop;

% evaluate fitness for 1:N_best domains
for m = 1:N_best
        [r, ~] = define_r(A_pop(m,:),B_pop(m,:));
        [a,b] = rescale_domain(A_pop(m,:),B_pop(m,:),r,t_col,dt_col);  
        A_pop(m,:)=a;
        B_pop(m,:)=b;
        [sum_n, ~, ~] = direct_problem(a,b,param,sum_length);
        sum_pop(m) = sum_n;
end

%% genetic algorithm loop

for j = 1:N_iter-1
    tic
    % evaluate fitness for N_best+1:N_pop domains
    for m = N_best+1:N_pop
        [r, ~] = define_r(A_pop(m,:),B_pop(m,:));
        [a,b] = rescale_domain(A_pop(m,:),B_pop(m,:),r,t_col,dt_col);  
        A_pop(m,:)=a;
        B_pop(m,:)=b;
        [sum_n, ~, ~] = direct_problem(a,b,param,sum_length);
        sum_pop(m) = sum_n;
    end

    % Selection
    [sorted_k_pop, sorted_indices] = sort(sum_pop);
    A_pop_new(1:N_best,:) = A_pop(sorted_indices(1:N_best),:);
    B_pop_new(1:N_best,:) = B_pop(sorted_indices(1:N_best),:);
    sum_pop_new(1:N_best) = sum_pop(sorted_indices(1:N_best));
    disp("Top 5 fitness values:")  
    disp(sum_pop_new(1:5))

    % Reproduction
    for k = 1:(N_pop-2*N_best)
        p_1 = randi(N_best); % index parent 1
        p_2 = randi(N_best); % index parent 2
        t = rand;   % interpolation parameter
        A_pop_new(N_best+k,2:end) = t*A_pop(p_1,2:end)+(1-t)*A_pop(p_2,2:end);
        B_pop_new(N_best+k,2:end) = t*B_pop(p_1,2:end)+(1-t)*B_pop(p_2,2:end);
    end

    % Add random domains
    A_temp = bsxfun(@times,normrnd(0,sigma,N_best,len_ab-1), weight_vector);
    B_temp = bsxfun(@times,normrnd(0,sigma,N_best,len_ab-1), weight_vector);
    A_rand_numbers = 0.1*rand(size(A_temp))-0.05;
    B_rand_numbers = 0.1*rand(size(B_temp))-0.05;
    A_temp(abs(A_pop(:,2:end))>c_cutoff) = A_rand_numbers(abs(A_pop(:,2:end))>c_cutoff);
    B_temp(abs(B_pop(:,2:end))>c_cutoff) = B_rand_numbers(abs(B_pop(:,2:end))>c_cutoff);
    A_pop_new(end-N_best+1:end,2:end) = A_temp;
    B_pop_new(end-N_best+1:end,2:end) = B_temp;

    % Take new population
    A_pop = A_pop_new;
    B_pop = B_pop_new;
    sum_pop = sum_pop_new;
    
    toc
    disp("iteration j = " + num2str(j) + " has finished")
end

% evaluate fitness for N_best+1:N_pop domains a last time
for m = N_best+1:N_pop
        [r, ~] = define_r(A_pop(m,:),B_pop(m,:));
        [a,b] = rescale_domain(A_pop(m,:),B_pop(m,:),r,t_col,dt_col);  
        A_pop(m,:)=a;
        B_pop(m,:)=b;
        [sum_n, ~, ~] = direct_problem(a,b,param,sum_length);
        sum_pop(m) = sum_n;
end

% final sorted population
[sorted_k_pop, sorted_indices] = sort(sum_pop);

A_pop = A_pop(sorted_indices,:);
B_pop = B_pop(sorted_indices,:);
sum_pop = sum_pop(sorted_indices);


