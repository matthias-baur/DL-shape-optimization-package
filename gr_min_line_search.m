%% Gradient descent line search
% search for minimum of objective function along negative gradient d

function beta_opt = gr_min_line_search(a,b,d,param, mult_number, beta_interval)
    len_ab = length(a);
    
    sum_beta = [0, 0, 0, 0];    % stores intermediate objective function values 
    
    gr = (sqrt(5) + 1) / 2;
    beta = [beta_interval, 0, 0];   % stores interval and midpoints on which golden ratio search works
    beta(3:4) = [beta(2) - (beta(2) - beta(1)) / gr, beta(1) + (beta(2) - beta(1)) / gr];
    sum_beta_needed = [0,0,1,1];    % tells which midpoints in GR search have to be newly calculated
    
    % golden ratio search
    while abs(beta(2) - beta(1)) > param.beta_rel*beta(1) || beta(1) == 0
        % calculate new value for new betas       
        if sum_beta_needed(3) == 1
            a3 = a - beta(3)*d(1:len_ab);
            b3 = b - beta(3)*d(len_ab+1:end);
            [sum_beta(3), ~, mult3] = direct_problem(a3,b3,param, mult_number);
        end 
        if sum_beta_needed(4) == 1
            a4 = a - beta(4)*d(1:len_ab);
            b4 = b - beta(4)*d(len_ab+1:end);
            [sum_beta(4), ~, mult4] = direct_problem(a4,b4,param, mult_number);
        end
        
        % reassign beta interval and midpoints
        if sum_beta(3) < sum_beta(4) 
            beta(2) = beta(4);
            beta(4) = beta(3);
            sum_beta(4) = sum_beta(3);
            mult4 = mult3;
            beta(3) = beta(2) - (beta(2) - beta(1)) / gr;
            sum_beta_needed = [0,0,1,0];  
        else
            beta(1) = beta(3);
            beta(3) = beta(4);
            sum_beta(3) = sum_beta(4);
            mult3 = mult4;
            beta(4) = beta(1) + (beta(2) - beta(1)) / gr;
            sum_beta_needed = [0,0,0,1];  
        end
        beta
    end
    
    beta_opt = (beta(2) + beta(1)) / 2;
   
end

