function dist_matrix = get_states_distance_matrix(transition_matrix,options)
    % dist_matrix returns distance matrix of all hidden states/position bin
    % Meaning of input:
    % ------------------------------------------------------------
    % transition_matrix: transition matrix of all states
    % options: options set by user,some fields used in this function are explained:
    %   options.option_y = 'linear':  abs(yi-yj); (paper equation (2.8))
    %   options.option_y = 'linear_bonded': distance based on transition matrix;l inear decay(paper  equation (2.9))
    %   options.option_y = 'exponential': distance based on transition matrix; exponential decay(paper  equation (2.10))
    %       Only in this circumstance, options.alpha, which controls decaying speed, is used
    % ------------------------------------------------------------
  
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    %%%%%%%%%%%%
    dmax = options.dmax;
    alpha = options.alpha;
    option_y = options.option_y;
    
    ty1 = strcmp(option_y,'linear');
    ty2 = strcmp(option_y,'exponential');
    ty3 = strcmp(option_y,'linear_bonded');
    N = size(transition_matrix,1);
    
    dist_matrix = zeros(N,N);
    
    %let diagnal elements of transition matrix = 0
    for i = 1:N
        transition_matrix(i,i)=0;
        transition_matrix(i,:) = transition_matrix(i,:) / sum(transition_matrix(i,:));
    end
    % get pmax after removing diagnal values
    pmax = max(transition_matrix(:));
    
    if ty1 ==0 && ty2 ==1 && ty3 == 0
        % options.option_y = 'exponential'
        for i = 1:N
            for j = 1:N
                if i == j
                    dist_matrix(i,j) = 0;
                else
                    p = (transition_matrix(i,j) + transition_matrix(j,i)) / 2;
                    dist_matrix(i,j) = dmax * exp(alpha*p/(p-1));
                end
            end
        end
    elseif ty1 ==1 && ty2 ==0 && ty3 ==0
        % options.option_y = 'linear'
        for i = 1:N
            for j = 1:N
                    dist_matrix(i,j) = abs(i-j);
                
            end
        end
    elseif ty1 ==0 && ty2 ==0 && ty3 ==1
        % options.option_y = 'linear_bonded'
        for i = 1:N
            for j = 1:N
                if i == j
                    dist_matrix(i,j) = 0;
                else
                    p = (transition_matrix(i,j) + transition_matrix(j,i)) / 2;
                    if p <= pmax
                        dist_matrix(i,j) = dmax*(1-p/(2*pmax));
                    else
                        dist_matrix(i,j) = 0;
                    end
                end
            end
        end
    end
end