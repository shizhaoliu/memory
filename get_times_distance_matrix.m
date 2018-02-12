function dist_time = get_times_distance_matrix(T,options)
    % get_times_distance_matrix returns distance matrix of every pair of time points 
    % in the candidate event(N_time bin * N_time bin matrix)
    
    % Meaning of input:
    % -------------------------------
    % T: number of time bin in this candidate event
    
    % options: options set by user, some fields used in this function are explained:
    %   options.option_x = 'linear' :abs(xi-xj); (paper equation (2.6))
    %           In this circumstance, options.delta_t is not used
    %   options.option_x = 'linear_bonded':abs(xi-xj)with a max distance (paper equation (2.7))
    %           In this circumstance, if the difference of two time points is larger than
    %           options.delta_t, then the distance between them is limited by dmax
  
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    
    %%%%%%%%%%%%%%
    option_x = options.option_x;
    delta_t = options.delta_t;
    dmax = options.dmax;
    tx1 = strcmp(option_x,'linear');
    tx2 = strcmp(option_x,'linear_bonded');
    dist_time = zeros(T,T);
    %linear
    if tx1 == 1 && tx2 == 0
        for i = 1:T
            for j = 1:T
                dist_time(i,j) = abs(i-j);
            end
        end
    end
    %linear_bonded
    if tx1 == 0 && tx2 ==1
        for i  = 1 : T
            for j = 1:T
                if abs(i - j) >= delta_t
                    dist_time(i,j) = dmax;
                else
                    dist_time(i,j) = dmax * abs(i-j) / delta_t;
                end
            end
        end
end
end