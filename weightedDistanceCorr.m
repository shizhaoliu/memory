function dcor = weightedDistanceCorr(data,dist_time,dist_matrix_states,options)
    % weightedDistanceCorr returns weighted-distance-correlation of an candidate event
    
    % Meaning of input
    % -----------------------------------
    % Data :  Posterior probability(time * state/position bin matrix) of the candidate event
    
    % dist_time : the distance matrix of every pair of time points of the candidate event(N_time bin * N_time bin matrix)
    
    % dist_states_matrix: distance matrix(number of position bin * number
    % of position bin) of all position bins
    
    % options: options set by user
    
  
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    %%%%%%%%%%%%%%
    dmax = options.dmax;
    % Nw states with largest posterior probability are used when calculating y-direction distance
    N_w = options.Nw;
    % Number of time bins
    T = size(dist_time,1);
    
    % Get the posterior probability of the Nw largest value at each t_bin
    weights = zeros(T,N_w);
    states = zeros(T,N_w);
    for t = 1:T
        [p,s] = sort(data(t,:),'descend');
        weights(t,:) = p(1:N_w) / sum(p(1:N_w));
        states(t,:) = s(1:N_w);
    end
    % Caculate the weighted-distance on y direction based on 
    % weights and distance_matrix of states
    dist_states = zeros(T,T);
    for i = 1:T
        for j = 1:T
            if i ~= j
                    w1 = weights(i,:);
                    w2 = weights(j,:);
                    s1 = states(i,:);
                    s2 = states(j,:);
                    dist_states(i,j) = 0;
                    for n = 1:N_w
                        for m = 1:N_w
                            dist_states(i,j)  = dist_states(i,j) + w1(n)*w2(m)*dist_matrix_states(s1(n),s2(m));
                        end
                    end
            end
        end
    end

% The following part is standard procudure for calculation distance correlation

% credits to Shen Liu
% Get from https://cn.mathworks.com/matlabcentral/fileexchange/39905-distance-correlation


a = dist_time;
mcol = mean(a);
mrow = mean(a,2);
ajbar = ones(size(mrow))*mcol;
akbar = mrow*ones(size(mcol));
abar = mean(mean(a))*ones(size(a));
A = a - ajbar - akbar + abar;

b = dist_states;
mcol = mean(b);
mrow = mean(b,2);
bjbar = ones(size(mrow))*mcol;
bkbar = mrow*ones(size(mcol));
bbar = mean(mean(b))*ones(size(b));
B = b - bjbar - bkbar + bbar;

% Calculate squared sample distance covariance and variances
dcov = sum(sum(A.*B))/(size(mrow,1)^2);

dvarx = sum(sum(A.*A))/(size(mrow,1)^2);
dvary = sum(sum(B.*B))/(size(mrow,1)^2);

% Calculate the distance correlation
dcor = sqrt(dcov/sqrt(dvarx*dvary));


end