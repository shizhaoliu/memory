function [prob,inferred_index,weights] = loglld(spikeObs,t_bin,position_firingRate)
    % loglld returns the following values:
    % ------------------------------------
    % prob: Posterior probability(time * state matrix) of event
    
    % inferred_index: index of selected position bin with
    % maximum posterior probability at each time point (1 * Time vector)
    
    % weights: weights of the above selected states.(i.e. their posterior probability)
    % (Time * 1 vector)
    % ------------------------------------
    
    % Meaning of input:
    % ------------------------------------
    % spikeObs: observation(neuron * time matrix)of this candidate event
    
    % t_bin: time bin size
    
    % position_firingRate:lmbda of all neurons(neuron * position bin/state matrix)
    
    % ------------------------------------
    
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    N_timeBin = size(spikeObs,2);
    N_posBin = size(position_firingRate,2);
    loglld = zeros(N_timeBin,N_posBin);
    prob = zeros(N_timeBin,N_posBin);
    inferred_index = zeros(1,N_timeBin);
    weights = zeros(N_timeBin,1);

    for t = 1 : N_timeBin
        sp = spikeObs(:,t);
        lamda = position_firingRate;
        lamda = lamda + eps * ones(size(lamda));
        %log likelihood
        loglld(t,:) = sp' * log(lamda * t_bin)  - sum(lamda*t_bin,1)-sum(log(factorial(sp)))*ones(1,N_posBin);
        % posterior matrix
        prob(t,:) = exp(loglld(t,:)) / sum(exp(loglld(t,:)));
        %max prob index and its weight
        inferred_index(t) = find(prob(t,:) == max(prob(t,:)));
        weights(t) = prob(t,inferred_index(t)) / sum(prob(t,:));

    end
end