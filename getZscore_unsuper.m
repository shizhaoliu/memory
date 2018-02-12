function [Rwd_null,Z_wd,thrs_wd] = getZscore_unsuper(spikeObs,distance_corr,dist_states_matrix,position_firingRate,options)

    % getZscore_unsuper returns the following values:
    % -----------------------------------------------
    % Rwd_null: all shuffled Rwd values (include three cells, also same as above)
    
    % Z_wd : Z-score of Rwd(include three values: shuffle both neuron and time, only shuffle neuron and only shuffle time)
    
    % thrs_wd: Threshold of Rwd based on empirical distribution of shuffled
    % Rwd(also includes three values, same as above Z_wd)
    % -----------------------------------------------
    
    % Meaning of input:
    % ----------------------------------------------
    % spikeObs: observation(neuron * time matrix)of this candidate event
    
    % distanc_corr: weighted-distance-correlation value of this event
    
    % dist_states_matrix: distance matrix(number of states * number of states) of all hidden states
    
    % position_firing_rate: lmbda of all neurons(neuron * state matrix)
    
    % options: options set by user
    % ----------------------------------------------
    
  
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    
    Nw = options.Nw;
    t_bin = options.t_bin;
    N_timeBin = size(spikeObs,2);
    Rwd_null = cell(3,1);
    %%
    % Shuffle both neuron and time
    Rwd_shuffled = zeros(1000,1);

    for trial = 1:1000
        % Shuffle neurons
        idx_row = randperm(size(spikeObs,1));
        shuffled_data1 = spikeObs(idx_row,:);
        % Shuffle time
        idx_column = randperm(size(spikeObs,2));
        shuffled_data1 = shuffled_data1(:,idx_column);
        % Caculate the weighted_distance_correlation of shuffled data
        [prob,~,~] = loglld(shuffled_data1,t_bin,position_firingRate);
        dist_time = get_times_distance_matrix(N_timeBin,options);
        Rwd_shuffled(trial) = weightedDistanceCorr(prob,dist_time,dist_states_matrix,options);


    end

    Rwd_null{1} = Rwd_shuffled;
    % Z score of real Rwd of this event
    [~,~,~,Z_wd(1)] = ztest(distance_corr,mean(Rwd_shuffled),std(Rwd_shuffled));
    % Threshold of Rwd(>=95% of the shuffled Rwd)
    % If Rwd of event is larger than the threshold, we consider this event
    % significant
    Rwd_shuffled = sort(Rwd_shuffled,'ascend');
    thrs_wd(1) = Rwd_shuffled(round(0.95*length(Rwd_shuffled)));


    %%
    %only shuffle neuron
    Rwd_shuffled_cell = zeros(1000,1);


    for trial = 1:1000
        idx_row = randperm(size(spikeObs,1));
        shuffled_data2 = spikeObs(idx_row,:);

        %caculate the weighted_correlation and distance_correlation of
        %shuffled data
        [prob,~,~] = loglld(shuffled_data2,t_bin,position_firingRate);
        dist_time = get_times_distance_matrix(N_timeBin,options);

        Rwd_shuffled_cell(trial) = weightedDistanceCorr(prob,dist_time,dist_states_matrix,options);


    end

    Rwd_null{2} = Rwd_shuffled_cell;
    % Z score of real Rwd of this event
    [~,~,~,Z_wd(2)] = ztest(distance_corr,mean(Rwd_shuffled_cell),std(Rwd_shuffled_cell));
    % Threshold of Rwd(>=95% of the shuffled Rwd)
    % If Rwd of event is larger than the threshold, we consider this event
    % significant
    Rwd_shuffled_cell = sort(Rwd_shuffled_cell,'ascend');
    thrs_wd(2) = Rwd_shuffled_cell(round(0.95*length(Rwd_shuffled_cell)));


    %%
    %only shuffle time
    Rwd_shuffled_time = zeros(1000,1);

    for trial = 1:1000
        idx_column = randperm(size(spikeObs,2));
        shuffled_data3 =spikeObs(:,idx_column);
        %caculate the weighted_correlation and distance_correlation of
        %shuffled data
        [prob,~,~] = loglld(shuffled_data3,t_bin,position_firingRate);
        dist_time = get_times_distance_matrix(N_timeBin,options);

        Rwd_shuffled_time(trial) = weightedDistanceCorr(prob,dist_time,dist_states_matrix,options);

    end

    Rwd_null{3} = Rwd_shuffled_time;
    % Z score of real Rwd of this event
    [~,~,~,Z_wd(3)] = ztest(distance_corr,mean(Rwd_shuffled_time),std(Rwd_shuffled_time));
    % Threshold of Rwd(>=95% of the shuffled Rwd)
    % If Rwd of event is larger than the threshold, we consider this event
    % significant
    Rwd_shuffled_time = sort(Rwd_shuffled_time,'ascend');
    thrs_wd(3) = Rwd_shuffled_time(round(0.95*length(Rwd_shuffled_time)));

end

