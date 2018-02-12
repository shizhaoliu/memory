function [Rw_null,Rwd_null,Z_w,Z_wd,thrs_w,thrs_wd] = getZscore_super(spikeObs,weighted_corr,distance_corr,dist_states_matrix,position_firingRate,options)
    % getZscore_unsuper returns the following values:
    % -----------------------------------------------
    % Rw_null: all shuffled Rw values (include three cells,shuffle both neuron and time, only shuffle neuron and only shuffle time)
    
    % Rwd_null: all shuffled Rwd values(include three cells,same as above)
    
    % Z_w : Z-score of |Rw|(include three values, same as above)
    
    % Z_wd : Z-score of Rwd(include three values, same as above)
    
    % thrs_w: Threshold of |Rw| based on empirical distribution of shuffled
    % |Rw|(also includes three values, same as above)
    
    % thrs_wd: Threshold of Rwd based on empirical distribution of shuffled
    % Rwd(also includes three values, same as above)
    % -----------------------------------------------
    
    % Meaning of input:
    % ----------------------------------------------
    % spikeObs: observation(neuron * time matrix)of this candidate event
    
    % weighted_corr: linear-weighted-correlation value of this event
    
    % distanc_corr: weighted-distance-correlation value of this event
    
    % dist_states_matrix: distance matrix(number of states * number of states) of all hidden states
    
    % position_firing_rate: lmbda of all neurons(neuron * state matrix)
    
    % options: options set by user
    % ----------------------------------------------
    
  
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    t_bin = options.t_bin;
    N_timeBin = size(spikeObs,2);
    Rw_null = cell(3,1);
    Rwd_null = cell(3,1);
    %%
    %shuffle both neuron and time
    Rwd_shuffled = zeros(1000,1);
    Rw_shuffled = zeros(1000,1);


    for trial = 1:1000
        idx_row = randperm(size(spikeObs,1));
        shuffled_data1 = spikeObs(idx_row,:);
        idx_column = randperm(size(spikeObs,2));
        shuffled_data1 = shuffled_data1(:,idx_column);
        %caculate the weighted_correlation and distance_correlation of
        %shuffled data
        [prob,inferred_index,weights] = loglld(shuffled_data1,t_bin,position_firingRate);
        dist_time = get_times_distance_matrix(N_timeBin,options);
        Rw = weightedcorrs([[1:1:N_timeBin];inferred_index]',weights);
        Rw_shuffled(trial) = Rw(1,2);
        Rwd_shuffled(trial) = weightedDistanceCorr(prob,dist_time,dist_states_matrix,options);


    end

    Rw_null{1} = Rw_shuffled;
    Rwd_null{1} = Rwd_shuffled;
    % Z score of real Rwd/|Rw| of this event
    [~,~,~,Z_w(1)] = ztest(abs(weighted_corr),mean(abs(Rw_shuffled)),std(abs(Rw_shuffled)));
    [~,~,~,Z_wd(1)] = ztest(distance_corr,mean(Rwd_shuffled),std(Rwd_shuffled));

    % Threshold of Rwd/|Rw|(>=95% of the shuffled Rwd/|Rw|)
    % If Rwd/|Rw| of event is larger than the threshold, we consider this event
    % significant
    Rw_shuffled = sort(abs(Rw_shuffled),'ascend');
    thrs_w(1) = Rw_shuffled(round(0.95*length(Rw_shuffled)));
    Rwd_shuffled = sort(Rwd_shuffled,'ascend');
    thrs_wd(1) = Rwd_shuffled(round(0.95*length(Rwd_shuffled)));


    %%
    %only shuffle neuron
    Rwd_shuffled_cell = zeros(1000,1);
    Rw_shuffled_cell = zeros(1000,1);


    for trial = 1:1000
        idx_row = randperm(size(spikeObs,1));
        shuffled_data2 = spikeObs(idx_row,:);

        %caculate the weighted_correlation and distance_correlation of
        %shuffled data
        [prob,inferred_index,weights] = loglld(shuffled_data2,t_bin,position_firingRate);
        dist_time = get_times_distance_matrix(N_timeBin,options);
        Rw = weightedcorrs([[1:1:N_timeBin];inferred_index]',weights);
        Rw_shuffled_cell(trial) = Rw(1,2);
        Rwd_shuffled_cell(trial) = weightedDistanceCorr(prob,dist_time,dist_states_matrix,options);


    end

    Rw_null{2} = Rw_shuffled_cell;
    Rwd_null{2} = Rwd_shuffled_cell;
    % Z score of real Rwd/|Rw| of this event
    [~,~,~,Z_w(2)] = ztest(abs(weighted_corr),mean(abs(Rw_shuffled_cell)),std(abs(Rw_shuffled_cell)));
    [~,~,~,Z_wd(2)] = ztest(distance_corr,mean(Rwd_shuffled_cell),std(Rwd_shuffled_cell));

    % Threshold of Rwd/|Rw|(>=95% of the shuffled Rwd/|Rw|)
    % If Rwd/|Rw| of event is larger than the threshold, we consider this event
    % significant
    Rw_shuffled_cell = sort(abs(Rw_shuffled_cell),'ascend');
    thrs_w(2) = Rw_shuffled_cell(round(0.95*length(Rw_shuffled_cell)));
    Rwd_shuffled_cell = sort(Rwd_shuffled_cell,'ascend');
    thrs_wd(2) = Rwd_shuffled_cell(round(0.95*length(Rwd_shuffled_cell)));


    %%
    %only shuffle the time
    Rwd_shuffled_time = zeros(1000,1);
    Rw_shuffled_time = zeros(1000,1);

    for trial = 1:1000
        idx_column = randperm(size(spikeObs,2));
        shuffled_data3 =spikeObs(:,idx_column);
        %caculate the weighted_correlation and distance_correlation of
        %shuffled data
        [prob,inferred_index,weights] = loglld(shuffled_data3,t_bin,position_firingRate);
        dist_time = get_times_distance_matrix(N_timeBin,options);
        Rw = weightedcorrs([[1:1:N_timeBin];inferred_index]',weights);
        Rw_shuffled_time(trial) = Rw(1,2);
        Rwd_shuffled_time(trial) = weightedDistanceCorr(prob,dist_time,dist_states_matrix,options);


    end

    Rw_null{3} = Rw_shuffled_time;
    Rwd_null{3} = Rwd_shuffled_time;
    % Z score of real Rwd/|Rw| of this event
    [~,~,~,Z_w(3)] = ztest(abs(weighted_corr),mean(abs(Rw_shuffled_time)),std(abs(Rw_shuffled_time)));
    [~,~,~,Z_wd(3)] = ztest(distance_corr,mean(Rwd_shuffled_time),std(Rwd_shuffled_time));

    % Threshold of Rwd/|Rw|(>=95% of the shuffled Rwd/|Rw|)
    % If Rwd/|Rw| of event is larger than the threshold, we consider this event
    % significant
    Rw_shuffled_time = sort(abs(Rw_shuffled_time),'ascend');
    thrs_w(3) = Rw_shuffled_time(round(0.95*length(Rw_shuffled_time)));
    Rwd_shuffled_time = sort(Rwd_shuffled_time,'ascend');
    thrs_wd(3) = Rwd_shuffled_time(round(0.95*length(Rwd_shuffled_time)));

end


