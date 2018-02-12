function results = caculateCorrandZ_super(data,dist_states_matrix,position_firingRate,options)
    % caculateCorrandZ_super returns an struct contains results of all
    % candidtate event.
    % This struct contains the following field:
    % --------------------------------------------------
    % obs: Observation of event(neuron * time matrix),removed no-spike bins
    
    % prob: Posterior probability(time * state matrix) of event
    
    % Rw: linear-weighted-correlation value of this event
    
    % Rwd: weighted-distance-correlation value of this event
    
    % Z_w: Z-score of |Rw|(include three values: shuffle both neuron and time, only shuffle neuron and only shuffle time)
    
    % Z_wd : Z-score of Rwd(include three values: shuffle both neuron and time, only shuffle neuron and only shuffle time)
    
    % thrs_w: Threshold of |Rw| based on empirical distribution of all
    % shuffled |Rw|(also includes three values, same as above Z_wd)

    % thrs_wd: Threshold of Rwd based on empirical distribution of all
    % shuffled Rwd(also includes three values, same as above Z_wd)
    
    % Rw_null: all shuffled Rw values (include three cells, also same as above)
    
    % Rwd_null: all shuffled Rwd values (include three cells, also same as above)
    
    % act_cell:index of active neurons
    % -------------------------------------------------------
    
    % Meaning of input:
    % ------------------------------------------------------
    % data: an struct containing observation(neuron * time matrix)of all candidate events
    
    % dist_states_matrix: distance matrix(number of position bin * number
    % of position bin) of all position bins
    
    % position_firing_rate: lmbda of all neurons(neuron * position bin matrix)
    
    % options: options set by user
    % -----------------------------------------------------
    
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    t_bin = options.t_bin;
    m = 1;
    
    for k = 1 :size(data,2)
        display(['running:',num2str(100*k/size(data,2)),'%'])
        spikeObs = data(k).obs;
        % Find time bins without any spike and remove them
        zeroSpike_bins = find(sum(spikeObs,1) == 0);
        spikeObs(:,zeroSpike_bins) = [];

        N_timeBin = size(spikeObs,2); %length of event after removing no-spikes bins
        act_cell = find(sum(spikeObs,2) >= options.numActCell_criterion); 
        
        if length(act_cell) <=options.leastNumCell || N_timeBin <= options.leastNumTimeBin %if event too short or active neuron too few, skip this event
            continue
        end
   
        [prob,inferred_index,weights] = loglld(spikeObs,t_bin,position_firingRate); %loglikelihood
        dist_time = get_times_distance_matrix(N_timeBin,options); % get time-distance matrix
        weighted_corr = weightedcorrs([[1:1:N_timeBin];inferred_index]',weights); %caculate linear-weighted correlation
        % Real value of Rwd and |Rw| of this event
        Rw = weighted_corr(1,2);
        Rwd = weightedDistanceCorr(prob,dist_time,dist_states_matrix,options);%caculate weighted-distance correlation
        % Caculate Zscore, threshold and null distribution of Rwd and |Rw|
        [Rw_null,Rwd_null,Z_w,Z_wd,thrs_w,thrs_wd] = getZscore_super(spikeObs,Rw,Rwd,dist_states_matrix,position_firingRate,options);%caculate Zscore of Rw and Rwd
        % Make a struct containing all results 
        results(m).obs = spikeObs;
        results(m).prob = prob;
        results(m).Rw = Rw;
        results(m).Rwd = Rwd;
        results(m).Z_w = Z_w;
        results(m).Z_wd = Z_wd;
        results(m).thrs_w = thrs_w;
        results(m).thrs_wd = thrs_wd;
        results(m).Rw_null = Rw_null;
        results(m).Rwd_null = Rwd_null;
        results(m).act_cell = act_cell;

        
        m = m+1;
    end