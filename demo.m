clear all
clc
close all
%%
load('data/position_firingRate_unsuper'); 
load('data/trans_matrix_unsuper');
load('data/state_idx')
load('data/position_firingRate_super'); 
load('data/trans_matrix_super');
load('data/example_data')
%%
options.t_bin = 0.02;%for sleep time bin is 20ms
%set options of caculating distance matrix
options.dmax = 1.6;%max distance
options.alpha = 8;  %exponential's decay speed, won't be used when option_y = 'linear_bonded'
options.delta_t = 3; %bond of time.if |t1-t2| > delta_t, d_time = dmax
options.Nw = 2;%Nw prob largest points are used in weighted distance correlation 
% option of x/y distance metric
options.option_x = 'linear_bonded';
options.option_y = 'linear_bonded';

%a neuron is considered active if it has at least 2 spike in a candidate event
options.numActCell_criterion = 2;
% To be analyzied, an event needs contain at least 7 active neurons and it's
% number of time bins is larger than 4
options.leastNumCell = 7;
options.leastNumTimeBin = 4;

%get distance between states based on transition matrix
dist_states_matrix_unsuper = get_states_distance_matrix(trans_matrix_unsuper,options);
dist_states_matrix_super = get_states_distance_matrix(trans_matrix_super,options);
%%
% run the experiment
example_results_unsuper =  caculateCorrandZ_unsuper(example_data,dist_states_matrix_unsuper,position_firingRate_unsuper,options);
example_results_super =  caculateCorrandZ_super(example_data,dist_states_matrix_super,position_firingRate_super,options);
%%
% Monte Carlo p-value of example data
[p_super_Rw(1,1),~] = SignificanceTest(abs(example_results_super(1).Rw), abs(example_results_super(1).Rw_null{1}));
[p_super_Rw(2,1),~] = SignificanceTest(abs(example_results_super(2).Rw), abs(example_results_super(2).Rw_null{1}));
[p_super_Rwd(1,1),~] = SignificanceTest(example_results_super(1).Rwd, example_results_super(1).Rwd_null{1});
[p_super_Rwd(2,1),~] = SignificanceTest(example_results_super(2).Rwd, example_results_super(2).Rwd_null{1});

[p_unsuper_Rwd(1,1),~] = SignificanceTest(example_results_unsuper(1).Rwd, example_results_unsuper(1).Rwd_null{1});
[p_unsuper_Rwd(2,1),~] = SignificanceTest(example_results_unsuper(2).Rwd, example_results_unsuper(2).Rwd_null{1});
T = table(p_super_Rw,p_super_Rwd,p_unsuper_Rwd,'Rownames',{'Example1';'Example2'})
%%
% Showing the figures
make_figure_example(example_results_super,example_results_unsuper,state_idx)
