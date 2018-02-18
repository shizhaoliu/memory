function [MonteCarloPval,H] = SignificanceTest(R, shuffleR, level)
    % SignificanceTest returns Monte Carlo p value of an event.
    % If p value is smaller than level, H = 1.Else H = 0
    
    % Meaning of input:
    % -------------------------------------
    % R : Real value of Rw/Rwd of this event
    
    % shuffle R: All  Rw/Rwd values of this event
    
    % level: Significance level. Default is 0.05
    % -------------------------------------
    % Note that R and shuffleR needs to be absolute values when doing test
    % fo linear-weighted-correlation
    
    % ----------------------------------------------------
    % Author: Shizhao Liu(ashleyshizhaoliu@126.com)
    % Last modified: 2018/02/11
    % ----------------------------------------------------
    if nargin < 3
        level = 0.05;
    end
    
    index = find(shuffleR >= R);
    MonteCarloPval = (length(index)+1)/(length(shuffleR)+1);
    if MonteCarloPval<=level
        H = 1;
    else
        H = 0;
    end
