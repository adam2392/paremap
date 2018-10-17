%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script scalePlots.m
%
% Description: To rescale clustering plots to match each other. Makes
% comparing seizure more doable.
%
% Input:    
%
% Output:   no output returned.
%
%
% NOTE: 
%
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 08/21/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create an interval
%% For plotting the clusters 
Patient = 'PY04N007';
numSZ = 3;
szindices = [61 74 81];

clusterdir = '/home/adamli/MATLAB/code_adam/burns_adapted/cluster';
addpath(clusterdir);
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/EVC');

addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted');

for i=1:numSZ
    cluster = load(strcat('clusterAll_', num2str(szindices(i))));
    cluster = cluster.clusterAll;
    
    dur = size(cluster, 1);          % store length of the seizure region
    interval = linspace(1, dur, 500);   % create an interval of 500 points
    %% Scale to 500 seconds
    if dur ~= 500    %%%% Scale up and interpolate
       newcluster = interp1(1:dur, cluster, interval, 'linear');
    end
end