%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script plotPreictalEVCpart2.m
%
% Description: part 2 analysis of preictal state firing
%
% Output:   no output returned.
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 09/22/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% 01: Create Settings 
Patient = 'PY04N013';
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
timeWindow = '7200'; 
% szindices = {'61' '74' '81'};
% szindices = {'124A', '124B', '125'}; %PY05N004
% szindices = {'762' '767' '777' '780' '783' '801'}; %PY04N013
% szindices = {'8', '20'}; %PY04N008 
szindices = {'767' '777' '780' '783' '801'};

%% Load in files
FreqBand = 'high';
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrode_analysis/7200');
addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrode_analysis/', FreqBand))
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/preictal_analysis (perseizure, focus electrodes)')
% loop through each seizure event
for i=5:length(szindices)
    seizurenum = szindices{i};
    load([Patient '_' seizurenum '_concatenateEVC_allpredistances']);
    
%     % Load in the EVC data from 7200 out
%     disp(['q' seizurenum 'pre_' Patient FreqBand])
%     eval(['load ' 'q' seizurenum 'pre_' Patient FreqBand ';'])
%     eval(['data = q' seizurenum 'pre_' Patient FreqBand ';'])  
% 
%     % load in EVC data of 600 seconds out
%     eval(['load q' seizurenum 'pre_' Patient '_600;'])
%     
%      %% Grab a period of 600 seconds from each time Window
%     % 7200 seconds out
%     qpre7200 = data(1:600,:);
% 
%     % 3600 seconds out
%     qpre3600 = data(4001:4600,:);
% 
%     % 1800 seconds out
%     qpre1800 = data(end-1800+1:end-1800+600,:);
% 
%     % original 600 seconds out
%     eval(['qpre600 = q' seizurenum 'pre;'])
%     
%     
%     qallpre = [qpre7200;qpre3600;qpre1800;qpre600];
    close all
    figure(6);
    imagesc(distancesall)
    colorbar;
    
    D = figure(6);
    imagesc(transpose(qallpre));
    colorbar;
    title([Patient ' seizure number: ' num2str(i) ' EVC Plot']);
    xlabel('time');
    ylabel('channels');
    
%     saveas(D, [Patient '_' num2str(i) 'seizure_concatenateEVCplot.jpg']);
    %% load in previously computed data
    % re-compute centroids
    centroids = centsall;

    changecents = mean(centroids([1 2 5],:));
    newcents = [centroids([3 4],:); changecents];

    % create new distance matrix
    distancesall = pdist(newcents, 'cosine');
    distancesall = squareform(distancesall);
    A = figure(1);
    imagesc(distancesall)
    colorbar;

    % load in EVCs and try with the 3 states
    % use KD-Tree to assign to nearest 3 centroids
    %%% Build a KD-tree to perform nearest neighbor search based on the
    %%% Unsupervised PreCentroids
    mdl = createns(newcents, 'NSMethod', 'kdtree', 'Distance', 'minkowski', 'BucketSize', 80);
    newallcluster = knnsearch(mdl, qallpre, 'K', 1);

    B = figure(2);
    plot(newallcluster)
    title('new cluster plot after merging centroids')

    % firing rate
    Z = figure(3);
    uniqueidx = unique(newallcluster);
    % Show firing rate for each preictal state
    for k=1:length(uniqueidx)
       showstate = uniqueidx(k);
       spiketrain = newallcluster == showstate;

       % compute the firing rate
       step = 5;
       for l=1:(length(spiketrain)-1)/step
           if l==1
               index = 1;
           else
               index = l*step;
           end
           showfiring(l) = sum(spiketrain(index:index+step))/step;
       end

       z = subplot(3, 1, k);
       plot(showfiring);
       title(['Firing rate for state: ' num2str(uniqueidx(k)) ' and seizure: ' seizurenum])
    end
    
    % save figures
    saveas(A, [Patient '_' num2str(i) 'seizure_concatenatedistancemat.jpg']);
    saveas(B, [Patient '_' num2str(i) 'seizure_concatenateclusterplot.jpg']);
    saveas(Z, [Patient '_' num2str(i) 'seizure_concatenatefiringrates.jpg']);
    saveas(D, [Patient '_' num2str(i) 'seizure_concatenateEVCplot.jpg']);
end
% firing rate, cluster plot