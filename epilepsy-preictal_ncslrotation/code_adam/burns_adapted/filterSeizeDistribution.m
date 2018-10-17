%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script filterSeizeDistribution.m
%
% Description: Filter the seizure distribution during seizure dynamics of
% each of the seizure. Then create plots to allow for computing it
%
% Then plot histogram over time
% 
% Output:   no output returned.
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 09/22/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
%% 01: Create Settings 
Patient = 'PY05N004';
FreqBand = 'high';
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
% szindices = [70 73 77 79];
% szindices = [61 74 81]; %PY04N007
szindices = {'124A', '124B', '125'}; %PY05N004
% szindices = {'762' '767' '777' '780' '783' '801'}; %PY04N013
% szindices = {'8', '20'}; %PY04N008 

szclusters = 6;
distanceFunc = 'cosine'; %or cityblock, correlation, sqeuclidean
options = statset('MaxIter', 3000);

mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_'
else
    mat = ''
end

eval(['Subject = ' 'Patient;'])

%% 02: Add Directory Paths and load related Mat files
% addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')
% addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster/', FreqBand))
% addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/EVC/', FreqBand))

% addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC/', FreqBand))
% addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster/', FreqBand))
% addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/MatFiles')

disp('Manually code in the index of the seizure centroids!')
disp('Look at 2nd for loop "load in centroids and idx values of seizure"')

%% 03: Load in Clu7ster data to do ictal analysis
% loop through FreqBands
for j=1:length(FreqBands)
    FreqBand = FreqBands{j}
    addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster'))
%     addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster/', FreqBand))
    for i=1:length(szindices)
        seizurenum = szindices{i};
        window = '600';
        %%%%% Load in preictal cluster files from GapStat.m
        eval(['load ', 'cluster', seizurenum, '_', Patient, '_', FreqBand, window])    

        % load in the centriods and idx values of seizure
        szcents = clusterAllCents;%clusterAllCents(2:2+5,:); % get the 2nd -> 7th column in centroid matrix
        szidx = clusterAllIDX; %clusterAllIDX(601:end-600);

        % rescale szidx to uniform length
        clusterdursz = size(szidx, 1);

        % create evenly spacedintervals 
        totallen = 250; % total length of the area of analysis (e.g. 500 seconds)
        intervalsz = transpose(linspace(1, clusterdursz, totallen));   % create an interval of 500 points

        %%% scale the seizure 
        if clusterdursz ~= totallen %%%% Scale up and interpolate
           newszidx = (round(interp1(1:clusterdursz, szidx, intervalsz, 'linear')));
        end

        A = figure(j);
        a = subplot(length(szindices), 2, 2*i-1);
        plot(newszidx)
        title([FreqBand ' band | time scaled seizure states for seizure#: ' seizurenum])
        ylim([0 9]);
        ax = gca;
        ax.YTick = [0 1 2 3 4 5 6 7 8 9];
        xlabel('seconds')
        ylabel('seizure states')

        %% 03b: Filter when the states change
        if strcmp(Patient, 'PY04N013')
            if strcmp(FreqBand, 'alpha')
                medfilteredidx = medfilt1(newszidx, 50);
            else
                medfilteredidx = medfilt1(newszidx, 40);
            end
        else
            medfilteredidx = medfilt1(newszidx, 40);
        end
            
        figure(j)
%         B = figure(length(FreqBands)+j);
        b = subplot(length(szindices),2,2*i);
        fig = gcf;
        plot(medfilteredidx)
        title([FreqBand ' band | median filtered states for seizure#: ' seizurenum])
        ylim([0 9]);
        ax = gca;
        ax.YTick = [0 1 2 3 4 5 6 7 8 9];
        xlabel('seconds')
        ylabel('seizure states')
        
        % store median filtered seizure distrib. in a struct
        field = strcat('seizure_', seizurenum);
        freq.(field).filteredidx = medfilteredidx;
        freq.(field).originalidx = newszidx;
        freq.(field).szcents = szcents;
    end
    % save figure as .fig
    saveas(A, strcat(Patient, '_', FreqBand, '.fig'))
    
    % store this corresponding frequency band's struct
    eval([FreqBand ' = freq;']);
end

filename = strcat(Patient, '_seizureDistribution');
save(filename, FreqBands{:});

%% Optional
% for i=1:length(szindices) 
%     % load in the preictal centroids from the +/- 10 min 
%     window = '600';
%     %%%%% Load in preictal EVC files from MakeQ2.m
%     eval(['load ', mat, 'q' num2str(szindices(i)) 'sz_', Patient, '_', window])
%     eval(['szdata = q' num2str(szindices(i)) 'sz;'])
%     
%     disp('Loaded ictal EVC')
%     
%     szclusters = 5; % 5 seizure states
%     [szidx, szcents] = kmeans(szdata, szclusters, 'Distance', distanceFunc, 'Options', options);  
%     szdistances = pdist(szcents, 'cosine');
%     szdistances = squareform(szdistances);
%     
%     %%% Get the most different state from the rest
%     [ii,jj] = sort(szdistances,2,'descend');
%     v = ii(:,1);    % get the most different state's distance
%     idx = jj(:,1);  % get the row/column it is in (state)
%     iter = 0;
%     % loop k-means alg. to make sure there is 1 "very" different state
%     while length(unique(idx)) > 2
%         iter = iter+1;
%         [preidx, szcents] = kmeans(predata, 5, 'Distance', distanceFunc, 'Options', options);  
%         szdistances = pdist(szcents, 'cosine');
%         szdistances = squareform(szdistances);
% 
%         [ii,jj] = sort(szdistances,2,'descend');
%         v = ii(:,1);
%         idx = jj(:,1);
%     end
%     unique(idx)
% %     
%     figure(length(szindices)+i);
%     imagesc(szdistances)
%     
%     %% 04: Determine which state corresponds to the most different state
%     % load in the preictal EVC's of 7200
%     % qpre
%     eval(['load ', mat, 'q' num2str(szindices(i)) 'pre_', Patient, '_', '7200;'])
%     eval(['predata = q' num2str(szindices(i)) 'pre;'])
%     
% %     predata = predata(end-3600:end,:);
%     
%     %%% Build a KD-tree to perform nearest neighbor search based on the
%     %%% Unsupervised PreCentroids
%     mdl = createns(szcents, 'NSMethod', 'kdtree', 'Distance', 'minkowski', 'BucketSize', 80);
%     newprecluster = knnsearch(mdl, predata, 'K', 1);
% 
%     % Create the spike train based on that
%     hypostate = mode(idx); % store the state that is most different
%     spiketrain = newprecluster == hypostate;
%     step = 7200/200;
%     for j=1:length(spiketrain)/step
%        firingrate(j) = sum(spiketrain(j:j+step))/step; 
%     end
%     
%     A = figure(i);
%     a = subplot(4, 1, 1);
%     plot(preidx)
%     title(['the original unsupervised clusters of preictal state seizure: ', num2str(szindices(i))])
%     
%     b = subplot(4, 1, 2);
%     plot(newprecluster)
%     title(['the new clusters of preictal state 7200 seconds out'])
%     
%     c = subplot(4, 1, 3);
%     plot(spiketrain)
%     title(['The spike train (How often "most different" state occurs) '])
%     
%     d = subplot(4, 1, 4);
%     plot(firingrate)
%     title(['the firing rate over window of 30 seconds '])
%     
%     saveas(A, [num2str(szindices(i)) 'spikeTrainCluster_' FreqBand '.bmp'])
% 
%     clear firingrate spiketrain hypostate newprecluster preidx
% end
