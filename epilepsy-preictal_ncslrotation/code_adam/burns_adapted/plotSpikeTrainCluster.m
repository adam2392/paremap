%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script plotSpikeTrainCluster.m
%
% Description: plots the spike train using the original centroids from the
% GapStatistic, and then uses KD-Tree NN classification on preictal
% datasets of up to 240 mins.
%
% Then plots a spike train of how many states go to the "farthest" preictal
% state.
% 
% For each frequency, each seizure for a patient
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
Patient = 'PY04N013';
FreqBand = 'high';
% szindices = {'124A', '124B', '125'}; %PY05N004
szindices = {'762' '767' '777' '780' '783' '801'}; %PY04N013
% szindices = {'8', '20'}; %PY04N008 
% szindices = {'70', '73', '77', '79'}; %PY04N012
% szindices = {'61', '74', '81'}; %PY04N007
preclusters = 6;
distanceFunc = 'cosine'; %or cityblock, correlation, sqeuclidean
options = statset('MaxIter', 3000);


mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_'
else
    mat = ''
end

%% 02: Add Directory Paths and load related Mat files
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')
addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster/', FreqBand))
addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/EVC/', FreqBand))
% addpath('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrodes_analysis')

addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC/', FreqBand))
addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster/', FreqBand))
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/MatFiles')

%% Optional:
for i=1:length(szindices)
    seizurenum = szindices{i};
    %% 04: Determine which state corresponds to the most different state
    % load in the preictal EVC's of 7200
    % qpre
    newwindow = '3600';
    eval(['load ', mat, 'q' seizurenum 'pre_', Patient, '_', newwindow ';'])
    eval(['predata = q' seizurenum 'pre;']) % the data

    % define # of clusters and perform k-means
    preclusters = 5;
    [preidx, precents] = kmeans(predata, preclusters, 'Distance', distanceFunc, 'Options', options);  
    predistances = pdist(precents, 'cosine');
    predistances = squareform(predistances);

    %%% Get the most different state from the rest
    [ii,jj] = sort(predistances, 2, 'descend');
    v = ii(:,1);    % get the most different state's distance
    idx = jj(:,1);  % get the row/column it is in (state)
    iter = 0;
    % loop k-means alg. to make sure there is 1 "very" different state
    while length(unique(idx)) > 2
        iter = iter+1;
        [preidx, precents] = kmeans(predata, 6, 'Distance', distanceFunc, 'Options', options);  
        predistances = pdist(precents, 'cosine');
        predistances = squareform(predistances);

        [ii,jj] = sort(predistances,2,'descend');
        v = ii(:,1);
        idx = jj(:,1);
    end

    % plot the imagesc of the distance matrix
    figure(10*i);
    imagesc(transpose(predistances));
    colorbar;

    save(strcat(Patient, '_predistances_', seizurenum), 'predistances', 'precents')

        %%%% Optional: show firing rate for each of the preictal states
%     Z = figure(2*length(szindices)+i);
%     testidx = unique(preidx);
%     % Show firing rate for each preictal state
%     for k=1:length(testidx)
%        showstate = testidx(k);
%        spiketrain = newprecluster == showstate;
%        step = 5;
%        for l=1:length(spiketrain)/step
%            showfiring(l) = sum(spiketrain(l:l+step))/step;
%        end
%        z = subplot(3, 2, k);
%        plot(showfiring);
%        title(['Firing rate for state: ' num2str(testidx(k)) ' and seizure: ' seizurenum])
%     end
%     saveas(Z, ['testFiringRate_' seizurenum '.bmp'])
    
    %% Show spike train and firing rate
     hypostate = mode(idx); % store the state that is most different

    % compute spike train
    spiketrain = preidx == hypostate;
    
    % compute firing rate
    step = 5;
    for j=1:(length(spiketrain)-1)/step
       if j==1
           index = 1;
       else
           index = j*step;
       end
       firingrate(j) = sum(spiketrain(index:index+step))/step; 
    end
    
    % plotting unsupervised clusters, spike train, firing rate
    A = figure(i);
    a = subplot(3, 1, 1);
    plot(preidx)
    title(['the original unsupervised clusters of preictal state seizure: ', seizurenum])
    
    c = subplot(3, 1, 2);
    plot(spiketrain)
    title(['The spike train (How often state ' num2str(hypostate) ' occurs) '])
    
    d = subplot(3, 1, 3);
    plot(firingrate)
    title(['the firing rate over window of 30 seconds '])
    
    % saving the figure
    saveas(A, [seizurenum 'spikeTrainCluster_' newwindow FreqBand '.bmp'])
end


%% set window and load in EVC data file
window='600';
eval(['load ', mat, 'qAllPre_', Patient, '_', window])
predata = qAllPre;

% define # of clusters and perform k-means
preclusters = 5;
[preidx, precents] = kmeans(predata, preclusters, 'Distance', distanceFunc, 'Options', options);  
predistances = pdist(precents, 'cosine');
predistances = squareform(predistances);

%%% Get the most different state from the rest
[ii,jj] = sort(predistances,2,'descend');
v = ii(:,1);    % get the most different state's distance
idx = jj(:,1);  % get the row/column it is in (state)
iter = 0;
% loop k-means alg. to make sure there is 1 "very" different state
while length(unique(idx)) > 2
    iter = iter+1;
    [preidx, precents] = kmeans(predata, 6, 'Distance', distanceFunc, 'Options', options);  
    predistances = pdist(precents, 'cosine');
    predistances = squareform(predistances);

    [ii,jj] = sort(predistances,2,'descend');
    v = ii(:,1);
    idx = jj(:,1);
end

% plot the imagesc of the distance matrix
figure(10);
imagesc(transpose(predistances));
colorbar;

disp(['the unique idx is: '])
disp(unique(idx))

save(strcat(Patient, '_predistances'), 'predistances', 'precents')

%% 03: Load in EVC data to do preictal analysis
for i=1:length(szindices)    
    disp('Loaded preictal EVC')
    seizurenum = szindices{i};
    %% 04: Determine which state corresponds to the most different state
    % load in the preictal EVC's of 7200
    % qpre
    newwindow = '3600';
    eval(['load ', mat, 'q' seizurenum 'pre_', Patient, '_', newwindow ';'])
    eval(['predata = q' seizurenum 'pre;'])
%     eval(['load ', mat, 'q' num2str(szindices(i)) 'pre_', Patient, '_', newwindow ';'])
%     eval(['predata = q' num2str(szindices(i)) 'pre;'])    
    
    %%% Build a KD-tree to perform nearest neighbor search based on the
    %%% Unsupervised PreCentroids
    mdl = createns(precents, 'NSMethod', 'kdtree', 'Distance', 'minkowski', 'BucketSize', 80);
    newprecluster = knnsearch(mdl, predata, 'K', 1);
    unique(newprecluster);
    
    %%%% Optional: show firing rate for each of the preictal states
    Z = figure(2*length(szindices)+i);
    testidx = unique(preidx);
    % Show firing rate for each preictal state
    for k=1:length(testidx)
       showstate = testidx(k);
       spiketrain = newprecluster == showstate;
       step = 5;
       for l=1:length(spiketrain)/step
           showfiring(l) = sum(spiketrain(l:l+step))/step;
       end
       z = subplot(3, 2, k);
       plot(showfiring);
       title(['Firing rate for state: ' num2str(testidx(k)) ' and seizure: ' seizurenum])
    end
    saveas(Z, ['testFiringRate_' seizurenum '.bmp'])

    %% 05: Create the spike train based on that
    hypostate = mode(idx); % store the state that is most different
%     disp(['state most different is: ' num2str(hypostate)])
    spiketrain = newprecluster == hypostate;
    step = 5;
    for j=1:(length(spiketrain)-1)/step
       if j==1
           index = 1;
       else
           index = j*step;
       end
       firingrate(j) = sum(spiketrain(index:index+step))/step; 
    end
    
    A = figure(i);
    a = subplot(4, 1, 1);
    plot(preidx)
    title(['the original unsupervised clusters of preictal state seizure: ', seizurenum])
    
    b = subplot(4, 1, 2);
    plot(newprecluster)
    title(['the new clusters of preictal state ' newwindow ' seconds out'])
    
    c = subplot(4, 1, 3);
    plot(spiketrain)
    title(['The spike train (How often state ' num2str(hypostate) ' occurs) '])
    
    d = subplot(4, 1, 4);
    plot(firingrate)
    title(['the firing rate over window of 30 seconds '])
    
    saveas(A, [seizurenum 'spikeTrainCluster_' newwindow FreqBand '.bmp'])

    clear firingrate spiketrain hypostate newprecluster
end
