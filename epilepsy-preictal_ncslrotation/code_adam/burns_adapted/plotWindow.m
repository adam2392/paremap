%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script plotClusterEVC.m
%
% Description: plots the separated EVC & corresponding cluster labels for
% the different freq. bands and time periods. FOR THE RESCALED CLUSTERS
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
Patient = 'PY04N007';
FreqBand = 'beta';
timeWindow = '3600'; 
szindices = [61 74 81];

doII = 0;
mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_'
else
    mat = ''
end

%% 02: Add Directory Paths and load related Mat files
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')
addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster/', FreqBand))

% Make separate plot for each seizure event
for i=1:length(szindices)
    cd('/home/adamli/MATLAB/code_adam/burns_adapted')
    %%%%%% Load in clustered files from 'separateClusters.m'
    if str2num(timeWindow) >= 3600
        eval(['cluster = load(''', mat, 'cluster', num2str(szindices(i)), '_', Patient, '_', FreqBand, timeWindow,''')'])
    else
        eval(['cluster = load(''', mat, num2str(szindices(i)), 'rescaled_clusterAllPre_', Patient, '_' FreqBand, timeWindow, ''')'])
    end
    cents = cluster.clusterAllCents;
    cluster = cluster.clusterAllIDX;
    disp('Loaded cluster')
%     %%% algorithm to re-assign cluster states based on timeWindow
%     maxpre = max(cluster(1:str2num(timeWindow)));
%     cluster(str2num(timeWindow)+1:end-str2num(timeWindow)) = cluster(str2num(timeWindow)+1:end-str2num(timeWindow)) + maxpre;
%     cluster(end-str2num(timeWindow):end) = cluster(end-str2num(timeWindow):end) + max(cluster(str2num(timeWindow)+1:end-str2num(timeWindow)));
    
    %%%%% Load in EVC files from MakeQ2.m
    eval(['load ', mat, 'q' num2str(szindices(i)) 'pre_', Patient, '_', timeWindow, ';'])
    eval(['predata = q' num2str(szindices(i)) 'pre;'])

    eval(['load ', mat, 'q' num2str(szindices(i)) 'sz_', Patient, '_', timeWindow, ';'])
    eval(['szdata = q' num2str(szindices(i)) 'sz;'])

    eval(['load ', mat, 'q' num2str(szindices(i)) 'post_', Patient, '_', timeWindow, ';'])
    eval(['postdata = q' num2str(szindices(i)) 'post;'])
    
    disp('Loaded EVC')
    
    % concatenate the EVC for one seizure event
    data = [predata;szdata;postdata];
    
    %%% Build a KD-tree to perform nearest neighbor search based on the GapStatistic
    mdl = createns(cents, 'NSMethod', 'kdtree', 'Distance', 'minkowski', 'BucketSize', 80);
    newcluster = knnsearch(mdl, data, 'K', 1);
    
    %% ** Store data, cluster, 'newcluster & TimeRescale SZ
    clusterpre = cluster(1:str2num(timeWindow));
    clusterpost = cluster(end-str2num(timeWindow)+1:end);
    clustersz = cluster(str2num(timeWindow)+1:end-str2num(timeWindow));
    
    window = 200;
    clusterpre = clusterpre(end-window:end);
    data = [predata(end-window:end,:); szdata];
    newclustersz = newcluster(str2num(timeWindow)+1:end-str2num(timeWindow));
    newclusterpre = newcluster(1:str2num(timeWindow));
    newclusterpre = newclusterpre(end-window:end);
    
    length(cluster)
    length(newcluster)
    cluster = [clusterpre;clustersz];
    newcluster = [newclusterpre; newclustersz];
    
    qclusters = 6; % the # of clusters we want to enforce
    distanceFunc = 'cosine'; %or cityblock, correlation, sqeuclidean
    options = statset('MaxIter', 3000);
    % Cluster k-means
    [kcluster, kcentroids] = kmeans(data, qclusters, 'Distance', distanceFunc, 'Options', options); 
    k = figure(i+10)
    plot(kcluster)
    xlim([0 length(kcluster)])
    ylim([-1 7])
    line([window window], [-1 7], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);
    
    %% 04: Plotting and Saving Figures
    a = figure(i);
    linewidth = 2;
%     set(gcf, 'PaperPositionMode', 'auto')
    
    clim = [0 0.2]
    %%%% plot of EVC Spectrogram
    b = subplot(3, 1, 1);
    axes(b)
    imagesc(transpose(data));
    b.YDir = 'normal';
    b.YTick = [1 20 40 60 75];
    colormap('jet');
    colorbar;
    ylabel('channels');
    xlabel('time (seconds)');
%     xlim([iitimes(1,1) iitimes(end,2)])
    title(['Seizure ' num2str(szindices(i)) ' EVC Spectrogram With ' FreqBand ' and time window: ' timeWindow])
    set(gca, 'CLim', clim);
    
    %%%% plot of clustered states
    c = subplot(3, 1, 2);
    axes(c)
    xrange = 1:length(cluster);
    plot(xrange, cluster)
    xlim([0 length(cluster)])
    ylabel('states');
    xlabel('time (seconds)');
    title(['Seizure ' num2str(szindices(i)) ' cluster states With ' FreqBand ' and time window: ' timeWindow])
    colorbar;

    %%%% plot of new clustered states after KD-tree nearest neighbor
    %%%% classificataion
    d = subplot(3, 1, 3);
    axes(d)
    xrange = 1:length(newcluster);
    plot(xrange, newcluster)
    xlim([0 length(newcluster)])
    ylabel('states');
    xlabel('time (seconds)');
    title(['Seizure ' num2str(szindices(i)) ' New Clusters With ' FreqBand ' and time window: ' timeWindow])
    colorbar;
 
    
%     cd(strcat('figures/', FreqBand))
%     saveas(a, [num2str(szindices(i)) 'clusterEVC_' FreqBand timeWindow '.fig'])
%     movefile('*.fig','figures');
end


