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

%%%%% Looking at Info mat And Getting the Seizure Times
infoevent = load('infoevent.mat');  % stores information about the seizure events
infotime = load('infotime.mat');    % carries information about the entire recording session

numpatients = length(fieldnames(infotime)); % should match w/ infoevent
eventname = strcat('event', Patient);
patseizevent = infoevent.(eventname);
patinfo = infotime.(Patient);
startrec = patinfo.time(1,1);

offset = str2num(timeWindow);

[seizfiles, TimesSZ, rectimes] = isolate_seizures(patseizevent, patinfo, offset);
iitimes = isolate_interictal(patseizevent, patinfo, offset);

%% 03: Save Which Seizures were Actual Seizures, and # seizures
if sum(Patient == 'PY04N012') == 8
  TimesSZ = TimesSZ([1 3 5 6],:);
  NumSZ = 4;
end

if sum(Patient == 'PY04N013') == 8
  TimesSZ = TimesSZ([1 12 71 72 73 81],:);
  NumSZ = 6;
end

if sum(Patient == 'PY04N015') == 8
  TimesSZ = TimesSZ([1 2 4 6],:);
  NumSZ = 4;
end

if sum(Patient == 'PY11N014') == 8
  TimesSZ = TimesSZ([3 5 6 7 11 12 14 15],:);
  NumSZ = 8;
end

if sum(Patient == 'PY05N004') == 8
  TimesSZ = TimesSZ([3 4 5 6],:);
  NumSZ = 4;
end

if sum(Patient == 'PY04N007') == 8
  TimesSZ = TimesSZ([4 5 6],:);
  seizfiles = seizfiles([4 5 6]);
  rectimes = rectimes([4 5 6],:);
  iitimes = iitimes([4 5 6 6+1],:);
  NumSZ = 3;
end

if sum(Patient == 'PY05N007') == 8
  TimesSZ = TimesSZ([1 2 3],:);
  NumSZ = 3;
end

%%%%%%%%%%%%% Rescaling the times %%%%%%%%%%%%%%%%%%%%%%%
TimesSZ = TimesSZ - TimesSZ(1,1)+1;  % vs. the beginning of the recording session
iitimes = iitimes - iitimes(1,1)+1;

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
    eval(['load ', mat, 'q' num2str(szindices(i)) 'pre_', Patient, '_', timeWindow])
    eval(['predata = q' num2str(szindices(i)) 'pre;'])

    eval(['load ', mat, 'q' num2str(szindices(i)) 'sz_', Patient, '_', timeWindow])
    eval(['szdata = q' num2str(szindices(i)) 'sz;'])

    eval(['load ', mat, 'q' num2str(szindices(i)) 'post_', Patient, '_', timeWindow])
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
    
    % store duration of the region
%     [X, Y] = meshgrid(1:500, 1:75);
    clusterdursz = size(clustersz, 1);
    
    % create evenly spacedintervals 
    totallen = 500; % total length of the area of analysis (e.g. 500 seconds)
    intervalsz = linspace(1, clusterdursz, totallen);   % create an interval of 500 points
    
    %%% Scale pre, seize, and post
    if clusterdursz ~= totallen-offset %%%% Scale up and interpolate
%        szdata = interp2(X, Y, szdata);
       clustersz = transpose(round(interp1(1:clusterdursz, clustersz, intervalsz, 'linear')));
    end
    
    %%% The vars in other plot
    rescaledCluster = [clusterpre; clustersz; clusterpost];
%     rescaleEVC{i} = [predata;szdata;postdata];
    EVC{i} = data;
    cstring{i} = strcat('seizure: ', num2str(szindices(i)));
    
    %%%%% Plot the rescaled clusters in separate plot
    A = figure(length(szindices)+1);
    B = subplot(4, 1, 4);
    hold on
    rescaledX = length(rescaledCluster);    
    plot(1:rescaledX, rescaledCluster);
    title(strcat('Clustering of States Pre, Ictal and Post for Patient: ', Patient));
    xlabel('time (sec)');
    ylabel('state number');
    xlim([0 length(cluster)])
    colorbar;
    
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
    line([str2num(timeWindow) str2num(timeWindow)], get(b, 'YLim'), 'LineStyle', '--', 'Color', 'black', 'LineWidth', linewidth);
    line([length(data)-str2num(timeWindow) length(data)-str2num(timeWindow)], get(b, 'YLim'), 'LineStyle', '--', 'Color', 'black', 'LineWidth', linewidth);
    
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
    line([str2num(timeWindow) str2num(timeWindow)], get(c, 'YLim'), 'LineStyle', '--', 'Color', 'black', 'LineWidth', linewidth);
    line([length(data)-str2num(timeWindow) length(data)-str2num(timeWindow)], get(c, 'YLim'), 'LineStyle', '--', 'Color', 'black', 'LineWidth', linewidth);
    
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
    line([str2num(timeWindow) str2num(timeWindow)], get(d, 'YLim'), 'LineStyle', '--', 'Color', 'black', 'LineWidth', linewidth);
    line([length(data)-str2num(timeWindow) length(data)-str2num(timeWindow)], get(d, 'YLim'), 'LineStyle', '--', 'Color', 'black', 'LineWidth', linewidth);
    
    cd(strcat('figures/', FreqBand))
    saveas(a, [num2str(szindices(i)) 'clusterEVC_' FreqBand timeWindow '.fig'])
%     movefile('*.fig','figures');
end

% also make separate plot of all seizures concatenated and all EVC charts
% TIME RESCALE
A = figure(length(szindices)+1);
legend(cstring);

for i=1:length(szindices)
    B = subplot(4, 1, i);
    hold on
    imagesc(transpose(EVC{i}));  
    xlim([0 length(newcluster)])
    B.YDir = 'normal';
    B.YTick = [1 20 40 60 75];
    set(gca, 'CLim', clim);
    colormap('jet');
    colorbar;
    title(['Seizure ' num2str(szindices(i)) ' EVC Spectrogram With ' FreqBand ' and time window: ' timeWindow])
    ylabel('channels');
    xlabel('time (seconds)');
end
saveas(A, ['allclusterEVC_' FreqBand timeWindow '.fig'])
