%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script unsupervisedCluster.m
%
% Description: performs unsupervised clustering (k-means) of the different
% EVC vectors of the different freq. bands and different time periods.
% 
% ** I should add in sampling of rest of the interictal and performing
% clustering and checking how the samples are close in euclidean distance.
%
% Output:   no output returned.
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 09/22/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
% pool = parpool;                      % Invokes workers
% stream = RandStream('mlfg6331_64');  % Random number stream
% options = statset('UseParallel',1,'UseSubstreams',1,...
%     'Streams',stream);

tic;
%% 01: Create Settings 
Patient = 'PY04N007';
FreqBand = 'high';
timeWindow = '7200'; 
szindices = [61 74 81];
qclusters = 12; % the # of clusters we want to enforce
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
addpath(['/home/adamli/MATLAB/code_adam/burns_adapted/EVC/', FreqBand])

%% 03: Extract information/data from Mat files
%%% Looking at Info mat And Getting the Seizure Times
infoevent = load('infoevent.mat');  % stores information about the seizure events
infotime = load('infotime.mat');    % carries information about the entire recording session

numpatients = length(fieldnames(infotime)); % should match w/ infoevent
eventname = strcat('event', Patient);
patseizevent = infoevent.(eventname);
patinfo = infotime.(Patient);
startrec = patinfo.time(1,1) 

offset = str2num(timeWindow);
[seizfiles, TimesSZ, rectimes] = isolate_seizures(patseizevent, patinfo, offset);
iitimes = isolate_interictal(patseizevent, patinfo, offset);

%% 03b: Save Which Seizures were Actual Seizures, and # seizures
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

%% 04: Load in EVC files from MakeQ2.m
for i=1:length(szindices)
    % qpre
    eval(['load ', mat, 'q' num2str(szindices(i)) 'pre_', Patient, '_', timeWindow])
    eval(['predata = q' num2str(szindices(i)) 'pre;'])

    % qsz
    eval(['load ', mat, 'q' num2str(szindices(i)) 'sz_', Patient, '_600'])
    eval(['szdata = q' num2str(szindices(i)) 'sz;'])

    % qpost
    eval(['load ', mat, 'q' num2str(szindices(i)) 'post_', Patient, '_', timeWindow])
    eval(['postdata = q' num2str(szindices(i)) 'post;'])
    
    % concatenate the EVC for one seizure event
    data = [predata;szdata];%;postdata];
    
    %% 05: K-Means Clustering & Distance Matrix
    [idx, cents] = kmeans(data, qclusters, 'Distance', distanceFunc, 'Options', options);    
%     preidx = idx(1:str2num(timeWindow));
    
    preclusters = 6;
    [preidx, precents] = kmeans(predata, 5, 'Distance', distanceFunc, 'Options', options);  
    predistances = pdist(precents, 'cosine');
    predistances = squareform(predistances);
    
    % compute a distance matrix
    distances = pdist(cents, 'cosine');
    distances = squareform(distances);
    
    %%%% Save k-means & distance matrix
    try
        cd(strcat('figures/', FreqBand))
    catch
        disp(['already in '  FreqBand ' folder'])
    end
    eval(['save unsuperCluster', num2str(szindices(i)), '_', Patient, '_', FreqBand, timeWindow, ' idx cents'])
    eval(['save distMat', num2str(szindices(i)), '_', Patient, '_', FreqBand, timeWindow, ' distances'])
    
    %% 06: Plot and Save Figure
    a = figure(i);
    
     %%%% ai) plot of unsupervised clusters with all periods
    b = subplot(2, 2, 1);
    axes(b);
    plot(preidx)
    ylabel('states');
    xlabel('time (seconds)');
    ylim([0 10])
    title([num2str(szindices(i)) ' Unsuper Preictal Clusters ' FreqBand ' and time window: ' timeWindow])

    %%%% aii) plot of distance matrix for a)
    c = subplot(2, 2, 3);
    imagesc(predistances)
    colorbar;
    set(c, 'YTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14])
    set(c, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14])
    ylabel('states');
    xlabel('states');
    title([num2str(szindices(i)) ' Distance Matrix of Pre-Clusters ' FreqBand ' and time window: ' timeWindow])
    
    %%%% bi) plot of unsupervised clusters with all periods
    d = subplot(2, 2, 2);
    axes(d);
    plot(idx)
    ylabel('states');
    xlabel('time (seconds)');
    ylim([0 15])
    title([num2str(szindices(i)) ' Unsuper Clusters ' FreqBand ' and time window: ' timeWindow])

    %%%% bii) plot of distance matrix for b)
    e = subplot(2, 2, 4);
%     h = axes(e);
    imagesc(distances)
    colorbar;
    set(e, 'YTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14])
    set(e, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14])
    ylabel('states');
    xlabel('states');
     
%     saveas(a, [num2str(szindices(i)) 'unsuperCluster_' FreqBand timeWindow '.fig'])
end
toc;
% delete(pool)
