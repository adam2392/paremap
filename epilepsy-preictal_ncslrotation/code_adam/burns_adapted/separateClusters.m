%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script separateClusters.m
%
% Description: separates clusters made from 'GapStat2.m' by seizure index
%
% Output:   no output returned.
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 09/21/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 01: Create Settings 
Patient = 'PY04N013';
FreqBand = 'high';
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
timeWindow = '600'; 
% szindices = [70 73 77 79];
% szindices = [61 74 81]; %PY04N007
if strcmp(Patient, 'PY05N004')
    szindices = {'124A', '124B', '125'}; %PY05N004
elseif strcmp(Patient, 'PY04N013')
    szindices = {'762' '767' '777' '780' '783' '801'}; %PY04N013
elseif strcmp(Patient, 'PY04N008')
    szindices = {'8', '20'}; %PY04N008 
end
    
doII = 0;
mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_'
else
    mat = ''
end

for f=1:length(FreqBands)
    FreqBand = FreqBands{f};
    
    %% 02: Add Directory Paths and load related Mat files
    addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')
    addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster'))%, FreqBand))
    % addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrode_analysis/', FreqBand))

    %%%%%% Load in clustered files from 'GapStat2.m'
    % eval(['clusterAllPost = load(''', mat, 'clusterAllPost_', Patient, '_', FreqBand, timeWindow,''')'])
    % eval(['clusterAllPre = load(''', mat, 'clusterAllPre_', Patient, '_', FreqBand, timeWindow,''')'])
    eval(['clusterAllSZ = load(''', mat, 'clusterAllSZ_', Patient, '_', FreqBand, '600'')'])

    % store centroids and IDX labels
    % clusterAllPostIDX = clusterAllPost.IDXAllPost;
    % clusterAllPreIDX = clusterAllPre.IDXAllPre;
    clusterAllSZIDX = clusterAllSZ.IDXAllSZ;
    % clusterAllPostCents = clusterAllPost.CentsAllPost;
    % clusterAllPreCents = clusterAllPre.CentsAllPre;
    clusterAllSZCents = clusterAllSZ.CentsAllSZ;

    clusterAllSZIDX = clusterAllSZIDX;% + max(clusterAllPreIDX);
    % clusterAllPostIDX = clusterAllPostIDX + max(clusterAllSZIDX);

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

    %% 04: Get Length of Seizures
    szlen = zeros(length(szindices), 1);
    for j=1:length(szindices)
        seizurenum = szindices{j};
        % load in q seizure
        qsz = load(strcat(mattype, seizurenum, 'sz_', Patient, '_', timeWindow));
        field = fields(qsz);
        qsz = qsz.(field{:});     % reset qsz

        [newszlen, dumby] = size(qsz);

        szlen(j) = newszlen;
    end

    disp(['The seizure lengths of patient: ' Patient ' were '])
    disp(szlen)

    % concatenate all centroid values from the GapStat.m calculation
    clusterAllCents = [clusterAllSZCents]; %[clusterAllPreCents; clusterAllSZCents]; % clusterAllPostCents];

    %% 05: Create Clustered regions For Each Seizure Event
    if str2num(timeWindow) >= 3600
        timeWindow = '3600';
    end

    % cd(strcat('cluster/', FreqBand))
    for i=1:length(szindices)
        seizurenum = szindices{i};
        disp(['On seizure: ' num2str(i)])

        try
            %%% taking the length of pre,sz, and post is +/- 10 minutes
            %%% ***If you want 2.5 minutes -> do +/150 indices
    %         clusterPreIDX = clusterAllPreIDX((i-1)*str2num(timeWindow)+1:i*str2num(timeWindow));    
            if i==1
                clusterSZIDX = clusterAllSZIDX(1:szlen(i));
            elseif i==length(szindices)
                clusterSZIDX = clusterAllSZIDX(totalszlen:end);
            else
    %             disp([num2str(szindices(i)) 'else statement'])
                clusterSZIDX = clusterAllSZIDX(totalszlen:sum(szlen(1:i)));
            end
            totalszlen = sum(szlen(1:i));
    %         clusterPostIDX = clusterAllPostIDX((i-1)*str2num(timeWindow)+1:i*str2num(timeWindow)); 
        catch error
            disp('Catching error on indicing through the clusters created by GapStat2.m!')
            throw(error)
        end

        % concatenate all the labels
        clusterAllIDX = [clusterSZIDX];%[clusterPreIDX; clusterSZIDX];% clusterPostIDX];

        save(['cluster' seizurenum '_' Patient '_' FreqBand timeWindow], 'clusterAllIDX', 'clusterAllCents');
    end
end