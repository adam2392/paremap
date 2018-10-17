%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script separateClusters02.m
%
% Description: separates clusters made from 'GapStat2.m' that were
% time-rescaled
%
% Output:   no output returned.
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 09/21/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 01: Create Settings 
Patient = 'PY04N007';
FreqBand = 'high';
timeWindow = '7200'; 
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

for i=1:length(szindices)
    %%%%%% Load in pre-clustered files from 'GapStat2.m'
    % eval(['clusterAllPost = load(''', mat, 'clusterAllPost_', Patient, '_', FreqBand, timeWindow,''')'])
    eval([num2str(szindices(i)) 'rescaled_clusterAllPre = load(''', mat, 'clusterAllPre_', Patient, '_', FreqBand, timeWindow,''')'])
    eval(['clusterAllSZ = load(''', mat, 'clusterAllSZ_', Patient, '_', FreqBand, '600'')'])

    % store centroids and IDX labels
    % clusterAllPostIDX = clusterAllPost.IDXAllPost;
    clusterAllPreIDX = clusterAllPre.IDXAllPre;
    clusterAllSZIDX = clusterAllSZ.IDXAllSZ;
    % clusterAllPostCents = clusterAllPost.CentsAllPost;
    clusterAllPreCents = clusterAllPre.CentsAllPre;
    clusterAllSZCents = clusterAllSZ.CentsAllSZ;

    clusterAllSZIDX = clusterAllSZIDX + max(clusterAllPreIDX);
    % clusterAllPostIDX = clusterAllPostIDX + max(clusterAllSZIDX);

    %% 04: Get Length of Seizures
    szlen = zeros(length(szindices), 1);
    for j=1:length(szindices)
        % load in q seizure
        qsz = load(strcat(mattype, num2str(szindices(j)), 'sz_', Patient));
        field = fields(qsz);
        qsz = qsz.(field{:});     % reset qsz

        [newszlen, dumby] = size(qsz);

        szlen(j) = newszlen;
    end

    disp(['The seizure lengths of patient: ' Patient ' were '])
    disp(szlen)

    % concatenate all centroid values from the GapStat.m calculation
%     clusterAllCents = [clusterAllPreCents; clusterAllSZCents; clusterAllPostCents];

    %% 05: Create Clustered regions For Each Seizure Event
    if str2num(timeWindow) >= 3600
        timeWindow = '3600';
    end

    cd(strcat('cluster/', FreqBand))
    for i=1:length(szindices)
        disp(['On seizure: ' num2str(i)])
    
        try
            %%% taking the length of pre,sz, and post is +/- 10 minutes
            %%% ***If you want 2.5 minutes -> do +/150 indices
            clusterPreIDX = clusterAllPreIDX((i-1)*str2num(timeWindow)+1:i*str2num(timeWindow));    
            if i==1
                clusterSZIDX = clusterAllSZIDX(1:szlen(i));
            elseif i==length(szindices)
                clusterSZIDX = clusterAllSZIDX(totalszlen:end);
            else
    %             disp([num2str(szindices(i)) 'else statement'])
                clusterSZIDX = clusterAllSZIDX(totalszlen:sum(szlen(1:i)));
            end
            totalszlen = sum(szlen(1:i));
            clusterPostIDX = clusterAllPostIDX((i-1)*str2num(timeWindow)+1:i*str2num(timeWindow)); 
        catch error
            disp('Catching error on indicing through the clusters created by GapStat2.m!')
            throw(error)
        end

        % concatenate all the labels
        clusterAllIDX = [clusterPreIDX; clusterSZIDX; clusterPostIDX];

        save(['cluster' num2str(szindices(i)) '_' Patient '_' FreqBand timeWindow], 'clusterAllIDX', 'clusterAllCents');
    end
end