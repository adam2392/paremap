%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script plotCluster.m
%
% Description: To help plot the clusters for each epileptic region: ictal,
% interictal, preictal and postictal
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
% clc
% clear
%% For plotting the clusters 
Patient = 'PY04N007';

doII = 0;

clusterdir = '/home/adamli/MATLAB/code_adam/burns_adapted/cluster';
addpath(clusterdir);
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/EVC');

addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster');

clusterAllPost = load(strcat('clusterAllPost_', Patient));
clusterAllPre = load(strcat('clusterAllPre_', Patient));
clusterAllSZ = load(strcat('clusterAllSZ_', Patient));
clusterAllPost = clusterAllPost.IDXAllPost;
clusterAllPre = clusterAllPre.IDXAllPre;
clusterAllSZ = clusterAllSZ.IDXAllSZ;

%% set clusters to distinct groups *** INHERENT ASSUMPTION THAT ALL
%%%% CLUSTERS ARE MUTUALLY EXCLUSIVE.... 
clusterAllSZ = clusterAllSZ + max(unique(clusterAllPre));
clusterAllPost = clusterAllPost + max(unique(clusterAllSZ));

%% get the clusters by seizure section...
%%%%% Only getting the first seizure right now!!!!!
numSZ = 3;
szindices = [61 74 81];

%% Looking at Info mat And Getting the Seizure Times
infoevent = load('infoevent.mat');  % stores information about the seizure events
infotime = load('infotime.mat');    % carries information about the entire recording session

numpatients = length(fieldnames(infotime)); % should match w/ infoevent
eventname = strcat('event', Patient);
patseizevent = infoevent.(eventname);
patinfo = infotime.(Patient);
startrec = patinfo.time(1,1)

%%%% Change offset matching the time frames of pre and post that you want
%%%% to see (e.g. 10 minutes -> offset = 600 (seconds))
offset = 150;  

[seizfiles, TimesSZ, rectimes] = isolate_seizures(patseizevent, patinfo, offset);
iitimes = isolate_interictal(patseizevent, patinfo, offset);

%% Save Which Seizures were Actual Seizures, and # seizures
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

% pointers to which index we left off on
prepointer = 0;
szpointer = 0;
postpointer = 0;
totalprelen = 0;
totalszlen = 0;
totalpostlen = 0;

%%%%% Settings for length of each zone
% 3xnumSZ array [ pre pre pre ..;
%                 sz   sz  sz ..;
%                 post post post..];
indexarray = zeros(3, numSZ);
for j=1:numSZ
    % load in pre, ictal and post EVC arrays to set indices
    % each q has their corresponding indexed seizure
    qpost = load(strcat('q', num2str(szindices(j)), 'post_', Patient));
    field = fields(qpost);
    qpost = qpost.(field{:});   % reset qpost
    qpre = load(strcat('q', num2str(szindices(j)), 'pre_', Patient));
    field = fields(qpre);
    qpre = qpre.(field{:});   % reset qpre
    qsz = load(strcat('q', num2str(szindices(j)), 'sz_', Patient));
    field = fields(qsz);
    qsz = qsz.(field{:});     % reset qsz
    
    [newprelen, dumby] = size(qpre);
    [newszlen, dumby] = size(qsz);
    [newpostlen, dumby] = size(qpost);
    
    indexarray(1,j) = newprelen;
    indexarray(2,j) = newszlen;
    indexarray(3,j) = newpostlen; 
end

% clear out variables not used
clear dumby newprelen newpostlen newszlen
clear qpost qpre qsz

%%%% Get all seizures
%%%% Get offset period of pre and post
%%%%%%% To understand logic, walk through pre,sz, and post separately
for i=1:numSZ
    disp(['On seizure: ' num2str(i)])
    % store original length of q files to be able to index through clusterAll

    % pointers to help navigate the indicing
    if i==1 
        prepointer = indexarray(1,i);
        szpointer = 1;
        postpointer = 1;  
    else
        prepointer = sum(indexarray(1, 1:i ));
        szpointer = sum(indexarray(2, 1:(i-1) ));    % point to start of seizure
        postpointer = sum(indexarray(3, 1:(i-1)));
    end
    
    try
        %%% taking the length of pre,sz, and post is +/- 10 minutes
        %%% ***If you want 2.5 minutes -> do +/150 indices
        clusterPre = clusterAllPre(prepointer - offset : prepointer - 1);                % get from end-offset : end
        clusterSZ = clusterAllSZ(szpointer : sum(indexarray(2,1:i)) - 1);                       % 1:length(seizure)
        clusterPost = clusterAllPost(postpointer : postpointer + offset - 1);             % 1:offset of post
    catch error
        disp('Catching error probably on indexing in "plotCluster.m"!')
        throw(error)
    end
    
    %%%%%% ******** Problem with data being too short (1 preictal is only
    %%%%%% 13 seconds?)*********
    if indexarray(1,i) <= offset
        clusterPre = clusterAllPre(prepointer - indexarray(1,i) : prepointer - 1);
        TimesSZ(i,1) = TimesSZ(i,1) + offset - indexarray(1,i) - 1;
        disp('Went inside if');
    end
    
    test1 = length(clusterPre);
    test2 = length(clusterSZ);
    test3 = length(clusterPost);
    disp([test1 test2 test3])
    length(TimesSZ(i,1):TimesSZ(i,2))
    
    %% concatenate clustering & Plotting
    clusterAll = [clusterPre; clusterSZ; clusterPost];

    A = figure(i);
    %%%% Do some fancy plotting in case off by an index of 1
    try
        plot(TimesSZ(i,1):TimesSZ(i,2), clusterAll);
    catch error
        disp(error)
        try
            plot(TimesSZ(i,1):TimesSZ(i,2), clusterAll(1:end-1));
        catch seconderror
            disp(['length of time segment, ' num2str(length(TimesSZ(i,1):TimesSZ(i,2)))])
            disp(['length of cluster segment, ' num2str(length(clusterAll))])
            throw(seconderror)
        end
    end
%     h = axes;
    title(strcat('Clustering of States Pre, Ictal and Post for Patient: ', Patient));
    xlabel('time (sec)');
    ylabel('state number');
    ylim([0 10]);
    saveas(A, strcat('cluster', num2str(szindices(i)), 'All_', Patient,'.fig'))
    
    save(strcat('clusterAll_', num2str(szindices(i))), 'clusterAll');
    movefile(strcat('clusterAll_', num2str(szindices(i)), '.mat'), 'cluster');
end

movefile('*.fig','figures');

