%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script plotEVC.m
%
% Description: To help plot the eigenvector centrality time series as seen
% in the epilepsy paper Figure 3 E/F. Plots 4 different graphs for each
% seizure. First version allows for 3 seizures, but change 'numSZ' and 'szindices'
% to allow for different seizures. 
%
% Input:    
%
% Output:   no output returned.
%
%
%
% NOTE: Only concatenated pre, seize, post have the correct time scale
% gathered from seizure event data
%
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 08/21/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
close all

%% For plotting the EVCs 
numSZ = 3;
dataname = 'AllII';
doII = 0;
szindices = [61 74 81];
Patient = 'PY04N007';

evcdir = '/home/adamli/MATLAB/code_adam/burns_adapted/EVC';
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/EVC');
addpath('/home/adamli/MATLAB/code_adam');

addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam');

%% Looking at Info mat And Getting the Seizure Times
infoevent = load('infoevent.mat');  % stores information about the seizure events
infotime = load('infotime.mat');    % carries information about the entire recording session

numpatients = length(fieldnames(infotime)); % should match w/ infoevent
eventname = strcat('event', Patient);
patseizevent = infoevent.(eventname);
patinfo = infotime.(Patient);
startrec = patinfo.time(1,1);

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

iitimes
%%%%%%%%%%%%% Rescaling the times %%%%%%%%%%%%%%%%%%%%%%%
TimesSZ = TimesSZ - TimesSZ(1,1)+1;  % vs. the beginning of the recording session
iitimes = iitimes - iitimes(1,1)+1;

%%%%%%%%%%%%%%%% Note: qpre, qpost are originally set to +/- 10 minutes
%%%%%%%%%%%%%%%% from seizure onset/ending

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

%% Plotting Seizure
% Make separate plot for each seizure
for i=1:numSZ
    % plot imagesc and the clustering
%     h = [];
%     h(1) = subplot(1, 3, 1);
%     h(2) = subplot(1, 3, 2);
%     h(3) = subplot(1, 3, 3);
    
    % load in pre, ictal and post EVC arrays
    qpost = load(strcat('q', num2str(szindices(i)), 'post_', Patient));
    field = fields(qpost);
    qpost = qpost.(field{:});   % reset qpost
    qpre = load(strcat('q', num2str(szindices(i)), 'pre_', Patient));
    field = fields(qpre);
    qpre = qpre.(field{:});   % reset qpre
    qsz = load(strcat('q', num2str(szindices(i)), 'sz_', Patient));
    field = fields(qsz);
    qsz = qsz.(field{:});     % reset qsz
    
    %%%% Set to +/- 2.5 minutes if necessary
    try
        qpre = qpre((end-150):end,:);
        qpost = qpost(1:150,:);
    catch error
        disp(error);
    end
    qcat = [qpre; qsz; qpost]; % concatenate the EVCs from different periods

    %%%%%% ******** Problem with data being too short (1 preictal is only
    %%%%%% 13 seconds?)*********
    if indexarray(1,i) <= offset
        TimesSZ(i,1) = TimesSZ(i,1) + offset - indexarray(1,i) - 1;
        disp('Went inside if');
        
        %%% the vertical lines
        SPre= TimesSZ(i,1) + indexarray(1,i);   % point goes here
        SPost=length(qsz)+indexarray(1,i) + TimesSZ(i,1);   % point goes here
    else
        SPre= TimesSZ(i,1) + offset;   % point goes here
        SPost=length(qsz) + offset + TimesSZ(i,1);   % point goes here
    end
%%% Plot all areas concatenated
    [x, y] = size(qcat);
    D = figure(i);
    h = axes;
    imagesc(TimesSZ(i,1):TimesSZ(i,2), 0, transpose(qcat));
    set(h, 'Ydir', 'normal')
    set(h, 'YTick', [1 20 40 60 75])
    colormap('jet');
    c = colorbar;
    ylabel('channels');
    xlabel('time (seconds)');
    xlim([TimesSZ(i,1) TimesSZ(i,2)])
    title(strcat('pre, seize and post EVC (index: ', num2str(szindices(i)),')'))
	clim = get(gca, 'CLim');
    %%% plot vertical lines at different zones
    line([SPre SPre], get(h, 'YLim'), 'Color', 'black', 'LineWidth', 4);
    line([SPost SPost], get(h, 'YLim'), 'Color', 'black', 'LineWidth', 3);
    
% %%% Plot 'preictal' areaest = test.name
%     A = figure(4*i);
%     h = axes;
%     imagesc(transpose(qpre));
%     set(h, 'Ydir', 'normal')
%     set(h, 'YTick', [1 20 40 60 75])
%     colormap('jet');
%     c = colorbar;
%     ylabel('channels');
%     xlabel('time (seconds)');
%     title(strcat('preictal EVC (index: ', num2str(szindices(i)),')'))
%     set(gca,'CLim',clim);
%     
% %%% Plot 'ictal' area
%     B = figure(4*i+1);
%     h = axes;
%     imagesc(transpose(qsz));
%     set(h, 'Ydir', 'normal')
%     set(h, 'YTick', [1 20 40 60 75])
%     colormap('jet');
%     c = colorbar;
%     ylabel('channels');
%     xlabel('time (seconds)');
%     title(strcat('ictal EVC (index: ', num2str(szindices(i)),')'))
%     set(gca,'CLim',clim);
%     
% %%% Plot 'postictal' area
%     C = figure(4*i+2);
%     h = axes;
%     imagesc(transpose(qpost));
%     set(h, 'Ydir', 'normal')
%     set(h, 'YTick', [1 20 40 60 75])
%     colormap('jet');
%     c = colorbar;
%     ylabel('channels');
%     xlabel('time (seconds)');
%     title(strcat('postictal EVC (index: ', num2str(szindices(i)),')'))
%     set(gca,'CLim',clim);
   
    
%     saveas(A, strcat('q', num2str(szindices(i)), 'pre_', Patient,'.bmp'))
%     saveas(B, strcat('q', num2str(szindices(i)), 'sz_', Patient,'.bmp'))
%     saveas(C, strcat('q', num2str(szindices(i)), 'post_', Patient,'.bmp'))
    saveas(D, strcat('qAll', num2str(szindices(i)), Patient, '.fig'))
end
%% Plotting Interictal
if doII
    dataname = 'AllII';              %array to cluster
    eval(['load q',dataname, '_', Patient])
    eval(['data = q',dataname,';'])  %set data = qAllII;
    % Plot all interictal graphs separately
    for j=1:numSZ+1
        datatemp = data(iitimes(j,1):iitimes(j,2), :);

        b = figure(15+j);
        h = axes;
        imagesc(iitimes(j,1):iitimes(j,2), 0, transpose(data));
        set(h, 'Ydir', 'normal')
        set(h, 'YTick', [1 20 40 60 75])
        colormap('jet');
        c = colorbar;
        ylabel('channels');
        xlabel('time (seconds)');
        xlim([iitimes(j,1) iitimes(j,2)])
        title(['Interictal EVC W/O induced Seizures On zone: ' num2str(j)])
        set(gca, 'CLim', clim); 
        saveas(b, strcat('iiZone', num2str(j), '_', Patient,'.fig'))
    end

    %%% Plot all areas concatenated
    [x, y] = size(data);
    a = figure(15);
    h = axes;
    imagesc(iitimes(1,1):iitimes(end,2), 0, transpose(data));
    set(h, 'Ydir', 'normal')
    set(h, 'YTick', [1 20 40 60 75])
    colormap('jet');
    c = colorbar;
    ylabel('channels');
    xlabel('time (seconds)');
    xlim([iitimes(1,1) iitimes(end,2)])
    title('Interictal Concatenated EVC W/O induced Seizures')
    set(gca, 'CLim', clim);

    %%% plot vertical lines at different zones
    for i=1:numSZ
        segment = iitimes(i,2);
        line([segment segment], get(h, 'YLim'), 'Color', 'black', 'LineWidth', 4);
    end
    saveas(a, strcat('iiAll','_', Patient,'.fig'))
end
movefile('*.fig','figures');