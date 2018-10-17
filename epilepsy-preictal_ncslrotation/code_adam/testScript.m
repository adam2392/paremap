%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Written By: Adam Li
% Contact: ali39@jhu.edu
% 
% Written to analyze the SVD vectors 1st eigenvector over time (in/around
% the seizure zone, or in the interictal periods) for the different
% frequency bands.
% 
% Calculating the index which to analyze in the array of matrix
% eigenvectors by splitting up the recording time over the length of that
% array... Depends on uniform samling rate throughout trial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%% Add file path dependencies:
% read_svd_matrix(filename,nch), read_svd_array(filename,nch), read_adj_matrix(filename,nch)
% The .mat information files (intoevent.mat, infotime.mat)
addpath('/home/adamli/MATLAB/code_santaniello/utility software/read_data/');
addpath('/home/adamli/MATLAB/format_mat/');

%% Read in Eigen Vectors 
% Folder and Directory setup
% ***These 3 can change regularly depending on patient, # channels, hard drive
maindir = '/media/ExtHDD04/';   % main file directory with all the data
patientnum = 'PY04N007';
patientdir = strcat(patientnum,'/rec01/');  % patient dir we're interested in
numchan = 75%length([1:113 115:126]); % number of channels gotten from 'code_santaniello/utility software/svd_eeg/shell_svdeeg.m' file

svdvecdir = 'svd_vectors/';     % SVD vector dir
svdvaldir = 'svd_values/';      % SVD value dir
adjmatdir = 'adj_matrices/';    % Adjacency matrix dir

% used if you want to read adj matrix
testhelp = fullfile(maindir, patientdir, 'eeg.hdr');
test = readhdr(testhelp);
testnum = test.file_fmt.Numb_chans;

%% Looking at Info mat
infoevent = load('infoevent.mat');  % stores information about the seizure events
infotime = load('infotime.mat');    % carries information about the entire recording session

numpatients = length(fieldnames(infotime)); % should match w/ infoevent
eventname = strcat('event', patientnum);
patseizevent = infoevent.(eventname);
patinfo = infotime.(patientnum);

[seizfiles, seiztimes, rectimes] = isolate_seizures(patseizevent, patinfo, 150);

seizfiles = seizfiles(4:end);   %first 3 seizures don't count (induced)
seiztimes = seiztimes(4:end,:);
rectimes = rectimes(4:end, :);

cd(fullfile(maindir,patientdir));

%% Plotting/Producing Graphs
% build file names and loop through each frequency band file
%freqband = {'alpha'; 'beta'; 'gamma'; 'high'};  % the different freq. bands
for j=1:length(seizfiles)
    for i=1:1%length(freqband)
        freqband = {'beta'};
        
        fileindex = strfind(seizfiles{j}, '_');
        fileindex = seizfiles{j}(fileindex+1:fileindex+5);
        
        disp(strcat('looking at freq. band: ', freqband{i}));
        svdvecfile = strcat('svd_l_chr_data_', fileindex, '_',freqband{i},'.dat'); % SVD vector filename
        svdvalfile = strcat('svd_v_chr_data_', fileindex,'_',freqband{i},'.dat');
        adjmatfile = strcat('adj_chr_data_', fileindex,'_', freqband{i}, '.dat');

        % the full file path to the data files
        evecfile = fullfile(maindir, patientdir, svdvecdir, svdvecfile);
        evalfile = fullfile(maindir, patientdir, svdvaldir, svdvalfile);
        adjfile = fullfile(maindir, patientdir, adjmatdir, adjmatfile);

        % Generate matrix of eigenvectors, eigenvalues and adjacency matrices
        A = read_svd_matrix(evecfile, numchan); % list of matrices of eigenvectors (mxnxn)
        B = read_svd_array(evalfile, numchan);  % array of eigenvalues (nxm)
        % C = read_adj_matrix(adjfile, testnum); % only available if their is folder with files

        % Now calculate index in eigenvector matrix, we want to analyze
        % E.g. analyze +/- 150 seconds from seiztimes. We need to split up rectimes
        % into length(A) entries, and take the entries from seiztimes(:,1) to
        % seiztimes(:, 2)
        step = (rectimes(1,2) - rectimes(1,1))/length(A); %**Currently looking at the first seizfile...**
        startindex = abs(rectimes(1,1) - seiztimes(1,1))/step;
        endindex = abs(rectimes(1,1) - seiztimes(1,2))/step;

        startindex = round(startindex); % round
        endindex = round(endindex);

        %% Clustering 
        %%% Done on the first eigenvector (e.g. 'A(:,:,1)') over time
        %%% A(time, matrix, eigenvector) -> A(1, :, 1) is looking at first
        %%% time point, over all matrices of the first eigenvector
%         areatoanalyze = A(startindex:endindex,:,1);
% 
%         %%%%%% ? look at the period of seizure +/- 150 seconds ?
%         % K-means with ~10-14 states?
%         idx = kmeans(areatoanalyze, 14);
% 
%         % plot imagesc and the clustering
%         figure(i);
%         h = [];
%         h(1) = subplot(1, 2, 1);
%         h(2) = subplot(1, 2, 2);
%         imagesc(transpose(areatoanalyze), 'Parent', h(1));
%         colormap(jet);
%         c = colorbar;
%         plot(idx, 'Parent', h(2));
%         ylabel('Clustered states');
%         xlabel('time in ?');
    end
end