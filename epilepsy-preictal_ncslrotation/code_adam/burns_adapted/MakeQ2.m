%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original File written by Samuel Burns
% Adapted by: Adam Li
%
% Ver.: 1.0 - Date: 08/14/2015
%
% Description: To generate array of eigenvectors(EVC) for time period. This
% is Network Analysis Part ii) and iii).
% 
% Dependencies: i) readhdr.m - read header file
%              ii) infoevent.mat - a mat file containing information about
%              seizure events during the recording sessions
%              iii) infotime.mat - a mat file containing information about
%              recording sessions
%
% Output: Concatenated arrays were produced:
% i.	Seizure EVCs called qAllSZ
% ii.	Preictal time periods called qAllPre
% iii.	Postictal called qAllPost
% iv.	Interictal called qAllII (in ClusterII2.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Pat* Info* Num* n2 *Length OverLap sfreq stimfreq
Patient = 'PY04N013';
debug_on = 1;

sfreq = 1000;
mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_'
else
    mat = ''
end

% set frequency bands for each patient
if     sum(Patient == 'PY04N007') == 8
  FreqBand = 'beta';
elseif sum(Patient == 'PY04N008') == 8
  FreqBand = 'beta';
elseif sum(Patient == 'PY04N012') == 8
  FreqBand = 'beta';
elseif sum(Patient == 'PY04N013') == 8
  FreqBand = 'beta';
elseif sum(Patient == 'PY04N015') == 8
  FreqBand = 'beta';
elseif sum(Patient == 'PY05N004') == 8
  FreqBand = 'alpha';
elseif sum(Patient == 'PY05N005') == 8
  FreqBand = 'beta';
elseif sum(Patient == 'PY05N006') == 8
  FreqBand = 'theta';
elseif sum(Patient == 'PY05N007') == 8
  FreqBand = 'beta';
end

if     sum(Patient == 'PY11N003') == 8
  FreqBand = 'gamma';
elseif sum(Patient == 'PY11N004') == 8
  FreqBand = 'gamma';
elseif sum(Patient == 'PY11N006') == 8
  FreqBand = 'beta';
elseif sum(Patient == 'PY11N009') == 8
  FreqBand = 'theta';
elseif sum(Patient == 'PY11N011') == 8
  FreqBand = 'theta';
  sfreq = 200; % why is this 200?
elseif sum(Patient == 'PY11N014') == 8
  FreqBand = 'beta';
end

% add utility software functions for helping read in data 
addpath('/home/adamli/MATLAB/code_santaniello/utility software/generate_pbs')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/read_data')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/svd_eeg')
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/')
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')
% addpath('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrode_analysis/')

% go to the directory with the mat files that hold information about
% seizure events, and recording sessions, and avg. Adj. Matrix and the std.
% Adj. matrix
homedir = '/home/adamli/MATLAB/code_adam/burns_adapted';
cd(homedir);
eval(['load(''infoevent.mat'',''event',Patient,''');'])
eval(['InfoEvent = event',Patient,';'])
eval(['load(''infotime.mat'',''',Patient,''');'])
eval(['InfoTime = ',Patient,';'])

%%% load in average, and standard deviations of each (i,j) entry in adj.
%%% matrix
% eval(['load mAllAdjMat_',Patient,'_',FreqBand])
% eval(['load mAllAdjSTD_',Patient,'_',FreqBand]) % original
%% Determine Patient Settings
clear All* NumWins NumChans* FPindex FreqBand
NumWins = 0;
sfreq = 1000;

% go through each PY11* patient and set:
% 1. number of channels
% 2. ???? FPindex
% 3. the main frequency band from R-spec peak freq. band
% 4. cd to that patient's directory
if     sum(Patient == 'PY04N007') == 8
  NumChans2 = 75;
  FPindex = 12:16;
  FreqBand = 'beta';
  eval(['cd /media/ExtHDD04/',Patient])
elseif sum(Patient == 'PY04N008') == 8
  NumChans2 = 40;
  FPindex = 12:16;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY04N012') == 8
  NumChans2 = 82;
  FPindex = 11:15;
  FreqBand = 'beta';
%   eval(['cd /media/ExtHDD04/',Patient])
eval(['cd /home/adamli/Desktop/',Patient])
elseif sum(Patient == 'PY04N013') == 8
  NumChans2 = 79;
  FPindex = 12:16;
  FreqBand = 'beta';
%   eval(['cd /media/ExtHDD02/',Patient])
elseif sum(Patient == 'PY04N015') == 8
  NumChans2 = 89;
  FPindex = 11:15;
  FreqBand = 'beta';
  eval(['cd /media/ExtHDD02/',Patient])
elseif sum(Patient == 'PY05N004') == 8
  NumChans2 = 94;
  FPindex = 11:15;
  FreqBand = 'alpha';
%   eval(['cd /media/ExtHDD01/',Patient])
elseif sum(Patient == 'PY05N005') == 8
  NumChans2 = 110;
  FPindex = 11:15;
  FreqBand = 'beta';
  eval(['cd /media/ExtHDD01/',Patient])
elseif sum(Patient == 'PY05N006') == 8
  NumChans2 = 22;
  FPindex = 11:15;
  FreqBand = 'theta';
  eval(['cd /media/ExtHDD01/',Patient])
elseif sum(Patient == 'PY05N007') == 8
  NumChans2 = 42;
  FPindex = 11:15;
  FreqBand = 'beta';
  eval(['cd /media/ExtHDD01/',Patient])
end

if     sum(Patient == 'PY11N003') == 8
  NumChans2 = 116;
  FPindex = 11:15;
  FreqBand = 'gamma';
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N004') == 8
  NumChans2 = 112;
  FPindex = 11:15;
  FreqBand = 'gamma';
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N006') == 8
  NumChans2 = 65;
  FPindex = 11:15;
  FreqBand = 'beta';
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N009') == 8
  NumChans2 = 99;
  FPindex = 11:15;
  FreqBand = 'theta';
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N011') == 8
  NumChans2 = 27;
  FPindex = 11:15;
  FreqBand = 'theta';
  sfreq = 200;
  eval(['cd /media/ExtHDD04/',Patient])
elseif sum(Patient == 'PY11N014') == 8
  NumChans2 = 104;
  FPindex = 11:15;
  FreqBand = 'beta';
  sfreq = 1000;
  eval(['cd /media/ExtHDD04/',Patient])
end

%%%%%%%%
% FreqBands = 'gamma'; % gamma band with high frequencies
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
% To Define timeWindow - 10 minutes = 600, 30 = 1800, 60 = 3600, 120 = 7200...
timeWindow = 7200;
%%%%%%%%

for f=1:length(FreqBands)
    FreqBand = FreqBands{f}
    
    % eval(['cd /media/',HDDrive,'/',Patient])
    eval(['cd /home/adamli/Desktop/', Patient])

    %% *** For looking at focus electrodes/channels only
    focuselectrodes = 0;
    if focuselectrodes
        % PY04N008:
        if strcmp(Patient, 'PY04N008')
            chans = [1, 2, 3, 9, 10, 11];
            NumChans2 = 6;
            focuselectrodes = 1;
        end

        % % PY04N013
        if strcmp(Patient, 'PY04N013')
            chans = [5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75];
            NumChans2 = 44;
            focuselectrodes = 1;
        end
        %   % PY05N004
        if strcmp(Patient, 'PY05N004')
            Patient
            chans = [80,81,82,83,84,85,88,89,90,91,92,93];
            NumChans2 = 12;
            focuselectrodes = 1;
        end
    end

    % load in mean/var adjacency matrices
    FreqBand
    eval(['load ', mat, 'mAllAdjMat_',Patient,'_',FreqBand])
    eval(['load ', mat, 'mAllAdjSTD_',Patient,'_',FreqBand])

    %% Power Spectrum Analysis
    WinLength = 3*sfreq;
    FFTlength = sfreq;
    Overlap = 1/2; 
    NumFFTs = WinLength/(Overlap*FFTlength) - 1;
    %NumFFTs = WinLength/FFTlength + (WinLength/FFTlength-1);
    stimfreq = sfreq*(0:FFTlength/2)/FFTlength;
    stimfreq = stimfreq(1:sfreq/2);

    PatTimes = InfoTime.time;       % load patient recording settings
    [NumFiles,n2] = size(PatTimes); % obtain size of recording time's array
    %eval(['cd /media/ExtHDD02/',Patient])
    %eval(['cd /home/sburns/Data/JHU/JHU2012/',Patient])

    clear NumSZ SZ* TimesSZ s2 TimesAll
    [NumSZ,s2] = size(InfoEvent.code);
    TimesAll = InfoTime.time;
    TimesSZ  = InfoEvent.time;

    if sum(Patient == 'PY11N014') == 8
      TimesAll(14:end,:) = TimesAll(14:end,:) + (24*3600)*19;
      TimesSZ(5:end,:)   = TimesSZ(5:end,:)   + (24*3600)*30;
    end

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
      NumSZ = 3;
    end

    if sum(Patient == 'PY05N007') == 8
      TimesSZ = TimesSZ([1 2 3],:);
      NumSZ = 3;
    end

    % set start and end times
    for kk =1:NumSZ
      SZstartPP(kk) = max(find(TimesSZ(kk,1) >= TimesAll(:,1)));
      SZendPP(kk)   = min(find(TimesSZ(kk,2) <= TimesAll(:,2)));
    end

    %% Generate EVC Arrays for Seizure, Pre, Post and Interictal
    clear SZstartQQ
    SZstartQQ = unique(SZstartPP);

    clear qNames*
    qNamesSZ   = [];
    qNamesPost = [];
    qNamesPre  = [];

    tic

    for hh = 1:length(SZstartQQ)
      pp = SZstartQQ(hh)
      clear EventNum SZnum OnsetNum SZevents All* NumWins
      NumWins = 0;

      clear FileName FileDir FilePath Data NumSamps NumChans NumSecs
      FilePath = char(InfoTime.filename(pp));
      FileDir  = FilePath(1:5);
      eval(['cd ',FileDir])

      hdr          = readhdr('eeg.hdr');
      num_channels = hdr.file_fmt.Numb_chans;
      format_file  = hdr.file_fmt.File_format;
      FileName     = FilePath(6:end);
      NumChans     = num_channels;

      %%% Read in the original part i) adjacency matrix to compute "deviation matrix"
      cd adj_matrices
      clear AdjMat

      if ~strcmp(mattype, 'q')
          disp('not q')
          eval(['AdjMat = read_adj_matrix(''adj_pwr_data_',FilePath(FPindex) ...
               ,'_',FreqBand,'.dat'',',num2str(NumChans),');'])
      else
          eval(['AdjMat = read_adj_matrix(''adj_chr_data_',FilePath(FPindex) ...
           ,'_',FreqBand,'.dat'',',num2str(NumChans),');'])
      end
      cd ..
      [NumSecs,n1,n2] = size(AdjMat);   % get array dimensions

      %%% Get the channels of interest 
      if sum(Patient == 'PY04N007') == 8 
        if     sum(FileDir(1:5) == 'rec01') == 5
          AdjMat2 = AdjMat(1:NumSecs,[1:6 8:21 24:27 29:31 33 35:40 42:51 54:56 58:63 65:79 81:87],[1:6 8:21 24:27 29:31 33 35:40 42:51 54:56 58:63 65:79 81:87]);
        end
      end

      if sum(Patient == 'PY04N008') == 8 
        if     sum(FileDir(1:5) == 'rec01') == 5
          AdjMat2 = AdjMat(1:NumSecs,[1:40],[1:40]);

            %%%%%% Only used when looking at focus electrodes
            if focuselectrodes == 1
              AdjMat2 = AdjMat(1:NumSecs,chans,chans);
            end 
        end
      end

      if sum(Patient == 'PY04N012') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:34 51:98],[1:34 51:98]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat;
        elseif sum(FileDir(1:4) == 'edf2') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:82],[1:82]);
        end
      end

      if strcmp(Patient, 'PY04N013')
        AdjMat2 = AdjMat(1:NumSecs,[1 4:10 12:61 63 65 70 72:89],[1 4:10 12:61 63 65 70 72:89]);

        %%%%%% Only used when looking at focus electrodes
        if focuselectrodes == 1
          AdjMat2 = AdjMat(1:NumSecs,chans,chans);
        end 
      end

      if sum(Patient == 'PY04N015') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:2 5:54 56:83 85:87 89:91 93:95],[1:2 5:54 56:83 85:87 89:91 93:95]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:2 5:54 56:83 85:87 89:91 93:95],[1:2 5:54 56:83 85:87 89:91 93:95]);
        end
      end

      if sum(Patient == 'PY05N004') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:12 14 17:37 39:63 65:82 84:100], ...
                                     [1:12 14 17:37 39:63 65:82 84:100]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:12 14 17:37 39:81 83:99], ...
                                     [1:12 14 17:37 39:81 83:99]);
        elseif sum(FileDir(1:4) == 'edf2') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:12 14 17:37 39:81 83:99], ...
                                     [1:12 14 17:37 39:81 83:99]);
        end

        %%%%%% Only used when looking at focus electrodes
        if focuselectrodes == 1
          AdjMat2 = AdjMat(1:NumSecs,chans,chans);
        end
      end

      if sum(Patient == 'PY05N005') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:21 23:38 41:83 89:118],[1:21 23:38 41:83 89:118]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:21 23:38 41:83 89:118],[1:21 23:38 41:83 89:118]);
        end
      end

      if sum(Patient == 'PY05N006') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:2 5:24],[1:2 5:24]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:2 5:24],[1:2 5:24]);
        end
      end

      if sum(Patient == 'PY05N007') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:2 5:7 10:14 17:48],[1:2 5:7 10:14 17:48]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:2 5:7 10:14 17:48],[1:2 5:7 10:14 17:48]);
        end
      end

      if sum(Patient == 'PY11N003') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:53 55:117],[1:53 55:117]);
      end

      if sum(Patient == 'PY11N004') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[2 4:15 17:31 33:106 108:117],[2 4:15 17:31 33:106 108:117]);
      end

      if sum(Patient == 'PY11N006') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:16 18:66],[1:16 18:66]);
      end

      if sum(Patient == 'PY11N009') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:22 24:30 32:81 85:90 92:105],[1:22 24:30 32:81 85:90 92:105]);
      end

      if sum(Patient == 'PY11N011') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[6:32],[6:32]);
      end

      if sum(Patient == 'PY11N014') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:21 23:40 42:46 48:52 55:68 70:81 83:111],[1:21 23:40 42:46 48:52 55:68 70:81 83:111]);
      end

      clear AdjMat
      AdjMat = AdjMat2;
      clear AdjMat2

      % initialize leading EVC vars
      clear Lead*
      LeadSingVals(1:NumSecs) = 0;
      LeadSingVecs(1:NumSecs,1:NumChans2) = 0;  %rows = seconds, cols = channels

      % loop through each time point (k) in part iii) EVC_k(i)
      for gg = 1:NumSecs
        clear U1 S1 V1 ADJtemp lambda
        lambda = log(2)/2;
        %%% Network Analysis Part ii)
        ADJtemp = squeeze(AdjMat(gg,1:NumChans2,1:NumChans2));
        ADJtemp = (ADJtemp - mAllAdjMat)./mAllAdjSTD;
        ADJtemp = exp(ADJtemp)./(1 + exp(ADJtemp)); % transform to [0,1] to retain "weighted connections" interpretation

        % set diagonals to 0 (electrodes don't connect to each other?)
        for bb = 1:NumChans2
          ADJtemp(bb,bb) = 0;
        end

        ADJtemp(isnan(ADJtemp)) = 0;    % set NaN to 0's

        %%%% perform singular value decomposition
        [U1,S1,V1] = svd(ADJtemp);
        LeadSingVals(gg) = S1(1,1);                             % leading eigenvalues
        LeadSingVecs(gg,1:NumChans2) = abs(U1(1:NumChans2,1));  % leading eigenvectors EVC_k(i)
      end

      clear AdjMat ADJtemp
      cd ..

      clear SZevents
      SZevents = find(SZstartPP == pp)

      clear(strcat('q',num2str(pp),'*'));
      % not really a loop, just finding the seizure
      for ff = 1:length(SZevents)
        %%% get the start and end times of the seizure event
        clear SZstart SZstop SZlength L1 L2 EventNum
        EventNum = SZevents(ff);
        SZstart  = floor((TimesSZ(EventNum,1) - TimesAll(pp,1)));
        SZstop   = ceil((TimesSZ(EventNum,2)  - TimesAll(pp,1)));

        if sum(Patient == 'PY11N011') == 8
          clear SZstart SZstop
          SZstart  = floor((TimesSZ(EventNum,1) - TimesAll(pp,1))/5);
          SZstop   = ceil((TimesSZ(EventNum,2)  - TimesAll(pp,1))/5);
        end

        % store seizure length
        SZlength = SZstop - SZstart + 1;
        [L1,L2] = size(LeadSingVecs);

        disp(['q' num2str(pp) ' is ' num2str(SZlength) ' seconds long'])
        %% Concatenate seizure, pre, and post ictal EVCs together from +/- 10 minutes
        %%% generating them in the form of q<pp>sz, q<pp>pre, q<pp>post where
        %%% pp represents the indice of the file with that seizure/post/pre

        if length(SZevents) == 1
          disp('Length of seize events is 1');
          eval(['q',num2str(pp),'sz = LeadSingVecs(SZstart:SZstop,:);'])
          qNamesSZ = [qNamesSZ [' q',num2str(pp),'sz''']];

          if (SZstart-timeWindow) >= 1
            disp('first if');
            eval(['q',num2str(pp),'pre  = LeadSingVecs(SZstart-timeWindow:SZstart-1,:);'])
            qNamesPre = [qNamesPre [' q',num2str(pp),'pre''']];
          elseif (SZstart-timeWindow) <  1
              disp('2 if');
            eval(['q',num2str(pp),'pre  = LeadSingVecs(1:SZstart-1,:);'])
            qNamesPre = [qNamesPre [' q',num2str(pp),'pre''']];
          end

          if (SZstop+timeWindow) <= L1
            eval(['q',num2str(pp),'post = LeadSingVecs(SZstop+1:SZstop+timeWindow,:);'])
            qNamesPost = [qNamesPost [' q',num2str(pp),'post''']];
          elseif (SZstop+timeWindow) >  L1
            eval(['q',num2str(pp),'post = LeadSingVecs(SZstop+1:end,:);'])
            qNamesPost = [qNamesPost [' q',num2str(pp),'post''']];
          end
        end

        if length(SZevents) > 1
          if ff == 1
            eval(['q',num2str(pp),'Asz = LeadSingVecs(SZstart:SZstop,:);'])
            qNamesSZ = [qNamesSZ [' q',num2str(pp),'Asz''']];

            if (SZstart-timeWindow) >= 1
              eval(['q',num2str(pp),'Apre  = LeadSingVecs(SZstart-timeWindow:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Apre''']];
            elseif (SZstart-timeWindow) <  1
              eval(['q',num2str(pp),'Apre  = LeadSingVecs(1:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Apre''']];
            end

            if (SZstop+timeWindow) <= L1
              eval(['q',num2str(pp),'Apost = LeadSingVecs(SZstop+1:SZstop+timeWindow,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Apost''']];
            elseif (SZstop+timeWindow) >  L1
              eval(['q',num2str(pp),'Apost = LeadSingVecs(SZstop+1:end,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Apost''']];
            end
          end

          if ff == 2
            eval(['q',num2str(pp),'Bsz = LeadSingVecs(SZstart:SZstop,:);'])
            qNamesSZ = [qNamesSZ [' q',num2str(pp),'Bsz''']];

            if (SZstart-timeWindow) >= 1
              eval(['q',num2str(pp),'Bpre  = LeadSingVecs(SZstart-timeWindow:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Bpre''']];
            elseif (SZstart-timeWindow) <  1
              eval(['q',num2str(pp),'Bpre  = LeadSingVecs(1:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Bpre''']];
            end

            if (SZstop+timeWindow) <= L1
              eval(['q',num2str(pp),'Bpost = LeadSingVecs(SZstop+1:SZstop+timeWindow,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Bpost''']];
            elseif (SZstop+timeWindow) >  L1
              eval(['q',num2str(pp),'Bpost = LeadSingVecs(SZstop+1:end,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Bpost''']];
            end
          end

          if ff == 3
            eval(['q',num2str(pp),'Csz = LeadSingVecs(SZstart:SZstop,:);'])
            qNamesSZ = [qNamesSZ [' q',num2str(pp),'Csz''']];

            if (SZstart-timeWindow) >= 1
              eval(['q',num2str(pp),'Cpre  = LeadSingVecs(SZstart-timeWindow:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Cpre''']];
            elseif (SZstart-timeWindow) <  1
              eval(['q',num2str(pp),'Cpre  = LeadSingVecs(1:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Cpre''']];
            end

            if (SZstop+timeWindow) <= L1
              eval(['q',num2str(pp),'Cpost = LeadSingVecs(SZstop+1:SZstop+timeWindow,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Cpost''']];
            elseif (SZstop+timeWindow) >  L1
              eval(['q',num2str(pp),'Cpost = LeadSingVecs(SZstop+1:end,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Cpost''']];
            end
          end

          if ff == 4
            eval(['q',num2str(pp),'Dsz = LeadSingVecs(SZstart:SZstop,:);'])
            qNamesSZ = [qNamesSZ [' q',num2str(pp),'Dsz''']];

            if (SZstart-timeWindow) >= 1
              eval(['q',num2str(pp),'Dpre  = LeadSingVecs(SZstart-timeWindow:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Dpre''']];
            elseif (SZstart-timeWindow) <  1
              eval(['q',num2str(pp),'Dpre  = LeadSingVecs(1:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Dpre''']];
            end

            if (SZstop+timeWindow) <= L1
              eval(['q',num2str(pp),'Dpost = LeadSingVecs(SZstop+1:SZstop+timeWindow,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Dpost''']];
            elseif (SZstop+timeWindow) >  L1
              eval(['q',num2str(pp),'Dpost = LeadSingVecs(SZstop+1:end,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Dpost''']];
            end
          end

          if ff == 5
            eval(['q',num2str(pp),'Esz = LeadSingVecs(SZstart:SZstop,:);'])
            qNamesSZ = [qNamesSZ [' q',num2str(pp),'Esz''']];

            if (SZstart-timeWindow) >= 1
              eval(['q',num2str(pp),'Epre  = LeadSingVecs(SZstart-timeWindow:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Epre''']];
            elseif (SZstart-timeWindow) <  1
              eval(['q',num2str(pp),'Epre  = LeadSingVecs(1:SZstart-1,:);'])
              qNamesPre = [qNamesPre [' q',num2str(pp),'Epre''']];
            end

            if (SZstop+timeWindow) <= L1
              eval(['q',num2str(pp),'Epost = LeadSingVecs(SZstop+1:SZstop+timeWindow,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Epost''']];
            elseif (SZstop+timeWindow) >  L1
              eval(['q',num2str(pp),'Epost = LeadSingVecs(SZstop+1:end,:);'])
              qNamesPost = [qNamesPost [' q',num2str(pp),'Epost''']];
            end
          end
        end

        %%% save the different sections of the EVC before concatenating 
        currentdir = pwd
        if length(SZevents) > 1
            disp('SZevents is > 1')
            switch ff
                case 1
                    qprelength = 0;
                    qpostlength = 0;
                    eval(['qprelength = size((q' num2str(pp) 'Apre))'])
                    eval(['qpostlength = size((q' num2str(pp) 'Apost))'])
                    numseizures = 1;
                case 2
                    eval(['qprelength = size((q' num2str(pp) 'Bpre))'])
                    eval(['qpostlength = size((q' num2str(pp) 'Bpost))'])
                    numseizures = numseizures + 1;
                case 3
                    eval(['qprelength = size((q' num2str(pp) 'Cpre))'])
                    eval(['qpostlength = size((q' num2str(pp) 'Cpost))'])
                    numseizures = numseizures + 1;
                case 4
                    eval(['qprelength = size((q' num2str(pp) 'Dpre))'])
                    eval(['qpostlength = size((q' num2str(pp) 'Dpost))'])
                    numseizures = numseizures + 1;
                case 5
                    eval(['qprelength = size((q' num2str(pp) 'Epre))'])
                    eval(['qpostlength = size((q' num2str(pp) 'Epost))'])
                    numseizures = numseizures + 1;
                otherwise
                    disp('WTF at switch statement?')
            end
        else 
            numseizures = 0;
            eval(['qprelength = size((q' num2str(pp) 'pre));'])
            eval(['qpostlength = size((q' num2str(pp) 'post));'])
        end

        index = 1;
        filename = {};
        while qprelength(1) < timeWindow
            cd(currentdir)

            [len, testfile] = helpSaveQ(pp-index, InfoTime, homedir, mattype, FPindex, FreqBand, Patient, NumChans2, mAllAdjMat, mAllAdjSTD, focuselectrodes);
            qprelength(1) = qprelength(1) + len; % add the size of each recording session
            
            filename = [filename; testfile];
            filename = filename(~cellfun(@isempty, filename));
            
            index = index+1;
        end

%         index = 1;
%         while qpostlength(1) < timeWindow
%             cd(currentdir)
% 
%             len = helpSaveQ(pp+index, InfoTime, homedir, mattype, FPindex, FreqBand, Patient, NumChans2, mAllAdjMat, mAllAdjSTD, focuselectrodes);
%             qpostlength(1) = qpostlength(1) + len; % add the size of each recording session
% 
%             index = index+1;
%         end

          %% Save the different qpre, sz and post files 
        cd(fullfile(homedir))    % cd to storage place of the EVCs

%          switch numseizures
%             case 0
%                 disp('case 0!')
%                 % save the separate EVC re)gions to plot separately
%                 save(strcat(mattype,num2str(pp),'post_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'post'))
%                 save(strcat(mattype,num2str(pp),'pre_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'pre'))
%                 save(strcat(mattype,num2str(pp),'sz_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'sz'))
%                 
%                 % store file names
%                 files = {{strcat(mattype,num2str(pp),'post_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'post')}; ...
%                             {strcat(mattype,num2str(pp),'pre_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'pre')}; ...
%                             {strcat(mattype,num2str(pp),'sz_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'sz')}};
%                 
%                 filename = [filename; files];
%                 filename = filename(~cellfun(@isempty, filename));
%             case 1
%                 % save the separate EVC re)gions to plot separately
%                 save(strcat(mattype,num2str(pp),'Apost_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Apost'))
%                 save(strcat(mattype,num2str(pp),'Apre_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Apre'))
%                 save(strcat(mattype,num2str(pp),'Asz_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Asz'))
%             case 2
%                 % save the separate EVC re)gions to plot separately
%                 save(strcat(mattype,num2str(pp),'Bpost_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Bpost'))
%                 save(strcat(mattype,num2str(pp),'Bpre_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Bpre'))
%                 save(strcat(mattype,num2str(pp),'Bsz_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Bsz'))
%             case 3
%                 % save the separate EVC re)gions to plot separately
%                 save(strcat(mattype,num2str(pp),'Cpost_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Cpost'))
%                 save(strcat(mattype,num2str(pp),'Cpre_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Cpre'))
%                 save(strcat(mattype,num2str(pp),'Csz_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Csz'))
%             case 4
%                 % save the separate EVC re)gions to plot separately
%                 save(strcat(mattype,num2str(pp),'Dpost_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Dpost'))
%                 save(strcat(mattype,num2str(pp),'Dpre_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Dpre'))
%                 save(strcat(mattype,num2str(pp),'Dsz_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Dsz'))
%             case 5
%                 % save the separate EVC re)gions to plot separately
%                 save(strcat(mattype,num2str(pp),'Epost_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Epost'))
%                 save(strcat(mattype,num2str(pp),'Epre_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Epre'))
%                 save(strcat(mattype,num2str(pp),'Esz_', Patient, '_', FreqBand, num2str(timeWindow)), strcat('q',num2str(pp),'Esz'))
%             otherwise
%                 disp('WTF at switch statement?')
%         end
        cd(currentdir)
      end
      
      %% Do mass saving and concatenating the matrices together
      cd(homedir)
      % loop through filenames and concatenate (only per freq. band)
      eval(['q' num2str(pp) 'pre_' Patient FreqBand ' = [];'])
      for ff=1:length(filename)
           load(filename{end-ff+1})
           eval(['q' num2str(pp) 'pre_' Patient FreqBand ' = [q' num2str(pp) 'pre_' Patient FreqBand ';' filename{end-ff+1} '];'])
%            clear(filename{end-ff+1})
      end
        
      % save the new timewindow qpre
      eval(['save q' num2str(pp) 'pre_' Patient FreqBand ' q' num2str(pp) 'pre_' Patient FreqBand])
      cd(currentdir)
    end
    

    
    %% Save Concatenated Arrays for All EVC Time Periods
    cd(homedir)

%     clear qAll*
%     eval(['qAllSZ   = [',qNamesSZ,']'';'])
%     eval(['qAllPre  = [',qNamesPre,']'';'])
%     eval(['qAllPost = [',qNamesPost,']'';'])
% 
%     eval(['save ', mattype, 'AllSZ_', Patient, '_', FreqBand, num2str(timeWindow), ' qAllSZ'])
%     eval(['save ', mattype, 'AllPre_', Patient, '_', FreqBand, num2str(timeWindow), ' qAllPre'])
%     eval(['save ', mattype, 'AllPost_', Patient, '_', FreqBand, num2str(timeWindow), ' qAllPost'])
end


%%%% If I want to add Freq Band to File names
%     cd(fullfile(homedir))    % cd to storage place of the EVCs
%     % save the separate EVC re)gions to plot separately
%     save(strcat(mattype,num2str(pp),'post_', Patient, '_', FreqBand), strcat('q',num2str(pp),'post'))
%     save(strcat(mattype,num2str(pp),'pre_', Patient, '_', FreqBand), strcat('q',num2str(pp),'pre'))
%     save(strcat(mattype,num2str(pp),'sz_', Patient, '_', FreqBand), strcat('q',num2str(pp),'sz'))
% 
%     cd(currentdir)
%   end
% end

% %% Save Concatenated Arrays for All EVC Time Periods
% cd(homedir)
% 
% clear qAll*
% eval(['qAllSZ   = [',qNamesSZ,']'';'])
% eval(['qAllPre  = [',qNamesPre,']'';'])
% eval(['qAllPost = [',qNamesPost,']'';'])
% 
% eval(['save ', mattype, 'AllSZ_', Patient, '_', FreqBand, ' qAllSZ'])
% eval(['save ', mattype, 'AllPre_', Patient, '_', FreqBand, ' qAllPre'])
% eval(['save ', mattype, 'AllPost_', Patient, '_', FreqBand, ' qAllPost'])
