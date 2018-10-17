%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original File written by Samuel Burns
% Adapted by: Adam Li
%
% Ver.: 1.0 - Date: 08/14/2015
%
% Description: To compute average interictal connectivity
% 
% Dependencies: i) readhdr.m - read header file
%              ii) infoevent.mat - a mat file containing information about
%              seizure events during the recording sessions
%              iii) infotime.mat - a mat file containing information about
%              recording sessions
%
% Output: Saving of standard dev. of adj. matrix for a certain patient 
% such as mAllAdjSTD_<Patient>_<FreqBand>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *** Used for Analyzing the Focus Only 
% PY04N008:
% chans = [1, 2, 3, 9, 10, 11];
% NumChans2 = 6;
%
% % PY04N013
% chans = [5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75];
% NumChans2 = 44;
%
% % PY05N004
% chans = [80,81,82,83,84,85,88,89,90,91,92,93];
% NumChans2 = 12;
%% 01: Settings for Running File and Adding Directory Paths and Loading in Mat Files
clear Pat* Info* Num* n2 *Length OverLap sfreq stimfreq FreqBand
Patient = 'PY05N004';
Patients = {'PY04N008', 'PY04N013', 'PY05N004'};

% add utility software functions for helping read in data 
addpath('/home/adamli/MATLAB/code_santaniello/utility software/generate_pbs')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/read_data')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/svd_eeg')
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')

for p=1:length(Patients)
    Patient = Patients{p};
    
    % go to the directory with the mat files that hold information about
    % seizure events, and recording sessions
    homedir = '/home/adamli/MATLAB/code_adam/burns_adapted'
    cd(homedir);
    eval(['load(''infoevent.mat'',''event',Patient,''');'])
    eval(['InfoEvent = event',Patient,';'])
    eval(['load(''infotime.mat'',''',Patient,''');'])
    eval(['InfoTime = ',Patient,';'])

    PatTimes = InfoTime.time;       % load patient recording settings
    [NumFiles,n2] = size(PatTimes); % obtain size of recording time's array
    %eval(['cd /media/ExtHDD01/',Patient])
    %eval(['cd /home/sburns/Data/JHU/JHU2012/',Patient])

    clear *SZ* s2 TimesAll
    [NumSZ,s2] = size(InfoEvent.code);
    TimesAll = InfoTime.time;
    TimesSZ  = InfoEvent.time;

    if sum(Patient == 'PY11N014') == 8
      TimesAll(14:end,:) = TimesAll(14:end,:) + (24*3600)*19;
      TimesSZ(5:end,:)   = TimesSZ(5:end,:)   + (24*3600)*30;
    end

    % loop through # of seizures there are
    for kk =1:NumSZ
        % set seizure start and end time
      SZstartPP(kk) = max(find(TimesSZ(kk,1) >= TimesAll(:,1)));
      SZendPP(kk)   = min(find(TimesSZ(kk,2) <= TimesAll(:,2)));
    end

    %% 02: Determine Patient Settings 
    clear All* NumWins NumChans2 AllAdjMat mAllAdjMat FPindex
    NumWins = 0;        %# of windows
    sfreq = 1000;       %sampling frequency

    % go through each PY11* patient and set:
    % 1. number of channels
    % 2. ???? FPindex
    % 3. the main frequency band from R-spec peak freq. band
    % 4. cd to that patient's directory
    if     sum(Patient == 'PY04N007') == 8
      NumChans2 = 75;
      FPindex = 12:16;
      FreqBand = 'beta';
      HDDrive = 'ExtHDD04';
    elseif sum(Patient == 'PY04N008') == 8
      NumChans2 = 40;
      FPindex = 12:16;
      FreqBand = 'beta';
    elseif sum(Patient == 'PY04N012') == 8
      NumChans2 = 82;
      FPindex = 11:15;
      FreqBand = 'beta';
      HDDrive = 'ExtHDD04';
    elseif sum(Patient == 'PY04N013') == 8
      NumChans2 = 79;
      FPindex = 12:16;
      FreqBand = 'beta';
    elseif sum(Patient == 'PY04N015') == 8
      NumChans2 = 89;
      FPindex = 11:15;
      FreqBand = 'beta';
    elseif sum(Patient == 'PY05N004') == 8
      NumChans2 = 94;
      FPindex = 11:15;
      FreqBand = 'alpha';
    elseif sum(Patient == 'PY05N005') == 8
      NumChans2 = 110;
      FPindex = 11:15;
      FreqBand = 'beta';
    elseif sum(Patient == 'PY05N006') == 8
      NumChans2 = 22;
      FPindex = 11:15;
      FreqBand = 'theta';
    elseif sum(Patient == 'PY05N007') == 8
      NumChans2 = 42;
      FPindex = 11:15;
      FreqBand = 'beta';
    end

    if sum(Patient == 'PY11N003') == 8
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
      NumChans2 = 101;
      FPindex = 11:15;
      FreqBand = 'beta';
      sfreq = 1000;
      eval(['cd /media/ExtHDD04/',Patient])
    end

    % go to that patient's directory
    % eval(['cd /media/',HDDrive,'/',Patient])
    eval(['cd /home/adamli/Desktop/', Patient])

    %%%%%%%%
    FreqBand = 'high'; % alpha, beta, gamma, high
    %%%%%%%%

      %% *** For looking at focus electrodes/channels only
      % PY04N008:
      focuselectrodes = 0;
    %   chans = [1, 2, 3, 9, 10, 11];
    %   NumChans2 = 6;
    %   focuselectrodes = 1;

    % %  % PY04N013
    %   chans = [5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75];
    %   NumChans2 = 44;
    %   focuselectrodes = 1;
    % %   
    % % %  PY05N004
%     chans = [80,81,82,83,84,85,88,89,90,91,92,93];
%     NumChans2 = 12;
%     focuselectrodes = 1;

    %% 03: FFT Settings
    WinLength = 3*sfreq;
    FFTlength = sfreq;
    Overlap = 1/3; 
    NumFFTs = WinLength/(Overlap*FFTlength) - 2;
    stimfreq = sfreq*(0:FFTlength/2)/FFTlength;
    stimfreq = stimfreq(1:sfreq/2);

    AllAdjMat(1:NumChans2,1:NumChans2) = 0;

    disp(['Running through patient ', Patient]);
    disp(['Main freq band: ', FreqBand]);
    disp(['Number of files: ', num2str(NumFiles)]);

    %% 03: Read in Adjacency Matrices and Perform Processing, to get them into averaged form
    for pp = 1:NumFiles 
      disp(strcat('On loop #: ', num2str(pp)));
      clear FileName FileDir FilePath Data NumSamps NumSecs
      FilePath = char(InfoTime.filename(pp));
      FileDir  = FilePath(1:5);
      eval(['cd ',FileDir])

      % Set up Headers
      hdr          = readhdr('eeg.hdr');
      num_channels = hdr.file_fmt.Numb_chans;
      format_file  = hdr.file_fmt.File_format;
      FileName     = FilePath(6:end);
      NumChans     = num_channels;

      cd adj_matrices
      clear AdjMat
      eval(['AdjMat = read_adj_matrix(''adj_chr_data_',FilePath(FPindex) ...
           ,'_',FreqBand,'.dat'',',num2str(NumChans),');'])
      cd ..
      [NumSecs,n1,n2] = size(AdjMat);


      %% reset adjmat
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

      if sum(Patient == 'PY04N013') == 8 
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
      end  %}

      if sum(Patient == 'PY11N005') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:11 13:33 35:74],[1:11 13:33 35:74]);
      end

      if sum(Patient == 'PY11N006') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:16 18:66],[1:16 18:66]);
      end

      %%%%%% added by Adam %
      if sum(Patient == 'PY11N007') == 8          % PY11N007
        AdjMat2 = AdjMat(1:NumSecs,[1:5 7:18 20:33 35:48],[1:5 7:18 20:33 35:48]);
      end
      if sum(Patient == 'PY11N008') == 8          % PY11N008
        AdjMat2 = AdjMat(1:NumSecs,[1:113 115:126],[1:113 115:126]);
      end
      %%%%%% end of Adam %

      if sum(Patient == 'PY11N009') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:22 24:30 32:81 85:90 92:105],[1:22 24:30 32:81 85:90 92:105]);
      end

      if sum(Patient == 'PY11N010') == 8          % PY11NO10
        AdjMat2 = AdjMat(1:NumSecs,[1:77 79:85 87:88],[1:77 79:85 87:88]);
      end

      if sum(Patient == 'PY11N011') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[6:32],[6:32]);
      end
      %%%%%% added by Adam %
      if sum(Patient == 'PY11N012') == 8          % PY11N012
        AdjMat2 = AdjMat(1:NumSecs,[1:55 57:89 91:110],[1:55 57:89 91:110]);
      end
      if sum(Patient == 'PY11N013') == 8          % PY11N013
        AdjMat2 = AdjMat(1:NumSecs,[1:111 113:119],[1:111 113:119]);
      end
      %%%%%% end of Adam % 
      if sum(Patient == 'PY11N014') == 8 
        AdjMat2 = AdjMat(1:NumSecs,[1:21 23:28 31:39 42:46 48:52 55:68 70:81 83:111],[1:21 23:28 31:39 42:46 48:52 55:68 70:81 83:111]);
        AdjMat2 = AdjMat(1:NumSecs,[1:21 23:40 42:46 48:52 55:68 70:81 83:111],[1:21 23:40 42:46 48:52 55:68 70:81 83:111]);
      end

      clear AdjMat
      AdjMat = AdjMat2;
      clear AdjMat2

      %% 04: Look at Seizures
      clear SZstarts SZstops TotalPP NumSZs
      NumSZs = sum(pp == SZstartPP);
      TotalPP(1:NumSecs) = 1;

      % if # of seizures > 0
      if NumSZs > 0
        SZind = find(SZstartPP == pp); % find start point
        for kk = 1:NumSZs % loop through each seizure
           SZstart(kk) = floor((TimesSZ(SZind(kk),1) - TimesAll(pp,1)));
           SZstop(kk)  = floor((TimesSZ(SZind(kk),2) - TimesAll(pp,1)));
           TotalPP(SZstart(kk):SZstop(kk)) = 0;
        end
      end

      for gg = 1:NumSecs   
        if TotalPP(gg) == 1
          %create temporary holder for adjacency matrix
          clear ADJtemp
          ADJtemp = squeeze(AdjMat(gg,1:NumChans2,1:NumChans2));

          for bb = 1:NumChans2
            ADJtemp(bb,bb) = 0;
          end

          if sum(sum(isnan(ADJtemp))) == 0
            AllAdjMat = AllAdjMat + ADJtemp; %update adj. matrix
            NumWins = NumWins + 1;           %update # windows
          end
        end
      end
      cd .. 
    end

    %% 05: Saving Average Adj. Matrix for Specific Freq. Band
    %Determine average adj. matrix
    mAllAdjMat = AllAdjMat/NumWins;

    cd(homedir) % go back to the main directory
    eval(['save mAllAdjMat_',Patient,'_',FreqBand,' mAllAdjMat'])
end