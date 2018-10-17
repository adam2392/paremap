%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original File written by Samuel Burns
% Adapted by: Adam Li
%
% Ver.: 1.0 - Date: 08/14/2015
%
% Description: To generate array of eigenvectors(EVC) of all interictal
% time perdiods in an array 'AllLeadSVD'.
% 
% Dependencies: i) readhdr.m - read header file
%              ii) infoevent.mat - a mat file containing information about
%              seizure events during the recording sessions
%              iii) infotime.mat - a mat file containing information about
%              recording sessions
%
% Output: To save the AllLeadSVD_<Patient> into a mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Pat* Info* Num* n2 *Length OverLap sfreq stimfreq
Patient = 'PY04N012';
sfreq = 1000;
mattype = 'q'; % q or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_';
else
    mat = '';
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
  sfreq = 200;
elseif sum(Patient == 'PY11N014') == 8
  FreqBand = 'beta';
end

% add utility software functions for helping read in data 
addpath('/home/adamli/MATLAB/code_santaniello/utility software/generate_pbs')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/read_data')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/svd_eeg')
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')
addpath('/home/adamli/MATLAB/code_adam/')

% go to the directory with the mat files that hold information about
% seizure events, and recording sessions, and avg. Adj. Matrix
homedir = '/home/adamli/MATLAB/code_adam/burns_adapted'
cd(homedir)
eval(['load(''infoevent.mat'',''event',Patient,''');'])
eval(['InfoEvent = event',Patient,';'])
eval(['load(''infotime.mat'',''',Patient,''');'])
eval(['InfoTime = ',Patient,';'])
eval(['load ', mat, 'mAllAdjMat_',Patient,'_',FreqBand])
eval(['load ', mat, 'mAllAdjSTD_',Patient,'_',FreqBand])

PatTimes = InfoTime.time;       % load patient recording settings
[NumFiles,n2] = size(PatTimes); % obtain size of recording time's array

%% Determine Patient Settings
clear All* NumWins NumChans* HDDrive FreqBand FPindex

% go through each patient and set:
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
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY04N012') == 8
  NumChans2 = 82;
  FPindex = 11:15;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD04';
elseif sum(Patient == 'PY04N013') == 8
  NumChans2 = 79;
  FPindex = 12:16;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD02';
elseif sum(Patient == 'PY04N015') == 8
  NumChans2 = 89;
  FPindex = 11:15;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD02';
elseif sum(Patient == 'PY05N004') == 8
  NumChans2 = 94;
  FPindex = 11:15;
  FreqBand = 'alpha';
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY05N005') == 8
  NumChans2 = 110;
  FPindex = 11:15;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY05N006') == 8
  NumChans2 = 22;
  FPindex = 11:15;
  FreqBand = 'theta';
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY05N007') == 8
  NumChans2 = 42;
  FPindex = 11:15;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD01';
end

if     sum(Patient == 'PY11N003') == 8
  NumChans2 = 116;
  FPindex = 11:15;
  FreqBand = 'gamma';
  HDDrive = 'ExtHDD03';
elseif sum(Patient == 'PY11N004') == 8
  NumChans2 = 112;
  FPindex = 11:15;
  FreqBand = 'gamma';
  HDDrive = 'ExtHDD03';
elseif sum(Patient == 'PY11N006') == 8
  NumChans2 = 65;
  FPindex = 11:15;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD03';
elseif sum(Patient == 'PY11N009') == 8
  NumChans2 = 99;
  FPindex = 11:15;
  FreqBand = 'theta';
  HDDrive = 'ExtHDD03';
elseif sum(Patient == 'PY11N011') == 8
  NumChans2 = 27;
  FPindex = 11:15;
  FreqBand = 'theta';
  sfreq = 200;
  HDDrive = 'ExtHDD04';
elseif sum(Patient == 'PY11N014') == 8
  NumChans2 = 101;
  FPindex = 11:15;
  FreqBand = 'beta';
  HDDrive = 'ExtHDD04';
end

%% Set Spectrum Settings
WinLength = 3*sfreq;
FFTlength = sfreq;
Overlap = 1/3; 
NumFFTs = WinLength/(Overlap*FFTlength) - 2;
stimfreq = sfreq*(0:FFTlength/2)/FFTlength;
stimfreq = stimfreq(1:sfreq/2);

clear *SZ* s2 TimesAll
[NumSZ,s2] = size(InfoEvent.code);
TimesAll = InfoTime.time;
TimesSZ  = InfoEvent.time;

if sum(Patient == 'PY11N014') == 8
  TimesAll(14:end,:) = TimesAll(14:end,:) + (24*3600)*19;
  TimesSZ(5:end,:)   = TimesSZ(5:end,:)   + (24*3600)*30;
end

for kk =1:NumSZ
  SZstartPP(kk) = max(find(TimesSZ(kk,1) >= TimesAll(:,1)));
  SZendPP(kk)   = min(find(TimesSZ(kk,2) <= TimesAll(:,2)));
end

% cd to the patient's directory
eval(['cd /media/',HDDrive,'/',Patient])
NumWins = 0;
disp(strcat('Running through patient ', Patient));
disp(strcat('Main freq band: ', FreqBand));
disp(strcat('Number of files: ', num2str(NumFiles)));

%% loop through each recording file...
for pp = 1:NumFiles 
  disp(['On file # ' num2str(pp) ' out of ' num2str(NumFiles)]);
  clear FileName FileDir FilePath Data NumSamps NumChans NumSecs
  FilePath = char(InfoTime.filename(pp));
  FileDir  = FilePath(1:5);
  eval(['cd ',FileDir])

  %set up headers
  hdr          = readhdr('eeg.hdr');
  num_channels = hdr.file_fmt.Numb_chans;
  format_file  = hdr.file_fmt.File_format;
  FileName     = FilePath(6:end);
  NumChans     = num_channels;

  %grab adjacency matrices from folder in patient data folder in hard drive
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
  cd ../..
  [NumSecs,n1,n2] = size(AdjMat);

  %%%% get the channels of interest
  if sum(Patient == 'PY04N007') == 8 
    if     sum(FileDir(1:5) == 'rec01') == 5
      AdjMat2 = AdjMat(1:NumSecs,[1:6 8:21 24:27 29:31 33 35:40 42:51 54:56 58:63 65:79 81:87],[1:6 8:21 24:27 29:31 33 35:40 42:51 54:56 58:63 65:79 81:87]);
    end
  end

  if sum(Patient == 'PY04N008') == 8 
    if     sum(FileDir(1:5) == 'rec01') == 5
      AdjMat2 = AdjMat(1:NumSecs,[1:40],[1:40]);
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
  
  if sum(Patient == 'PY11N005') == 8
    AdjMat2 = AdjMat(1:NumSecs,[1:11 13:33 35:74],[1:11 13:33 35:74]);
  end
  
  if sum(Patient == 'PY11N006') == 8 
    AdjMat2 = AdjMat(1:NumSecs,[1:16 18:66],[1:16 18:66]);
  end
  
  %%%%%%% added by Adam %
  if sum(Patient == 'PY11N007') == 8          % PY11N007
    AdjMat2 = AdjMat(1:NumSecs,[1:5 7:18 20:33 35:48],[1:5 7:18 20:33 35:48]);
  end
  if sum(Patient == 'PY11N008') == 8          % PY11N008
    AdjMat2 = AdjMat(1:NumSecs,[1:113 115:126],[1:113 115:126]);
  end
  %%%%%%% end of Adam %
  
  if sum(Patient == 'PY11N009') == 8 
    AdjMat2 = AdjMat(1:NumSecs,[1:22 24:30 32:81 85:90 92:105],[1:22 24:30 32:81 85:90 92:105]);
  end
  
  if sum(Patient == 'PY11N010') == 8          % PY11NO10
    AdjMat2 = AdjMat(1:NumSecs,[1:77 79:85 87:88],[1:77 79:85 87:88]);
  end
  
  if sum(Patient == 'PY11N011') == 8 
    AdjMat2 = AdjMat(1:NumSecs,[6:32],[6:32]);
  end
  
  %%%%%%% added by Adam %
  if sum(Patient == 'PY11N012') == 8          % PY11N012
    AdjMat2 = AdjMat(1:NumSecs,[1:55 57:89 91:110],[1:55 57:89 91:110]);
  end
  if sum(Patient == 'PY11N013') == 8          % PY11N013
    AdjMat2 = AdjMat(1:NumSecs,[1:111 113:119],[1:111 113:119]);
  end
  %%%%%%% end of Adam % 
  
  if sum(Patient == 'PY11N014') == 8 
    AdjMat2 = AdjMat(1:NumSecs,[1:21 23:28 31:39 42:46 48:52 55:68 70:81 83:111],[1:21 23:28 31:39 42:46 48:52 55:68 70:81 83:111]);
  end

  % reset adjacency matrix to only get the channel's we are interested in
  clear AdjMat
  AdjMat = AdjMat2;
  clear AdjMat2
  
  % get the leading singular values and vectors for all interested channels
  clear Lead*
  LeadSingVals(1:NumSecs) = 0;
  LeadSingVecs(1:NumSecs,1:NumChans2) = 0;

  % loop through each time point in adjacency matrix
  for gg = 1:NumSecs
    clear U1 S1 V1 ADJtemp lambda
    lambda = log(2)/2;  %????
    
    %%% Network analysis part ii)
    ADJtemp = squeeze(AdjMat(gg,1:NumChans2,1:NumChans2));
    ADJtemp = (ADJtemp - mAllAdjMat)./mAllAdjSTD;
    %ADJtemp = exp(lambda.*ADJtemp);
    %ADJtemp = tanh(ADJtemp) + 1;
    ADJtemp = exp(ADJtemp)./(1 + exp(ADJtemp)); % transformation onto [0,1]

    % set diagonals to 0
    for bb = 1:NumChans2
      ADJtemp(bb,bb) = 0;
    end

    ADJtemp(isnan(ADJtemp)) = 0;
    
    %%% perform SVD calculation
    [U1,S1,V1] = svd(ADJtemp);
    LeadSingVals(gg) = S1(1,1);
    LeadSingVecs(gg,1:NumChans2) = abs(U1(1:NumChans2,1));
  end

  clear AdjMat ADJtemp
  clear SZstarts SZstops TotalPP NumSZs
  [NumSecs,NumChans2] = size(LeadSingVecs);  
  NumSZs = sum(pp == SZstartPP);
  TotalPP(1:NumSecs) = 1;

  % find the seizu  eval(['AdjMat = read_adj_matrix(''adj_chr_data_',FilePath(FPindex) ...
  %     ,'_',FreqBand,'.dat'',',num2str(NumChans),');'])re zones in this specific file... to get SVDs around it
  if NumSZs > 0
    SZind = find(SZstartPP == pp);
    for kk = 1:NumSZs
       SZstart(kk) = floor((InfoEvent.time(SZind(kk),1) - InfoTime.time(pp,1)));
       SZstop(kk)  = floor((InfoEvent.time(SZind(kk),2) - InfoTime.time(pp,1)));
       TotalPP(SZstart(kk):SZstop(kk)) = 0;
    end
  end

  % go through entire recording session, and grab the interictal time
  % periods
  for gg = 1:NumSecs   
    % save All the Leading SVD vectors
    if TotalPP(gg) == 1
      clear SVDtemp s1 s2
      NumWins = NumWins + 1;
      AllLeadSVD(NumWins,1:NumChans2) = squeeze(LeadSingVecs(gg,1:NumChans2));
    end
  end
end
 
%% Save the AllLeadSVD_<Patient>
cd(homedir)
qAllII = AllLeadSVD;
% eval(['save ', mat, 'AllLeadSVD_',Patient,' AllLeadSVD'])
eval(['save ', mat, 'qAllII_', Patient, ' qAllII'])