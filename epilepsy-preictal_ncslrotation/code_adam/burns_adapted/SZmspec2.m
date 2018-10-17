%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original File written by Samuel Burns
% Adapted by: Adam Li
%
% Ver.: 1.0 - Date: 08/16/2015
%
% Description: To generate seizure power spectrum as well as compute
% R-spectrum
% 
% Dependencies: i) readhdr.m - read header file
%              ii) infoevent.mat - a mat file containing information about
%              seizure events during the recording sessions
%              iii) infotime.mat - a mat file containing information about
%              recording sessions
%
% Output: Saving of ictal power spectrum for a certain patient such as
% SZmspec_<Patient>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 00: Setup Directories and Vars
clear Pat* Info* Num* n2 *Length OverLap sfreq stimfreq
Patient = 'PY11N009';

% go to the directory with the mat files that hold information about
% seizure events, and recording sessions
homedir = '/home/adamli/MATLAB/code_adam/burns_adapted'
cd(homedir)
eval(['load(''infoevent.mat'',''event',Patient,''');'])
eval(['InfoEvent = event',Patient,';'])
eval(['load(''infotime.mat'',''',Patient,''');'])
eval(['InfoTime = ',Patient,';'])

% add utility software functions for helping read in data 
addpath('/home/adamli/MATLAB/code_santaniello/utility software/generate_pbs')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/read_data')
addpath('/home/adamli/MATLAB/code_santaniello/utility software/svd_eeg')

PatTimes = InfoTime.time;
[NumFiles,n2] = size(PatTimes); 

clear *SZ* s2 TimesAll
[NumSZ,s2] = size(InfoEvent.code);
TimesAll = InfoTime.time;
TimesSZ  = InfoEvent.time;

if sum(Patient == 'PY11N014') == 8
  TimesAll(14:end,:) = TimesAll(14:end,:) + (24*3600)*19;
  TimesSZ(5:end,:)   = TimesSZ(5:end,:)   + (24*3600)*30;

  TimesSZ = TimesSZ([5 6 7 11 12 14 15],:);
  NumSZ = 7;
end

% loop through # of seizures and detect start and end times
for kk =1:NumSZ
  SZstartPP(kk) = max(find(TimesSZ(kk,1) >= TimesAll(:,1)));
  SZendPP(kk)   = min(find(TimesSZ(kk,2) <= TimesAll(:,2)));

  %clear S1temp S2temp
  %S1temp = find(TimesSZ(kk,1) >= TimesAll(:,1));
  %S2temp = find(TimesSZ(kk,2) <= TimesAll(:,2));
  %SZstartPP(kk) = intersect(S1temp,S2temp);
  %SZendPP(kk)   = SZstartPP(kk);
end

%% 01: Determine Patient Settings  
clear All* NumWins NumChans2 AllAdjMat mAllAdjMat FPindex
NumWins = 0;        %# of windows for spectral analysis
sfreq = 1000;       %set default sampling freq.

% go through each PY11* patient and set:
% 1. number of channels
% 2. sampling frequency
% 3. cd to that patient's directory
if     sum(Patient == 'PY04N007') == 8
  NumChans2 = 75; %NumChans2 = 86;
  HDDrive = 'ExtHDD04';
elseif sum(Patient == 'PY04N008') == 8
  NumChans2 = 31;
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY04N012') == 8
  NumChans2 = 82;
  HDDrive = 'ExtHDD04';
elseif sum(Patient == 'PY04N013') == 8
  NumChans2 = 79;
  HDDrive = 'ExtHDD02';
elseif sum(Patient == 'PY04N015') == 8
  NumChans2 = 89;
  HDDrive = 'ExtHDD02';
elseif sum(Patient == 'PY05N004') == 8
  NumChans2 = 94;
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY05N005') == 8
  NumChans2 = 110;
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY05N006') == 8
  NumChans2 = 22;
  HDDrive = 'ExtHDD01';
elseif sum(Patient == 'PY05N007') == 8
  NumChans2 = 42;
  HDDrive = 'ExtHDD01';
end

if sum(Patient == 'PY11N003') == 8              % PY11N003
  NumChans2 = 116;
  sfreq = 1000;
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N004') == 8          % PY11N004
  NumChans2 = 112;
  sfreq = 1000;
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N005') == 8          % PY11N005
  NumChans2 = 72;
  sfreq = 1000;
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N006') == 8          % PY11N006
  NumChans2 = 65;
  sfreq = 1000;
  eval(['cd /media/ExtHDD03/',Patient])
%%%%%%% added by Adam %
elseif sum(Patient == 'PY11N007') == 8          % PY11N007
  NumChans2 = length([1:5 7:18 20:33 35:48]);
  sfreq = 1000;
  eval(['cd /media/ExtHDD03/',Patient])
elseif sum(Patient == 'PY11N008') == 8          % PY11N008
  NumChans2 = length([1:113 115:126]);
  sfreq = 1000;
  eval(['cd /media/ExtHDD03/',Patient])
%%%%%%% end of Adam %
elseif sum(Patient == 'PY11N009') == 8          % PY11NOO9
  NumChans2 = 99;
  sfreq = 1000;
  HDDrive = 'ExtHDD03';
%%%%%ADAM %
elseif sum(Patient == 'PY11N010') == 8          % PY11NO10
  NumChans2 = length([1:77 79:85 87:88]);
  sfreq = 1000;
  eval(['cd /media/ExtHDD04/',Patient])
%%%%%END %
elseif sum(Patient == 'PY11N011') == 8          % PY11N011
  NumChans2 = 27;
  sfreq = 200;
  eval(['cd /media/ExtHDD04/',Patient])
  eval(['cd /media/ExtHDD03/',Patient])
%%%%%%% added by Adam %
elseif sum(Patient == 'PY11N012') == 8          % PY11N012
  NumChans2 = length([1:55 57:89 91:110]);
  sfreq = 1000;
  eval(['cd /media/ExtHDD04/',Patient])
elseif sum(Patient == 'PY11N013') == 8          % PY11N013
  NumChans2 = length([1:111 113:119]);
  sfreq = 1000;
  eval(['cd /media/ExtHDD04/',Patient])
%%%%%%% end of Adam % 
elseif sum(Patient == 'PY11N014') == 8          % PY11N014
  NumChans2 = 104;
  sfreq = 1000;
  eval(['cd /media/ExtHDD04/',Patient])
elseif sum(Patient == 'PY11N015') == 8          % PY11N015
  NumChans2 = 27;
  sfreq = 1000;
  HDDrive = 'ExtHDD04';
end

% go to that patient's directory 
eval(['cd /media/',HDDrive,'/',Patient])
%% 02: Power Spectrum Analysis & Settings for FFT

WinLength = 3*sfreq;        % window length for seizure period
FFTlength = sfreq;          % length of FFT
Overlap = 1/3;              % percent overlap
NumFFTs = round(WinLength/(Overlap*FFTlength) - 2);
stimfreq = sfreq*(0:FFTlength/2)/FFTlength;
stimfreq = stimfreq(1:sfreq/2);             % the freq. "bands" account for the '-1' since freq. starts at 0, but MATLAB indices at 1

clear mDataFFTall DataFFTBB FFTcountBB
mDataFFTall = [];
DataFFTBB(NumChans2,FFTlength/2) = zeros;
FFTcountBB(NumChans2) = zeros;
disp(strcat('Running through patient ', Patient));
disp(strcat('NumSZ is: ', num2str(NumSZ)));

if sum(Patient == 'PY04N013') == 8
  clear SZkk NumSZ
  SZkk = [1 12 71 73 76 81];
  NumSz = 6;
end

%% 03: loop through each file for Seizure
%SZkk = [5 6 11 13 14];
%for kk = 1:5
for kk = 1:NumSZ
  disp(strcat('On loop #: ', num2str(kk)));
  
  %pp = SZstartPP(SZkk(kk));
  pp = SZstartPP(kk);
  clear FileName FileDir FilePath Data 
  FilePath = char(InfoTime.filename(pp));
  FileDir  = FilePath(1:5);
  eval(['cd ',FileDir])     % cd into file directory of the hard drive w/ patient
  disp(strcat('FileDir is: ', FileDir));
  
  % Set up Headers
  hdr          = readhdr('eeg.hdr');
  num_channels = hdr.file_fmt.Numb_chans;
  format_file  = hdr.file_fmt.File_format;
  offset       = hdr.file_fmt.Data_offset;

  FileName = FilePath(6:end);
  fid  = fopen([FileName],'rb'); 
  fseek(fid,offset,'bof');
  Data = fread(fid,[num_channels inf],format_file);
  fclose(fid);

  clear NumChans NumSamps a b DataLP
  [NumChans,NumSamps] = size(Data);

  %% 03b: check through each patient and only get the channels we are interested
  if sum(Patient == 'PY04N007') == 8 
    if     sum(FileDir(1:5) == 'rec01') == 5
      Data2 = Data([1:6 8:21 24:27 29:31 33 35:40 42:51 54:56 58:63 65:79 81:87],:);
    end
  end

  if sum(Patient == 'PY04N008') == 8 
    if     sum(FileDir(1:5) == 'rec01') == 5
      Data2 = Data([1:15 18 20:28 33:38],:);
    end
  end

  if sum(Patient == 'PY04N012') == 8 
    if     sum(FileDir(1:4) == 'edf0') == 4
      Data2 = Data([1:34 51:98],:);
    elseif sum(FileDir(1:4) == 'edf1') == 4
      Data2 = Data;
    elseif sum(FileDir(1:4) == 'edf2') == 4
      Data2 = Data([1:82],:);
    end
  end

  if sum(Patient == 'PY04N013') == 8 
    Data2 = Data([1 4:10 12:61 63 65 70 72:89],:);
  end

  if sum(Patient == 'PY04N015') == 8 
    if     sum(FileDir(1:4) == 'edf0') == 4
      Data2 = Data([1:2 5:54 56:83 85:87 89:91 93:95],:);
    elseif sum(FileDir(1:4) == 'edf1') == 4
      Data2 = Data([1:2 5:54 56:83 85:87 89:91 93:95],:);
    end
  end

  if sum(Patient == 'PY05N004') == 8 
    if     sum(FileDir(1:4) == 'edf0') == 4
      Data2 = Data([1:12 14 17:37 39:63 65:82 84:100],:);
    elseif sum(FileDir(1:4) == 'edf1') == 4
      Data2 = Data([1:12 14 17:37 39:81 83:99],:);
    elseif sum(FileDir(1:4) == 'edf2') == 4
      Data2 = Data([1:12 14 17:37 39:81 83:99],:);
    end
  end

  if sum(Patient == 'PY05N005') == 8 
    if     sum(FileDir(1:4) == 'edf0') == 4
      Data2 = Data([1:21 23:38 41:83 89:118],:);
    elseif sum(FileDir(1:4) == 'edf1') == 4
      Data2 = Data([1:21 23:38 41:83 89:118],:);
    end
  end

  if sum(Patient == 'PY05N006') == 8 
    if     sum(FileDir(1:4) == 'edf0') == 4
      Data2 = Data([1:2 5:24],:);
    elseif sum(FileDir(1:4) == 'edf1') == 4
      Data2 = Data([1:2 5:24],:);
    end
  end

  if sum(Patient == 'PY05N007') == 8 
    if     sum(FileDir(1:4) == 'edf0') == 4
      Data2 = Data([1:2 5:7 10:14 17:48],:);
    elseif sum(FileDir(1:4) == 'edf1') == 4
      Data2 = Data([1:2 5:7 10:14 17:48],:);
    end
  end
  
  if sum(Patient == 'PY11N003') == 8 
    Data2 = Data([1:53 55:117],1:NumSamps);
  end

  if sum(Patient == 'PY11N004') == 8 
    Data2 = Data([2 4:15 17:31 33:106 108:117],1:NumSamps);
  end

  if sum(Patient == 'PY11N005') == 8 
    Data2 = Data([1:11 13:33 35:74],1:NumSamps);
  end

  if sum(Patient == 'PY11N006') == 8 
    Data2 = Data([1:16 18:66],1:NumSamps);
  end
  %%%%%%% added by Adam %
  if sum(Patient == 'PY11N007') == 8          % PY11N007
    Data2 = Data([1:5 7:18 20:33 35:48],1:NumSamps);
  end
  if sum(Patient == 'PY11N008') == 8          % PY11N008
    Data2 = Data([1:113 115:126],1:NumSamps);
  end
  %%%%%%% end of Adam %
  if sum(Patient == 'PY11N009') == 8 
    Data2 = Data([1:22 24:30 32:81 85:90 92:105],1:NumSamps);
  end
  %%%%%ADAM %
  if sum(Patient == 'PY11N010') == 8          % PY11NO10
    Data2 = Data([1:77 79:85 87:88], 1:NumSamps);
  end
  %%%%%END %
  if sum(Patient == 'PY11N011') == 8 
    Data2 = Data([6:32],1:NumSamps);
  end
  %%%%%%% added by Adam %
  if sum(Patient == 'PY11N012') == 8          % PY11N012
    Data2 = Data([1:55 57:89 91:110], 1:NumSamps);
  end
  if sum(Patient == 'PY11N013') == 8          % PY11N013
    Data2 = Data([1:111 113:119], 1:NumSamps);
  end
  %%%%%%% end of Adam % 
  if sum(Patient == 'PY11N014') == 8 
    Data2 = Data([1:21 23:28 31:39 42:46 48:52 55:68 70:81 83:111],1:NumSamps);
    %Data2 = Data([1:21 23:40 42:46 48:52 55:68 70:81 83:111],1:NumSamps);
  end
   
  %reset data var
  clear Data
  Data = Data2;
  clear Data2
  
  %% 04: Look at Seizures
  clear SZstart SZstop SZtime
  SZstart = floor((TimesSZ(kk,1) - TimesAll(pp,1))*sfreq);
  SZstop  = ceil((TimesSZ(kk,2)  - TimesAll(pp,1))*sfreq);
  %SZstart = floor((TimesSZ(SZkk(kk),1) - TimesAll(pp,1))*sfreq);
  %SZstop  = ceil((TimesSZ(SZkk(kk),2)  - TimesAll(pp,1))*sfreq);
  SZtime  = SZstop - SZstart;

  clear DataFFT FFTcount mDataFFT
  NumWins = floor(SZtime/sfreq)-WinLength/sfreq;
  DataFFT(NumChans2,NumWins,FFTlength/2) = zeros;
  FFTcount = 1;

  for gg = SZstart:sfreq:SZstop-WinLength
    clear FFTtemp GGindex
    FFTtemp(1:NumFFTs,1:FFTlength) = zeros;
    GGindex = gg:(gg-1)+WinLength;

    for bb = 1:NumChans2
      clear DataTemp mFFTtemp
      DataTemp = squeeze(Data(bb,GGindex));

      for hh = 1:NumFFTs
        clear HHindex DataTemp2
        HHindex = (hh-1)*floor(FFTlength*Overlap)+1:(hh-1)*floor(FFTlength*Overlap) + FFTlength;
        DataTemp2 = detrend(DataTemp(HHindex));
        FFTtemp(hh,1:FFTlength) = abs(fft(DataTemp2)).^2;
        %FFTtemp(hh,1:FFTlength) = abs(fft(hamming(sfreq)'.*DataTemp(HHindex))).^2;
        %FFTtemp(hh,1:FFTlength) = FFTtemp(hh,1:FFTlength)/sum(FFTtemp(hh,1:FFTlength));
      end

      mFFTtemp = mean(FFTtemp); 

      if sum(isnan(mFFTtemp)) == 0
        DataFFT(bb,FFTcount,1:FFTlength/2) = mFFTtemp(1:FFTlength/2);
        DataFFTBB(bb,1:FFTlength/2) = DataFFTBB(bb,1:FFTlength/2)+mFFTtemp(1:FFTlength/2);
        FFTcountBB(bb) = FFTcountBB(bb) + 1;
      end
    end

    FFTcount = FFTcount + 1;  
  end
  
  % update variables
  FFTcount  = FFTcount - 1;
  mDataFFT  = squeeze(mean(DataFFT));
  mDataFFTall = [mDataFFTall mDataFFT'];
  clear Data
  cd ..
end

%% 05: Post Processing
% set seizure fast fourier transform?
clear SZfftBB
for bb = 1:NumChans2
  SZfftBB(bb,:) = DataFFTBB(bb,:)/FFTcountBB(bb);
end 

% set mean fast fourier transform?
clear mSZfft
mSZfft = mean(mDataFFTall');

%% 06: Save Seizure Spectrum & Computed R-spectrum
cd /home/adamli/MATLAB/code_adam
eval(['load IImspec_',Patient])

%calculate average R-spectrum
for hh = 1:NumChans2
  RspecBB(hh,:) = SZfftBB(hh,:)./IIfftBB(hh,:);
end

mRspecBB = mean(RspecBB); % avg. across each channel
eval(['save SZmspec_',Patient,' mSZfft SZfftBB'])
eval(['save Rspec_',Patient,' RspecBB mRspecBB'])