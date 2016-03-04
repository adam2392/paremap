% RAW DATA SECTION 
% rawdata.waveletFreqs = waveletFreqs;    % store the wavelet freqs.
% rawdata.eegWaveV = eegWaveV;            % store the voltage time series per event
% rawdata.resampledrate = resampledrate;  % store resampling rate
% rawdata.waveletWidth = waveletWidth;    % width of wavelets
% rawdata.uniqueTrigType = uniqueTrigType;% the unique trigger types (e.g. brick, clock)
% rawdata.trigType = trigType;            % trigger per event
% rawdata.waveT = waveT;                  % the time points per eeg voltage
% rawdata.chanNum = channel_num;
% rawdata.chanStr = chanStr;
% rawdatafreqBandYtick = freqBandYticks;
% rawdata.freqBandYlabel = freqBandYtickLabels;

data = load('48_MST2-global_rawdata.mat');
data = data.rawdata;

%- For raw data
eegWaveV = data.eegWaveV;            % get the voltage time series for all events
trigType = data.trigType;            % get the trigger vector
fs = data.resampledrate;
waveletFreqs = data.waveletFreqs;
waveT = data.waveT;

% robust spect parameters
alpha = 10;
window = 41;
% temp indicies for saving into 4D power Matrix
iT  = 1:size(eegWaveV,2); % # of time points
iF  = 1:length(trigType); % # of freqs.
numChanPrealloc = 1;

powerMat  = zeros(length(trigType), length(waveletFreqs), length(waveT));

%%- loop through all events
for iTrig = 1:length(trigType)
    trigger = trigType(iTrig);      % get current trigger
    eegWave = eegWaveV(iTrig,:);    % eegWave time series for this event
    
    % robust spectrograms
    % Spectrotemporal Pursuit Algorithm
    % Input:
    % data: time series (1D vector)
    % fs: sampling rate (in Hz)
    % window: size of window (equal to number of frequency bands) 
    % alpha: model parameter
    tic;
    [powerEst, freq, tWin, iter] = specPursuit(eegWaveV, fs, window, alpha);
    toc;
    
    
    % chan X event X freq X time
    % make power 10*log(power)
%     powerMat(iTrig,iF,iT) = 10*log10(powerEst);
    figure
    imagesc(waveT,freq,10*log10(abs(powerEst)));axis xy;colorbar;
end


%  chan X event X freq X time
% make power 10*log(power)
% powerMat(iChanSave,iEv,iF,iT) = 10*log10(rawPow);
% phaseMat(iChanSave,iEv,iF,iT) = rawPhase;
% 
% % for each eegfile stem, z-score each channel and frequency
% fprintf(' [%.1f sec] --> z-score', toc);  tic;
% stemList = unique({eventTrigger.eegfile});
% for iStem=1:length(stemList),
%     fprintf('.');
%     iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));
%     for iF = 1:length(waveletFreqs),
%         allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
%         mu = mean(allVal); stdev = std(allVal);
% 
%         % create the power matrix
%         powerMatZ(iChanSave,iEvStem,iF,iT) = (powerMat(iChanSave,iEvStem,iF,iT)-mu)/stdev;
% 
%         if sum(isnan(powerMatZ(iChanSave,iEvStem,iF,iT)))>0
%             keyboard;
%         end
%     end
% end