function preProcessChannelCheck()
    load('tempworkspace');
    
    powerMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    powerMatZ = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);

    thisChan = chanList(iChan,:);   % the channel to use in this loop (e.g. 48)
    thisChanStr = chanStr{iChan};
    strStart    = sprintf('\n STEP 5.%d -- Grab %d/%d: %s', iChan, iChan, numChannels, thisChanStr );  strStart(end+1:35)=' '; %buffer length so everything lines up
    fprintf('%s', strStart);   tic;
    
    %%- gete_ms: get the eegWaveV
    % eegwaveform for each event over the duration of time for a certain channel
    eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);

    % notch filter to eliminate 60 Hz noise
    fprintf(' [%.1f sec] --> notch filt\n', toc); tic;
    eegWaveV = buttfilt(eegWaveV,[59.5 60.5],resampledrate,'stop',1); %-filter is overkill: order 1 --> 25 dB drop (removing 5-15dB peak)
    
    if ~ROBUST_SPEC % OPTION 1: perform wavelet spectral analysis
        %%- multiphasevec3: get the phase and power
        % power, phase matrices for events x frequency x duration of time for each channel
        fprintf(' [%.1f sec] --> freq decomp', toc); tic;
        [rawPhase,rawPow] = multiphasevec3(waveletFreqs,eegWaveV,resampledrate,waveletWidth);
        fprintf(' [%.1f sec] --> save', toc);  tic;
        fprintf('\n');

        %%- REMOVE LEADING/TRAILING buffer areas from power, phase,
        %%eegWave, timeVector
        rawPow   = rawPow(:,:,BufferMS+1:end-BufferMS);
        rawPhase = rawPhase(:,:,BufferMS+1:end-BufferMS);
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
        eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
        if length(eegWaveT)<size(eegWaveV,2), % error check on time vs. voltage length
            fprintf('wave time length off'); 
            eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
        end
        % x-axis of time series
        waveT = eegWaveT;

        % temp indicies
        iEv = 1:length(eventTrigger); % # of events
        iT  = 1:size(eegWaveV,2); % # of time points
        iF  = 1:length(waveletFreqs); % # of freqs.
        iChanSave = 1;

        % chan X event X freq X time
        % make power 10*log(power)
        powerMat(iChanSave,iEv,iF,iT) = 10*log10(rawPow);
        phaseMat(iChanSave,iEv,iF,iT) = rawPhase;

    %     for each eegfile stem, z-score each channel and frequency
        fprintf(' [%.1f sec] --> z-score', toc);  tic;
        stemList = unique({eventTrigger.eegfile});
        
        % indices of the powerMat to Z-score wrt
        for iStem=1:length(stemList),
            fprintf('.');
            iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));
            for iF = 1:length(waveletFreqs),
                allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
%                 allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,fixOnToOff)),length(iEvStem)*length(fixOnToOff),1); % normalize wrt fixation period
                mu = mean(allVal); stdev = std(allVal);

                % create the power matrix
                powerMatZ(iChanSave,iEvStem,iF,iT) = (powerMat(iChanSave,iEvStem,iF,iT)-mu)/stdev;

                if sum(isnan(powerMatZ(iChanSave,iEvStem,iF,iT)))>0
                    keyboard;
                end
            end
        end
        % set two paramters from robust spectrotemp to 0 
        tWin = 0;
        freq = 0;
        
        fprintf(' [%.1f sec]', toc); tic;
        clear rawPow rawPhase
        disp('powerMatZ, powerMat and phaseMat are created')
    else % OPTION 2: Perform robust spectrotemporal pursuit instead 
        %%- REMOVE LEADING/TRAILING BUFFER REGION FOR EEGWAVE
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
        eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
        if length(eegWaveT)<size(eegWaveV,2), % error check on time vs. voltage length
            fprintf('wave time length off'); 
            eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
        end
        
        %%- ROBUST SPECTROTERMPORAL PURSUIT PARAMETERS
        % initialize arrays for storing 
        powerMat = zeros(length(eventTrigger), WINDOW/2, length(eegWaveV)/WINDOW); % #EVENTS X #FREQ X #TIME
        powerMatZ = zeros(length(eventTrigger), WINDOW/2, length(eegWaveV)/WINDOW);
        
        % generate robust spectrogram for each event
        tic;
        for iEv=1:length(eventTrigger),
            [xEst,freq,tWin,iter] = specPursuit(eegWaveV(iEv,:), FS, WINDOW, ALPHA);
            size(xEst)
            xEst = 20*log10(abs(xEst)); % convert to power spectrum
            xEst = reshape(xEst, 1, size(xEst,1), size(xEst,2));
            powerMat(iEv,:,:) = xEst;
        end
        fprintf(' [%.1f sec] --> robust spect pursuit', toc);
        
        % reset tWin var to reflect over seconds
        tWin = tWin - timeZero*1000;
        
        %%- Zscore all power for robust spectral pursuit
        fprintf(' [%.1f sec] --> z-score robust spec pursuit', toc); tic;
        for iEvent =1:size(powerMat,1),
            for iF = 1:size(powerMat,2),
                fixedVal = powerMat(iEvent, iF, 1:4); %allVal for particular chan and freq
                mu = mean(fixedVal); stdev = std(fixedVal);

                % create the power matrix
                powerMatZ(iEvent, iF, :) = (powerMat(iEvent, iF, :)-mu)/stdev;

                if sum(isnan(powerMatZ(iEvent,iF,:)))>0
                    keyboard;
                end
            end
        end
        fprintf(' [%.1f sec]', toc);
    end
%     clear powerMat
    
    % normalized to 1 so all can be shifted to fit on a single plot
    eegWaveMeanSub  = eegWaveV-mean(mean(eegWaveV));   %double mean and double max because multiple events from same channel should be normalized together
    eegWaveShift    = eegWaveMeanSub./max(max(abs(eegWaveMeanSub)));
    wavesSft = eegWaveShift;
    
    % create vector of the actual seconds in time axis for the powerMat
    % (since its time binned)...
%     LOWERTIME = 1001;
%     UPPERTIME = 6000;
    OVERLAP = 100;
    WINSIZE = 500;
%     FS = 1000;
%     TIMEZERO = 2000;
    if tWin == 0, % if not set yet
        tWin = (LOWERTIME) :OVERLAP/FS: (UPPERTIME)-WINSIZE/FS;
    end
    timeZero = abs(0-(LOWERTIME))/(OVERLAP/FS);
    
    %%- Plot spectrograms
    powPlot = mean(powerMatZ, 2);
    size(powPlot)
    titleStr = sprintf('mean power: chan %s, %d events', chanStr{iChan}, size(powerMatZ,2));
    powPlot = squeeze(powPlot);
    
    typeSync = 'ProbeWord Locked';
    
    figure
    % first plot eeg trace
    subplot(411)
    eegTitle = sprintf('%s %s: %s : channel %s number(%s) over %d events', subj, typeSync, session, chanStr{iChan}, num2str(chanList(iChan)), size(eegWaveV, 1));
    hold on; plot(eegWaveT, mean(wavesSft, 1));
    set(gca,'tickdir','out','YDir','normal');
    set(gca,'fontsize',figFontAx)
    title(eegTitle, 'fontsize',20)
    xlabel('time (s)')
    ylabel('voltage value (normalized)')
    
    subplot(412)
    hIMg = imagesc(waveT,log10(waveletFreqs),powPlot);
    hold on; colormap(jet);
    hCbar = colorbar('east');
    set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right') 
    set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
    set(gca,'tickdir','out','YDir','normal');
    set(gca,'fontsize',figFontAx)
    title(titleStr, 'fontsize',20)
    
    % remake powerMatZ to the points that we want (before probe on -> 3.5
    % seconds later
    powerMatZ = squeeze(powerMatZ); % only get the powerMatZ time points we want... (1001 - 2000+3500) = -1.0 seconds -> 3.5 seconds
    %% TIME BIN POWERMATZ WITH WINDOWSIZE AND OVERLAP
    addpath('../m_oldAnalysis_anovaANDsinglechannel/');
    WINDOWSIZE = 500; % in milliseconds
    OVERLAP = 100;    % in milliseconds
    powerMatZ = timeBinSpectrogram(powerMatZ, WINDOWSIZE, OVERLAP);
    
    powPlot = squeeze(mean(powerMatZ, 1));
    size(powPlot)
    titleStr = sprintf('mean power: chan %s, %d events', chanStr{iChan}, size(powerMatZ,1));
    
    subplot(413)
    hIMg = imagesc(tWin,log10(waveletFreqs),powPlot);
    hold on; colormap(jet);
    hCbar = colorbar('east');
    set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right')   
    set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
    set(gca,'tickdir','out','YDir','normal');
    set(gca,'fontsize',figFontAx)
    title(titleStr, 'fontsize',20)
    
    %% FREQUENCY BIN WITH FREQUENCY BANDS
    rangeFreqs = reshape([freqBandAr.rangeF], 2, 7)';
    waveletFreqs = waveletFreqs;
    powerMatZ = freqBinSpectrogram(powerMatZ, rangeFreqs, waveletFreqs);
    
    powPlot = squeeze(mean(powerMatZ, 1));
    size(powPlot)
    titleStr = sprintf('mean power: chan %s, %d events', chanStr{iChan}, size(powerMatZ,1));
    
    subplot(414)
    hIMg = imagesc(tWin,1:7,powPlot);
    hold on; colormap(jet);
    hCbar = colorbar('east');
    set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right', 'fontsize',14) 
    set(gca, 'xtick');
    set(gca,'ytick',1:7,...
        'yticklabel',{'delta', 'theta', 'alpha', 'beta', 'low gamma', 'high gamma', 'HFO'}, ...
        'fontsize',20)
    set(gca,'tickdir','out','YDir','normal');
    set(gca,'fontsize',figFontAx)
    title(titleStr, 'fontsize',20)
    
    figureDir = strcat('../Figures/', subj, '_SpectCheck/');
    figureFile = strcat(figureDir, chanStr{iChan});
    if ~exist(figureDir)
        mkdir(figureDir)
    end
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    pos = [20.9722    2.2917   13.8750   15.5694];
    fig.PaperPosition = pos;

    %%- Save Image
    print(figureFile, '-dpng', '-r0')
    
    pause(0.5)
	close all
end