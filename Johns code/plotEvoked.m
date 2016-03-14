function [statsStruct, responseStruct, mainFigHandle] = plotEvoked(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Function PlotEvoked(data)
%
% plot the evoked response, broken down by trigger type
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%-- createobject that stiches together all data for saving and/or analysis
% %- figure output parameters
% data.subj           = subj;
% data.session        = session;
% data.figFontAx      = figFontAx;
% data.FIG_OFFSET     = FIG_OFFSET;
% data.HIDE_FIGURES   = HIDE_FIGURES;
%
% %- extraction parameters (filter and wavelet info)
% data.waveletFreqs   = waveletFreqs;
% data.waveletWidth   = waveletWidth;
% data.preFiltFreq    = preFiltFreq;
% data.preFiltType    = preFiltType;
% data.preFiltOrder   = preFiltOrder;
% data.freqBandAr     = freqBandAr;
% data.freqBandYticks = freqBandYticks;
% data.freqBandYtickLabels = freqBandYtickLabels;
% data.chanList       = chanList;
% data.chanStr        = chanStr;
% data.thisChan       = thisChan; %applies to last channel processed, important when processing is sequential
% data.thisChanStr    = thisChanStr;

% %- extracted EEG: evoked waveforms
% data.waveT          = waveT;      %- time
% data.wavesRaw       = wavesRaw;   %- raw wave
% data.wavesMsb       = wavesMsb;   %- mean subtracted
% data.wavesSft       = wavesSft;   %- mean shifted so multiple traces do not overlap
% data.wavesSftG      = wavesSftG;  %- mean shifted and scaled to global max to represent absolute voltage
% data.instPow        = instPow;    %- instantaneous power
% data.instPowSft     = instPowSft; %- power mean shifted for non-overlapping traces
%
% %- extracted EEG: power and phase
% data.powerMat       = powerMat;
% data.phaseMat       = phaseMat;
% data.powerMatZ      = powerMatZ;  %- each wavelet frequency z-scored against all extracted data
%
% %- meta events variables
% data.eventsMeta     = eventsMeta;
% data.metaYval       = metaYval;
% data.metaZeroMS     = metaZeroMS;
% data.metaYstr       = metaYstr;
% data.metaYtransition = metaYtransition;
%
% %- trigger event variables
% data.THIS_TRIGGER   = THIS_TRIGGER;
% data.eventsTrig     = eventsTrig;
% data.trigType       = trigType;
% data.trigZeroMS     = trigZeroMS;
% data.trigTypeStr       = trigTypeStr;
% data.eventsTriggerXlim = eventsTriggerXlim;
% data.eventsAveWindowMS = eventsAveWindowMS;

%%%%% pull out individual variables from global/passed structure %%%%%%

subj            = data.subj;
session         = data.session;
figFontAx       = data.figFontAx;
FIG_OFFSET      = data.FIG_OFFSET;
HIDE_FIGURES    = data.HIDE_FIGURES;

chanStr         = data.chanStr;
thisChanStr     = data.thisChanStr;  %applies to last channel processed, important when processing is sequential
waveletFreqs    = data.waveletFreqs;
freqBandAr      = data.freqBandAr;
freqBandYticks  = data.freqBandYticks;
freqBandYtickLabels = data.freqBandYtickLabels;

THIS_TRIGGER    = data.THIS_TRIGGER;
trigType        = data.trigType;
eventsTrig      = data.eventsTrig;
trigTypeStr     = data.trigTypeStr;
trigZeroMS      = data.trigZeroMS;
metaYstr        = data.metaYstr;
eventsTriggerXlim = data.eventsTriggerXlim;

waveT           = data.waveT;
wavesSft        = data.wavesSft;
wavesSftG       = data.wavesSftG;
powerMatZ       = data.powerMatZ;
powerMat        = data.powerMat;
phaseMat        = data.phaseMat;
    

numChan = size(powerMatZ,1);
if numChan==1 && length(chanStr)~=1, %- numChan~=length(chanStr) means processed sequentially... make sure title refects correct channel
    %keyboard
    chanStr{1} = thisChanStr;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-----------------------   Evoked responses: comparing different conditions      ------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  raw         mean-sub    +shift&scale    +shft&scaleGlobal    rawInstPow    shift&scale
% wavesRaw     wavesMsb      wavesSft      wavesSftG            instPow      instPowSft
%wavePlotSelect = wavesRaw;  waveYlabel = 'Evoked Response (uV)';
wavePlotSelect = wavesSft;  waveYlabel = 'Evoked Response (norm)';


% identify and trim high-count events
uniqueTrigType = unique(trigType);
numUniqueTrig  = length(uniqueTrigType);
for thisTrig=1:numUniqueTrig,  numTrig(thisTrig) = length(find(trigType==uniqueTrigType(thisTrig))); end
numTrigSort = sort(numTrig,'descend');
iTrigTrim = -1;
if length(numTrigSort)>1 & numTrigSort(1)>1.5*numTrigSort(2),
    iKeep=unique(round(1:numTrigSort(1)/numTrigSort(2):numTrigSort(1)));
    iTrigTrim = find(numTrig==max(numTrig));
    %fprintf('\n note: category %d has way more events than the others (%d events), will trim to %d events',iTrigTrim,numTrigSort(1),length(iKeep));
end


if numUniqueTrig>1,
    for chanNum = 1:numChan,  % loop over channels --- should always be just 1 when calling this function
        
        %-- Select Figure and open/hide --%
        thisFigNum = 1000+chanNum+FIG_OFFSET;
        if ~ishghandle(thisFigNum), figure(thisFigNum);
        else                    set(0,'CurrentFigure', thisFigNum); end
        if HIDE_FIGURES==1,     set(thisFigNum,'visible','off');    end
        
        clf
        set(gcf,'color','w')
        cTrigAr = 'rkbgmcy';
        
        
        for thisTrig=1:numUniqueTrig,
            iTrig = find(trigType==uniqueTrigType(thisTrig));
            if thisTrig==iTrigTrim, iTrig=iTrig(iKeep); end    %... trim excess events
            %iTrig = iTrig(1:3);                               %... look at single events by limiting iTrig here...
            
            cTrig = cTrigAr(thisTrig);
            if strcmp(THIS_TRIGGER(1:6),'sample') | strcmp(THIS_TRIGGER(1:4),'test') | strcmp(THIS_TRIGGER,'asterisk') ,
                clnResultStr = regexprep(eventsTrig(iTrig(1)).resultStr,'_f','');
                clnResultStr = regexprep(clnResultStr,'*','');
                if     strcmp(clnResultStr,'n/a') & length(find([eventsTrig(iTrig).responseCorrect]==0))==length(iTrig),  clnResultStr = 'E';      %-all errors (incorrect)
                elseif strcmp(clnResultStr,'n/a') & length(find([eventsTrig(iTrig).responseCorrect]==1))==length(iTrig),  clnResultStr = 'C';  end %-all correct
                %thisTrigTypeStr = sprintf('%s (%s) [%d ev]',eventsTrig(iTrig(1)).asteriskStr, clnResultStr, length(iTrig));
                thisTrigTypeStr = sprintf('%s (%s)',eventsTrig(iTrig(1)).asteriskStr, clnResultStr);  %-add event counts below
            else
                thisTrigTypeStr = sprintf('%s %d', trigTypeStr, uniqueTrigType(thisTrig));
            end
            
            
            
            %%----------------- plot event types--------------------------------------------------------
            subplot(2+numUniqueTrig,1,1)
            axTask = gca;
            
            events   = eventsTrig(iTrig);
            yValues  = thisTrig; % increment sequentially... don't use value of trigType, which can jump
            tOffset  = trigZeroMS(iTrig);
            
            
            %-- event start and stop time
            eTimeOnS  = ([events.mstime]-tOffset)/1000;
            eTimeOffS = ([events.mstime]+[events.msDuration]-tOffset)/1000;
            hP = plot(eTimeOnS,yValues,'b.','MarkerSize',10); hold on;               % all events get a blue + to indicate start
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',4); % all events get a line indicating duration
            %-- next (meta) event start time
            eTimeEnxtS = ([events.msMetaEventNext]-tOffset)/1000;   % time that next meta event starts
            hP = plot(eTimeEnxtS,yValues-.1,'k>','MarkerSize',6);
            
            %-- other event markers that can be used with trigger event or meta event
            % time cross is presented
            eTimeOnS  = ([events.msCrossStart]-tOffset)/1000;
            eTimeOffS = ([events.msCrossStop] -tOffset)/1000;
            hP = plot(eTimeOnS,yValues-.15,'k+','MarkerSize',10);  % time cross is presented
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',2); % all events get a line indicating duration
            % time sample/test text is presented
            eTimeOnS  = ([events.msTextStart]-tOffset)/1000;
            eTimeOffS = ([events.msTextStop] -tOffset)/1000;
            hP = plot(eTimeOnS,yValues-.15,'ks','MarkerSize',6);
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',2); % all events get a line indicating duration
            % time asterisk is presented
            eTimeOnS  = ([events.msAsteriskStart]-tOffset)/1000;
            eTimeOffS = ([events.msAsteriskStop] -tOffset)/1000;
            hP = plot(eTimeOnS,yValues-.15,'kp','MarkerSize',12);
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',2); % all events get a line indicating duration
            
            %-- label trigger type
            cleanTrigTypeStr = regexprep(thisTrigTypeStr,' \[.*\]', ''); %- remove bracketed expression, just for top panel
            hText = text(min(eventsTriggerXlim)-.05*abs(min(eventsTriggerXlim)), mean(yValues), cleanTrigTypeStr);
            set(hText,'HorizontalAlignment','right','FontSize',figFontAx);
            
            set(gca,'ydir','reverse','XAxisLocation','bottom','tickdir','out','fontsize',figFontAx)
            set(gca,'ylim',[1-.2 length(uniqueTrigType)+.2],'ytick',[1:length(uniqueTrigType)],'xlim',[eventsTriggerXlim], 'ydir','reverse')
            axis off
            box off
            
            
            %-- figure title
            titleStr = sprintf('%s : %s : channel %s', subj, session, chanStr{chanNum});
            titleStr(find(titleStr=='_'))=' ';
            title(titleStr, 'fontsize',figFontAx+7)
            
            
            
            %%----------------- plot evoked responses--------------------------------------------------------
            subplot(2+numUniqueTrig,1,2)
            axWave = gca;
            
            %-compute the evoked response for each trigger type
            wavePlot = wavePlotSelect(:,iTrig,:);
            clear evokedRespMu evokedRespSEM;
            evokedRespMu(1:size(wavePlot,3))  = mean(wavePlot(chanNum,:,:),2);
            evokedRespSEM(1:size(wavePlot,3)) = std(wavePlot(chanNum,:,:),0,2)./sqrt(size(wavePlot,2));
            
            hEvoked = plot(waveT, evokedRespMu, 'k-'); hold on
            set(hEvoked,'Color',cTrig,'LineWidth',2)
            
            ii    = length(evokedRespMu):-1:1;
            fillX = [waveT waveT(ii)];
            fillY = [evokedRespMu+evokedRespSEM evokedRespMu(ii)-evokedRespSEM(ii)];
            hFill = fill(fillX,fillY,cTrig);
            set(hFill,'LineStyle','none','FaceAlpha',0.5);
            
            %%- Legend label
            hLegIn(thisTrig)  = hFill;
            hLegStr{thisTrig} = sprintf('%s [%d ev] .',thisTrigTypeStr,size(wavePlot,2));  %
            
            axis tight
            box off
            set(gca,'tickdir','out','fontsize',figFontAx-3)
            %set(gca,'ydir','reverse','XAxisLocation','bottom'); %only necessary when multiple channels plotted at once... single channel doesnt need reverse
            
            ylabel(waveYlabel,'fontsize',figFontAx)
            
            if numChan==1,
                axis tight;
                ax = axis;
                set(gca,'ylim',[-1 1]*max(abs(ax(3:4)))); %make y-limits symetric
                %    set(gca,'ylim',waveYlim);
            end
            
            %%- following assumes multiple channels... but this function designed to process only 1 channel at a time
            %skipChan = round(size(wavePlot,1)/10); if skipChan<1, skipChan=1; end
            %set(gca,'ytick',[0:skipChan:size(wavePlot,1)-1]*2,'yticklabel',[1:skipChan:size(wavePlot,1)])
            %set(gca,'ytick',[0:skipChan:size(wavePlot,1)-1]*2,'yticklabel', chanStr([1:skipChan:size(wavePlot,1)])); %
            
            
            %%----------------- plot the spectrograms--------------------------------------------------------
            subplot(2+numUniqueTrig,1,2+thisTrig)
            axSpec(thisTrig)=gca;
            
            % use iTrig to select the appropriate events
            thisPowMat = powerMatZ(chanNum,iTrig,:,:); if thisTrig==1, cAx = [-.8 .8]; if length(iTrig)<4,cAx=cAx*3; elseif length(iTrig)>100, cAx=cAx/2; end;   end;   % z-scored power... set cAx for first panel and keep same for all others
            %thisPowMat = powerMat(chanNum,iTrig,:,:); if thisTrig==1, cAx = [max(max(mean(powerMat,2)))-60  max(max(mean(powerMat,2)))];  end;  % power
            
          
            % reduce the dimensionality of the matrix for plotting
            powPlot = mean(thisPowMat,2); titleStr = sprintf('mean power: chan %s, %d events', chanStr{chanNum}, size(thisPowMat,2));
            powPlot = squeeze(powPlot); % squeeze out the singleton dimension
            hImg    = imagesc(waveT,log10(waveletFreqs),powPlot); hold on;  colormap(jet);
            
            %- a few different options for y-tick labels
            %set(gca,'ytick',log10(waveletFreqs(1:skipFreqLabel:end)),'yticklabel',freqBandLabels(1:skipFreqLabel:end))
            set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
            %set(gca,'YTick',[])
            
            set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
            set(gca,'fontsize',figFontAx-3)
            set(gca,'XTick',[],'Box','off');
            ylabel( sprintf('%s\n[%d ev]', thisTrigTypeStr, size(thisPowMat,2)),'fontsize',figFontAx)
            
            
            caxis(cAx)
            hCbar = colorbar('east');
            set(hCbar,'ycolor',[1 1 1]*.1, 'fontsize', figFontAx-3, 'YAxisLocation', 'right')
            %ylabel(hCbar, 'morl 10*log10 power', 'fontsize',figFontAx, 'color',[1 1 1]*.1)
            box off
            
            % label the ranges for each band
            xlim=get(gca, 'xlim');
            for iFB=1:length(freqBandAr),
                hB = plot(xlim, log10([1 1]*freqBandAr(iFB).rangeF(1)), '--', 'linewidth', 2, 'color', [1 1 1]*.7);
                hB = plot(xlim, log10([1 1]*freqBandAr(iFB).rangeF(2)), '--', 'linewidth', 2, 'color', [1 1 1]*.7);
                text(xlim(1)+range(xlim)*.02, log10(freqBandAr(iFB).centerF) , freqBandAr(iFB).name, 'FontSize', 16, 'horizontalalignment','left');
            end
            
            
            %- create structure that saves the average evoked and spectrogram
            responseStruct.subj          = subj;
            responseStruct.session       = session;
            responseStruct.thisChanStr   = thisChanStr;
            responseStruct.chanSubjStr   = sprintf('%s [%s]', thisChanStr, subj);
            responseStruct.time = waveT;
            responseStruct.freq = waveletFreqs;
            responseStruct.thisTrigStr{thisTrig}  = hLegStr{thisTrig};
            responseStruct.evokedRespMu{thisTrig} = evokedRespMu;
            responseStruct.averagePower{thisTrig} = powPlot;
            thisITFCplot = squeeze(abs(mean(exp(1i*phaseMat(chanNum,iTrig,:,:)),2)));  %- inter-trial phase coherence
            responseStruct.ITFC{thisTrig}         = thisITFCplot; 
        end
        
        
        
        %%%%---- Compute sliding Window difference between evoked responses -----%%%%
        
        %%- first create a set of unique strings for each trigger type (very similar to loop above, but slightly different, so copy here.... -%%
        clear trigTypeStrCells; trigMembers = [];
        for thisTrig=1:numUniqueTrig,
            iTrig = find(trigType==uniqueTrigType(thisTrig));
            if thisTrig==iTrigTrim, iTrig=iTrig(iKeep); end    %... trim excess events
            
            if strcmp(THIS_TRIGGER(1:6),'sample') | strcmp(THIS_TRIGGER(1:4),'test') | strcmp(THIS_TRIGGER,'asterisk') ,
                clnResultStr = regexprep(eventsTrig(iTrig(1)).resultStr,'_f','');
                clnResultStr = regexprep(clnResultStr,'*','');
                if     strcmp(clnResultStr,'n/a') & length(find([eventsTrig(iTrig).responseCorrect]==0))==length(iTrig),  clnResultStr = 'E';      %-all errors (incorrect)
                elseif strcmp(clnResultStr,'n/a') & length(find([eventsTrig(iTrig).responseCorrect]==1))==length(iTrig),  clnResultStr = 'C';  end %-all correct
                thisTrigTypeStr = sprintf('%s (%s) [%d ev]',eventsTrig(iTrig(1)).asteriskStr, clnResultStr, length(iTrig));
            else
                thisTrigTypeStr = sprintf('%s %d', trigTypeStr, uniqueTrigType(thisTrig));
            end
            
            trigMembers = [trigMembers iTrig];
            for iTr=iTrig, trigTypeStrCells{iTr} = thisTrigTypeStr;  end
        end
        trigMembers = unique(trigMembers);
        
        
        %%- now loop over windowed times, compute ave of each trial's trace, then test for differences between groups
        winWidth = 50;
        winDT    = 10;
        slidingWindow = [[eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'-winWidth/2 [eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'+winWidth/2];  %-new way... centered for both window sizes
        
        for iTimeWindow=1:size(slidingWindow,1),
            
            thisTimeWindow = slidingWindow(iTimeWindow,1:2);
            iTime          = find(waveT*1000>=thisTimeWindow(1) & waveT*1000<=thisTimeWindow(2));
            
            
            %-compute time-averaged evoked response
            timeAve = mean(wavePlotSelect(chanNum,:,iTime),3);
            
            
            %- statistical test comparing distributions
            [pVal, table, stats] = anova1(timeAve(trigMembers),trigTypeStrCells(trigMembers),'off');
            SSvals = cell2mat([table(2,2) table(4,2)]);
            varExpl = 100*SSvals(1)/SSvals(2);
            [comparison, means, hFig, groupNames] = multcompare(stats,'display','off'); %- use this to define sigPairs below... 
            
            %- create binary vector & t-test matrix indicating significantly different pairs...
            labelPairs = {};
            sigPairs   = zeros(size(comparison,1),1);   clear ttPvalPairs tTstat;
            for iComp=1:size(comparison,1),
                if ~(comparison(iComp,3)<=0 & comparison(iComp,5)>=0),
                    sigPairs(iComp)=1;
                end
                labelPairs{iComp} = sprintf('%s  vs  %s', groupNames{comparison(iComp,1)}, groupNames{comparison(iComp,2)});
                
                
                %- actually do t-test to get pairwise p-values and t-statistics
                iTrigA = find(trigType==uniqueTrigType(comparison(iComp,1)));
                iTrigB = find(trigType==uniqueTrigType(comparison(iComp,2)));
                
                [ttH,ttPval,ttCI,ttSTATS] = ttest2(timeAve(iTrigA),timeAve(iTrigB));
                ttPvalPairs(iComp) = ttPval;
                ttTstat(iComp)     = ttSTATS.tstat;
                
            end
            sigEpocPairs(iTimeWindow,1:size(comparison,1))     = sigPairs;
            ttPvalEpocPairs(iTimeWindow,1:size(comparison,1))  = ttPvalPairs;
            ttTstatEpocPairs(iTimeWindow,1:size(comparison,1)) = ttTstat;
            
            pValEpocs(iTimeWindow)    = pVal;
            varExplEpocs(iTimeWindow) = varExpl;
            labelEpocs{iTimeWindow}   = sprintf('%d to %d ms', thisTimeWindow);
            
        end %-time window loop
        
        pValEpocAdj = pAdjust(pValEpocs); %- adjust for multiple comparisons (BH correction)
        
        set(gcf,'CurrentAxes',axWave);
        yVal  = min(get(axWave,'Ylim'));
        pSigX = mean(slidingWindow,2)/1000;
        
        
        %%
        pThresh = 0.05;
        
        
        %%- grab the pVals and threshold on the fly
        sigEpocs = pValEpocs;  % metaPval(chan,band,time);
        sigEpocs(find(sigEpocs>=pThresh))=nan;
        sigEpocs(find(~isnan(sigEpocs)))=1;
        pSigY = sigEpocs*yVal;
        hAst = plot(pSigX,pSigY,'k*','MarkerSize',8); %%- plot the asterisks
        
        
        %%- grab the pVals and threshold on the fly
        sigEpocs = pValEpocAdj;  % metaPval(chan,band,time);
        sigEpocs(find(sigEpocs>=pThresh))=nan;
        sigEpocs(find(~isnan(sigEpocs)))=1;
        pSigY = sigEpocs*yVal;
        hAst = plot(pSigX,pSigY,'r*','MarkerSize',8); %%- plot the asterisks
        
        
        %%- plot line indicating window width (don't bother... window sooo small)
        %iWin=min(find(mean(slidingWindow,2)/1000>eventsTriggerXlim(2)-0.250));
        %hL = line(slidingWindow(iWin,:)/1000,[1 1]*yVal);
        %set(hL,'Color','r','LineWidth',6)
        %hP = plot(mean(slidingWindow(iWin,:))/1000,[1 1]*yVal);
        %set(hP,'Color','r','Marker','+','MarkerSize',8)
        hT = text(eventsTriggerXlim(1),yVal*1.2,sprintf('win %d ms\nstep %d ms\nblk p<%.03f\nred pAdj',winWidth,winDT,pThresh));
        set(hT,'HorizontalAlignment','Right')
            
        
        %%- Save the stats test results into output structure
        statsStruct.subj          = subj;
        statsStruct.session       = session;
        
        statsStruct.thisChanStr   = thisChanStr;
        statsStruct.chanSubjStr   = sprintf('%s [%s]', thisChanStr, subj);
        %statsStruct.labelBands    = {'Evoked Response'};
        statsStruct.labelEpocs    = labelEpocs;
        statsStruct.labelPairs    = labelPairs;
        
        statsStruct.pValEpocs     = pValEpocs;
        statsStruct.varExplEpocs  = varExplEpocs;
        statsStruct.sigEpocPairs  = sigEpocPairs;
        statsStruct.ttPvalEpocPairs  = ttPvalEpocPairs;
        statsStruct.ttTstatEpocPairs = ttTstatEpocPairs;
        
        statsStruct.eventsAveWindowMS = slidingWindow;
        statsStruct.eventsTriggerXlim = eventsTriggerXlim;
        
        
        
        
        % tweak some properties of the evokedResponse plot
        set(gcf,'CurrentAxes',axWave);
        %axes(axWave); %make hidden figure visible
        
        hLeg = legend(hLegIn,hLegStr);
        %set(hLeg,'fontsize',figFontAx-2)
        xlabel(sprintf('Time from %s (s)',metaYstr),'fontsize',figFontAx+5)
        grid on
        
        % link and resize the axes
        linkaxes([axTask, axWave, axSpec],'x')
        axStimPos = get(axWave,   'position');  %position: [left, bottom, width, height]
        set(axTask, 'position', [axStimPos(1) .90 axStimPos(3) .05], 'tickdir','out')
        set(axWave, 'position', [axStimPos(1) .68 axStimPos(3) .20], 'tickdir','out')
        ht = .62/length(axSpec);
        gap = 0.01;
        step = gap;
        for iSpec=length(axSpec):-1:1,
            set(axSpec(iSpec), 'position', [axStimPos(1) step axStimPos(3) ht-gap], 'tickdir','out');
            step=step+ht;
        end
        
        
    end %loop over channels (each channel gets its own fig)
end %if more than one trigger type...


mainFigHandle = gcf;


