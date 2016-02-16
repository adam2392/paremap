
function [statsStruct, figList] = plotTimeFreqAve(data, slidingWindow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Function plotTimeFreqAve(data)
%
% compute the mean power for different conditions and time-frequency windows
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
% data.thisChanStr    = thisChanStr;  %applies to last channel processed, important when processing is sequential
%
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
thisChanStr     = data.thisChanStr;
waveletFreqs    = data.waveletFreqs;
freqBandAr      = data.freqBandAr;
%freqBandYticks  = data.freqBandYticks;
%freqBandYtickLabels = data.freqBandYtickLabels;

THIS_TRIGGER    = data.THIS_TRIGGER;
% trigType        = data.trigType;
eventsTrig      = data.eventsTrig;
trigTypeStr     = data.trigTypeStr;
trigZeroMS      = data.trigZeroMS;
metaYstr        = data.metaYstr;
eventsTriggerXlim = data.eventsTriggerXlim;
eventsAveWindowMS = data.eventsAveWindowMS;

waveT           = data.waveT;
%wavesSftG       = data.wavesSftG;
powerMatZ       = data.powerMatZ;
powerMat        = data.powerMat;


numChan = size(powerMatZ,1);
if numChan==1 && length(chanStr)~=1,  %- numChan~=length(chanStr) means processed sequentially... make sure title refects correct channel
    %keyboard
    chanStr{1} = thisChanStr;
end
if numChan>1, fprintf('\nWARNING: output statsStruct not setup for multiple channels at once'); keyboard; end


PLOT_SEPARATE_MULTICOMPARE = 0;  %-


SLIDING_PLOT=0;
if ~isempty(slidingWindow), SLIDING_PLOT=1; eventsAveWindowMS = slidingWindow; end


figList  = [];
statsStr = '';
sigEpocs  = zeros(length(freqBandAr),size(eventsAveWindowMS,1));
pValEpocs = zeros(length(freqBandAr),size(eventsAveWindowMS,1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-----------------------   Average Time-Frequency Computation, Plot, Comparison  ------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


% loop over channels --- should never be more than one, but just in case....
for chanNum = 1:numChan,  
    
    %-- Select Figure and open/hide --%
    if SLIDING_PLOT==0,
        thisFigNum = 2000+chanNum+FIG_OFFSET;
        if ~ishghandle(thisFigNum), figure(thisFigNum);
        else                    set(0,'CurrentFigure', thisFigNum); end
        if HIDE_FIGURES==1,     set(thisFigNum,'visible','off');    end
        
        clf
        set(gcf,'color','w')
        cTrigAr = 'rkbgmcy';
        
        figList = gcf;
    end
    
    %-- figure title
    titleStr = sprintf('%s : %s : channel %s', subj, session, chanStr{chanNum});
    titleStr(find(titleStr=='_'))=' ';
    
    
    %- Loop over frequency bands
    for iFreqBand = 1:length(freqBandAr),
        
        thisBand = freqBandAr(length(freqBandAr)-iFreqBand+1);
        iWavelet = find(waveletFreqs>=thisBand.rangeF(1) & waveletFreqs<=thisBand.rangeF(2));
        %fprintf('\n %s, %.0f to %.0f Hz: %d frequency bins included', thisBand.name, thisBand.rangeF, length(iFreq));
        
        
        %- Loop over time windows
        for iTimeWindow=1:size(eventsAveWindowMS,1),
            
            thisTimeWindow = eventsAveWindowMS(iTimeWindow,1:2);
            iTime          = find(waveT*1000>=thisTimeWindow(1) & waveT*1000<=thisTimeWindow(2));
            
            %- compute the timeFreqAverage metric (z-scored power)
            timeFreqAve    = mean(mean(powerMatZ(chanNum,:,iWavelet,iTime),4),3); %mean over samples, then frequencies
            
            
            %- create the subplot
            if SLIDING_PLOT==0,
                iPlot = iTimeWindow + (iFreqBand-1)*size(eventsAveWindowMS,1);
                subplot(length(freqBandAr),size(eventsAveWindowMS,1),iPlot);
                
                %- Loop over trigger types
                clear trigTypeStrCells ocAll; trigMembers = [];
                for thisTrig=1:numUniqueTrig,
                    iTrig = find(trigType==uniqueTrigType(thisTrig));
                    if thisTrig==iTrigTrim, iTrig=iTrig(iKeep); end    %... trim excess events
                    cTrig = cTrigAr(thisTrig);
                    
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
                    
                    
                    %- histogram of observed values
                    if numTrigSort(2)>100,
                        dZ = .1;
                    else
                        dZ = 0.2;
                    end
                    dZ = dZ/log10(thisBand.centerF);  %attempt to balance against frequency differences
                    zPowEdges = [-inf -5.0:dZ:5.0 inf];
                    edgeClean = [zPowEdges(2)-dZ zPowEdges(2:end-1) zPowEdges(end-1)+dZ];
                    
                    oc = histc(timeFreqAve(iTrig),zPowEdges);
                    ocNorm = oc./max(oc); %if normalizing by sum must use same bin width for each trace
                    
                    iStep = []; for i=1:length(oc), iStep = [iStep i i]; end % use this index to create outline of bar plot
                    
                    histX = edgeClean(iStep);
                    histY = [0 ocNorm(iStep(1:end-2)) 0];
                    alphaPlot = 0.5;
                    
                    %- add these three lines to separate the plots vertically
                    histY = histY./numUniqueTrig/1.05 + ((numUniqueTrig-thisTrig)/numUniqueTrig) ; %separate plots
                    hLine = line(histX,zeros(size(histX))+((thisTrig-1)/numUniqueTrig)); set(hLine,'color','k');
                    alphaPlot = 1.0;
                    
                    hFill = fill(histX, histY,cTrig); hold on;
                    set(hFill,'LineStyle','none','FaceAlpha',alphaPlot);
                    
                    %- Legend label
                    hLegIn(thisTrig)  = hFill;
                    hLegStr{thisTrig} = thisTrigTypeStr;
                    
                    ocAll(thisTrig,1:length(oc)) = oc;
                    
                end %-trigger type loop
                
                %-set x-limit
                zNonZero = edgeClean(find(sum(ocAll,1)>0));
                xLimVal  = [min(zNonZero)-1*dZ max(zNonZero)+2*dZ];  %asymetric because using edges
                if length(xLimVal)<2, xLimVal = [min(edgeClean) max(edgeClean)]; end
                set(gca,'xlim', xLimVal);
                
                
                if iFreqBand==1 & iTimeWindow==round(size(eventsAveWindowMS,1)/2), title(titleStr, 'fontsize',figFontAx);   end
                if iFreqBand==1 & iTimeWindow==size(eventsAveWindowMS,1),   legend(hLegIn, hLegStr, 'Location', 'EastOutside'); end
                if iTimeWindow==1,                  ylabel(thisBand.name, 'fontsize', figFontAx-1);                          end
                if iFreqBand==length(freqBandAr),   xlabel(sprintf('%d to %d ms',thisTimeWindow), 'fontsize', figFontAx-3);     end
                
                set(gca,'tickdir','out');
                set(gca,'fontsize',figFontAx-8)
                set(gca,'YTick',[],'Box','off');
                %set(gca,'Ylim',[0 1.05], 'Xlim', [-1 1]*max(abs(zPowBins)))
                %axis tight
                ax = axis;
                axis([ax(1) ax(2) 0 1.28]);
            
            else
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
            end %- SLIDING_PLOT
           
            
            %- statistical test comparing distributions
            trigMembers = unique(trigMembers);
            [pVal, table, stats] = anova1(timeFreqAve(trigMembers),trigTypeStrCells(trigMembers),'off');
            SSvals = cell2mat([table(2,2) table(4,2)]);
            varExpl = 100*SSvals(1)/SSvals(2);

            if stats.df==0,
                keyboard;
            end
            [comparison, means, hFig, groupNames] = multcompare(stats,'display','off');
            pStr = 'n.s.'; fS = 12;
            if pVal<0.1, fS=20; pStr = ''''; if pVal<0.05,  pStr = '*'; if pVal<0.01, pStr = '**'; if pVal<0.001, pStr = '***'; end; end; end; end
            
            
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
                    
                [ttH,ttPval,ttCI,ttSTATS] = ttest2(timeFreqAve(iTrigA),timeFreqAve(iTrigB));
                ttPvalPairs(iComp) = ttPval;
                ttTstat(iComp)     = ttSTATS.tstat;
                
            end
            sigEpocPairs(iFreqBand,iTimeWindow,1:size(comparison,1))     = sigPairs;
            ttPvalEpocPairs(iFreqBand,iTimeWindow,1:size(comparison,1))  = ttPvalPairs;
            ttTstatEpocPairs(iFreqBand,iTimeWindow,1:size(comparison,1)) = ttTstat;
            
            
            %- save info for output structure
            if pVal<0.05,
                thisStr  = sprintf('%s %s, %d to %d ms', pStr, thisBand.name, thisTimeWindow);
                if isempty(statsStr), statsStr = thisStr;
                else                  statsStr = sprintf('%s\n%s',statsStr,thisStr); end
                %statsStr = sprintf('%s\n%s %s, %d to %d ms', statsStr, pStr, thisBand.name, thisTimeWindow);
                sigEpocs(iFreqBand, iTimeWindow) = 1;
            end
            warnStr = '';
            if sum(sigPairs)>0  & pVal>=0.05, warnStr = sprintf('\n   heads up: pair-wise difference but no overall difference! pVal=%.3f, %s, %d to %d ms', pVal, thisBand.name, thisTimeWindow); end        
            if sum(sigPairs)==0 & pVal< 0.05, warnStr = sprintf('\n   heads up: overall difference but no pairwise differences! pVal=%.3f, %s, %d to %d ms', pVal, thisBand.name, thisTimeWindow); end        
            
            statsStr                           = sprintf('%s%s',statsStr,warnStr);
            pValEpocs(iFreqBand,  iTimeWindow) = pVal;
            varExplEpocs(iFreqBand, iTimeWindow) = varExpl;
            labelEpocs{iFreqBand, iTimeWindow} = sprintf('%s, %d to %d ms', thisBand.name, thisTimeWindow);
            labelBands{iFreqBand}              = thisBand.name;
            
            
            %- output graphical stats and multiple comparison when time windows are selected strategically (instead of using sliding window)
            if SLIDING_PLOT==0,
                %- output overall comparison test result
                hStat = text(ax(2)-.1*range(ax(1:2)),0.95, pStr);
                set(hStat,'fontsize',fS);

                %- integrated multiple-comparison plot... just plot scaled version of mean estimates
                meanRange = [min(means(:,1)-means(:,2)) max(means(:,1)+means(:,2))];
                xScales = [xLimVal./meanRange]; xScale = 0.90*min(xScales(find(xScales>0))); if isempty(xScale), xScale=1; end;
                for thisTrig=1:numUniqueTrig,
                    yVal = 1.30 - thisTrig*.05;
                    xVal = means(thisTrig,1)+[-1 1]*means(thisTrig,2);
                    xVal = xVal*xScale; %%scale up so differences are visible
                    
                    hL = line(xVal,[1 1]*yVal);
                    set(hL,'LineWidth',3,'color',cTrigAr(thisTrig));
                    hP = plot(mean(xVal),yVal,'o');
                    set(hP,'MarkerSize',2,'MarkerEdgeColor',cTrigAr(thisTrig),'MarkerFaceColor',cTrigAr(thisTrig));
                end
                hL = line([0 0],[yVal 1.25],'color','k');
                
                
                %- seperate multiple comparisons plot
                if pVal<0.05 & PLOT_SEPARATE_MULTICOMPARE,
                    %- create new figure for multcompare output
                    if HIDE_FIGURES==1, figure('visible','off');
                    else                figure;                 end
                    pause(0.2);
                    
                    [comparison, means, hFig, groupNames] = multcompare(stats);
                    title(sprintf('%s :: %d to %d ms',thisBand.name, thisTimeWindow));
                    figList = [figList hFig];
                    
                    %- go back to the main figure
                    set(0,'CurrentFigure', thisFigNum); %-back to main figure... use this method to keep hidden and on bottom
                    pause(0.2);                         %-give the fig a chance to come to front... otherwise wrong fig can get updated
                end
            end

            
        end %-time window loop
        
    end %-freq band loop
    
end %-channel loop (each channel gets its own fig)

if length(statsStr)==0, statsStr = 'no sig differences'; end


%- compute x and y values to label times of significance; valid even if just doing epocs (instead of sliding window)
pSigX = mean(eventsAveWindowMS,2)/1000;
pSigY = sigEpocs;
pSigYtick = [];
for iFreqBand = 1:length(freqBandAr),
    pSigYtick = [pSigYtick iFreqBand];
    pSigY(iFreqBand,:) = sigEpocs(iFreqBand,:)*iFreqBand;
end
pSigY(find(pSigY==0)) = nan;


%- stats Structure outputs.
statsStruct.subj          = subj;
statsStruct.session       = session;

statsStruct.thisChanStr   = thisChanStr;
statsStruct.chanSubjStr   = sprintf('%s [%s]', thisChanStr, subj);
statsStruct.labelBands    = labelBands;
statsStruct.labelEpocs    = labelEpocs;
statsStruct.labelPairs    = labelPairs;

statsStruct.statsStr      = statsStr;
statsStruct.pValEpocs     = pValEpocs;
statsStruct.varExplEpocs  = varExplEpocs;
statsStruct.sigEpocs      = sigEpocs;
statsStruct.sigEpocPairs  = sigEpocPairs;
statsStruct.ttPvalEpocPairs  = ttPvalEpocPairs;
statsStruct.ttTstatEpocPairs = ttTstatEpocPairs;

statsStruct.plotSigX      = pSigX;
statsStruct.plotSigY      = pSigY;
statsStruct.plotSigYtick  = pSigYtick;
statsStruct.plotSigYtickL = labelBands;

statsStruct.waveletFreqs  = waveletFreqs;
statsStruct.freqBandAr    = freqBandAr;
statsStruct.eventsAveWindowMS = eventsAveWindowMS;
statsStruct.eventsTriggerXlim = eventsTriggerXlim;  %- equal to predefined windows or sliding window 


if SLIDING_PLOT==1,
    thisFigNum = 5000+diff(eventsAveWindowMS(1,1:2))+FIG_OFFSET;
    if ~ishghandle(thisFigNum), figure(thisFigNum);
    else                    set(0,'CurrentFigure', thisFigNum); end
    if HIDE_FIGURES==1,     set(thisFigNum,'visible','off');    
    else                    set(thisFigNum,'visible','on');     end
    
    clf
    set(gcf,'color','w');
    set(gcf,'position',[1480   500+diff(eventsAveWindowMS(1,1:2))        1080        120]); %
                
    
    plot(pSigX,pSigY,'b*'); hold on;
    
    title(titleStr, 'fontsize',figFontAx);
    set(gca,'ytick',pSigYtick,'yticklabel',labelBands,'YDir','reverse');
    set(gca,'tickdir','out','Box','off');
    set(gca,'fontsize',figFontAx-8)
    %axis tight
    axis([eventsTriggerXlim(1) eventsTriggerXlim(2) min(pSigYtick) max(pSigYtick)]) 
    grid on
    %set(gca,'xlim',eventsTriggerXlim);
    xlabel(sprintf('Time (s) :: [%d ms sliding windows]', diff(eventsAveWindowMS(1,1:2))));
    
    figList = gcf;
    
    pause(0.5); %- require a small delay so that the next function call doesn't accidentally use this figure handel... weird!
end
