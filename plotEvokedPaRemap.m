function [statsStruct, responseStruct, mainFigHandle] = plotEvokedPaRemap(data)
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
% trigType        = data.trigType;
eventsTrig      = data.eventsTrig;
trigTypeStr     = data.trigTypeStr;
trigZeroMS      = data.trigZeroMS;
% metaYstr        = data.metaYstr;
eventsTriggerXlim = data.eventsTriggerXlim;

waveT           = data.waveT;
wavesSft        = data.wavesSft;
wavesSftG       = data.wavesSftG;
powerMatZ       = data.powerMatZ; % numchannels x numevents x numfreqs x duration(ms)
powerMat        = data.powerMat;
phaseMat        = data.phaseMat;

numChan = size(powerMatZ, 1); % get the number of channels used

% store the channel string we're using (e.g. MST2-global)
if numChan==1 && length(chanStr)~=1, %- numChan~=length(chanStr) means processed sequentially... make sure title refects correct channel
    %keyboard
    chanStr{1} = thisChanStr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-----------------------   Evoked responses: comparing different conditions      ------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavePlotSelect = wavesSft;             % 
waveYlabel = 'Evoked Response (norm)'; % y axis label for ER plot

% WHAT DOES THIS PART MEAN? identify and trim high-count events?


% loop over channels
if numUniqueTrig > 1
    for chanNum = 1:numChan % numChan always = 1 when calling function - plot evoked one at a time
        %-- Select Figure and open/hide --%
        %set figure number based on offset handed in
        thisFigNum = 1000 + chanNum + FIG_OFFSET;
        
        if ~ishghandle(thisFigNum) % if no figure created, create it
            figure(thisFigNum);
        else
            set(0, 'CurrentFigure', thisFigNum);
        end
        if HIDE_FIGURES == 1
            set(thisFigNum, 'visible', 'off');
        end
        
        clf % clear the figure
        set(gcf, 'color', 'w')
        cTrigAr = 'rkbgmcy';
        
        for thisTrig=1:numUniqueTrig
            iTrig = find(trigType==uniqueTrigType(thisTirg)); % indices of the triggers
%             if thisTrig == iTrigTrim   % trim excess events add here
            
            cTrig = cTrigAr(thisTrig);
            if strcmp(THIS_TRIGGER(1:6), 'sample')
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
            
            events = eventsTrig(iTrig);
            yValues = thisTrig;
            tOffSet = trigZeroMS(iTrig); % time offset (ms)
            
            % event start and stop time
            eTimeOnS  = ([events.mstime]-tOffset)/1000;
            eTimeOffS = ([events.mstime]+[events.msDuration]-tOffset)/1000;
            
            hP = plot(eTimeOnS,yValues,'b.','MarkerSize',10); hold on; % all events with a blue '+' to indicate start
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',4); % all events get a line indicating duration
            
            % next (meta) event start 
            eTimeEnxtS = ([events.msMetaEventNext]-tOffset)/1000;   % time that next meta event starts
            hP = plot(eTimeEnxtS,yValues-.1,'k>','MarkerSize',6);
            
            %-- other event markers that can be used with trigger event or meta event
            % time cross is presented
            
            
            %-- label trigger type
            cleanTrigTypeStr = regexprep(thisTrigTypeStr,' \[.*\]', ''); %- remove bracketed expression, just for top panel
            hText = text(min(eventsTriggerXlim)-.05*abs(min(eventsTriggerXlim)), mean(yValues), cleanTrigTypeStr);
            set(hText,'HorizontalAlignment','right','FontSize',figFontAx);
            
            set(gca,'ydir','reverse','XAxisLocation','bottom','tickdir','out','fontsize',figFontAx)
            set(gca,'ylim',[1-.2 length(uniqueTrigType)+.2],'ytick',[1:length(uniqueTrigType)],'xlim',[eventsTriggerXlim], 'ydir','reverse')
            axis off
            box off
            
            %-- figure title
            % titleStr is like: NIH034 : Meta Session[0-3] : MST-034 Global
            titleStr = sprintf('%s : %s : channel %s', subj, session, chanStr{chanNum});
            titleStr(find(titleStr=='_'))=' ';
            title(titleStr, 'fontsize', figFontAx)
            
            %%----------------- plot evoked responses--------------------------------------------------------
            subplot(2+numUniqueTrig,1,2)
            axWave = gca;
            
            %-compute the evoked response for each trigger type
            wavePlot = wavePlotSelect(:,iTrig,:);
            clear evokedRespMu evokedRespSEM;
            evokedRespMu(1:size(wavePlot,3))  = mean(wavePlot(chanNum,:,:),2);
            evokedRespSEM(1:size(wavePlot,3)) = std(wavePlot(chanNum,:,:),0,2)./sqrt(size(wavePlot,2));
            
            % plot the evoked response average
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
            
            ylabel(waveYlabel,'fontsize',figFontAx)
            
            if numChan==1,
                axis tight;
                ax = axis;
                set(gca,'ylim',[-1 1]*max(abs(ax(3:4)))); %make y-limits symetric
                %    set(gca,'ylim',waveYlim);
            end
            
            
        end
    end
end