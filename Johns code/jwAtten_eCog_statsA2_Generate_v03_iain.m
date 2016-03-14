%function jwAtten_eCog_statsA_Generate( subj, sessNum, THIS_TRIGGER, THIS_REF_TYPE, FIG_OFFSET )
%
%   jwAttn_eCog_v01a.m
%
%         -- select subset of behavioral events
%         -- get waveform, power, phase for a time window around the event of interest
%         --
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%%% How to choose sessions?
% 1. Can clump vs. look at individual sesssions
% How to selectr trigger type?

%%- SUBJECT AND BLOCK SELECTION
% subj='NIH016';  sessNum = [0:4];  %[0:4] forced choice for all sessions.  * vs non-* only sig in session [2:3]
subj='NIH034';  sessNum = [0:3];  %[4:9 [0:10] session 10 not OK for physio; [sess 0-3,10 FC]; [sess 4-9 RECOG]
% subj='NIH019';  sessNum = [1:4 6:9];  %[0:9] forced choice all sessions. poor performance sess 0 and 5
% subj='NIH020';  sessNum = [0:3];  %[0:3] all look great. all sessions were recognition
% subj='NIH022';  sessNum = [0:3];  %[0:3] all look OK. all sessions were recognition
%subj='NIH024';  sessNum = [7];  %[0:3] all look OK. all sessions were recognition

%%- EXTRACTION OPTIONS
TRIGGER_TYPES   = {'sampleStart' 'testStart' 'testEnd' 'asterisk' 'sampleVSastStart' 'sampleVSastEnd' 'fixation' 'blockStart'}; %-fixation and blockStart not ready for physio analysis
THIS_TRIGGER    = TRIGGER_TYPES{1};   %%%%%%%%%% SELECT TRIGGER TYPE FOR EXTRACTION AND ANALYSIS


REF_TYPES       = {'noreref', 'bipolar', 'global'};
THIS_REF_TYPE   = REF_TYPES{3}; % (1) noreref, (2) bipolar, (3) laplacian


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    SET GLOBAL VARIABLES AND PARAMETERS                          ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HIDE_FIGURES    = 0;
USE_CHAN_SUBSET = 1; %0=all channels (not the subset); >=1 means process than many of the subset
% channelIWant = [1, 2, 3];

%%- FILTERING OPTIONS
BP_FILTER_RAW                 = 1;  %-0 or 1: apply a bandpass filter to the raw traces (1-499 hz)
PROCESS_CHANNELS_SEQUENTIALLY = 1;  %0 or 1:  0 means extract all at once, 1 means sequentially


%%- PLOT PARAMETERS
figFontAx       = 18;    if ispc, figFontAx = figFontAx-5; end
if ~exist('FIG_OFFSET','var'), FIG_OFFSET = 0; end %- default to 0, but if called as function then allow that value to persist
FIG_EVENTS      = 100+FIG_OFFSET;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n***********  SUBJECT: %s [%d-%d]  *****  TRIGGER: "%s"  *******  REFERENCE: "%s"  ***********\n', subj,min(sessNum),max(sessNum),THIS_TRIGGER,THIS_REF_TYPE); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end


FIG_TEMP_ROOT = '.';                 %- place the output figures in a local directory (eats up dropbox space)
FIG_TEMP_ROOT = [eegRootDir(1:end-4) '_jwAnalysis/'];  %- place the output figures near the FIGS_KEEP folder 


if sessNum == -1 | length(sessNum)>1,
    session = 'Meta Session [all]';
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/attentionTask/');
    sessStr = '[all]';
else
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/attentionTask/', session);
    sessStr = sprintf('[%d]',sessNum);
end

subjDir = fullfileEEG(eegRootDir,subj);
docsDir = fullfileEEG(subjDir,'docs');
talDir  = fullfileEEG(subjDir,'tal');
defaultEEGfile = fullfileEEG('/Volumes/Shares/FRNU/data/eeg/',subj,'/eeg.reref/');  % default event eegfile fields point here... switch to local before loading

events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('\nLoaded %d events from %s\n', length(events), behDir);

if length(sessNum)>1,
    session = sprintf('Meta Session [%d-%d]', min(sessNum), max(sessNum));
    sessStr = sprintf('[%d-%d]', min(sessNum), max(sessNum));
    cutSess = setdiff([events.sessionNum],sessNum);
    if length(cutSess)>0,
        keepSess = ismember([events.sessionNum],sessNum);
        events = events(find(keepSess));
        fprintf('--heads up: eliminated %d sessions (%d total events) from events structure--\n', length(cutSess), length(keepSess)-sum(keepSess));
    end
elseif sessNum==-1,
    session = sprintf('Meta Session [%d-%d]', min([events.sessionNum]), max([events.sessionNum]));
    sessStr = sprintf('[%d-%d]', min([events.sessionNum]), max([events.sessionNum]));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------     Create the meta and trigger events and timing variables     ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iEv=1:length(events), events(iEv).iOriginal = iEv; end;
eventsUnfilt = events; % creat unedited copy of loaded events structure


%-------- trigger on entire session, or set of blocks
events          = eventsUnfilt;
%iSrt            = find( [events.blockCount]>=1 & [events.blockCount]<=2);
iSrt            = 1:length(events); % uncomment this line to grab all blocks
blockEventsMeta = events(iSrt);
cutSessStart    = blockEventsMeta( ~strcmp({blockEventsMeta.type},'SESS_START') );
iBlockStarts    = [1 find(diff([cutSessStart.blockCount]))+1];
blockEventsTrig = cutSessStart(iBlockStarts);

for iBlk=1:length(blockEventsTrig)-1,
    blockEventsTrig(iBlk).msDuration = blockEventsTrig(iBlk+1).mstime - blockEventsTrig(iBlk).mstime;
end
blockEventsTrig(end).msDuration =  blockEventsMeta(end).mstime - blockEventsTrig(end).mstime;

% extract blocks separately or as one long event
CONCATINATE_BLOCKS = 0;
if length(blockEventsTrig)>1 & CONCATINATE_BLOCKS>0,
    for iBlk=1:length(blockEventsMeta),
        blockEventsMeta(iBlk).blockCount = 1;
        blockEventsMeta(iBlk).msBlockStart = blockEventsMeta(1).msBlockStart;
    end
    blockEventsTrig = blockEventsTrig(1);
    blockEventsTrig.msDuration = blockEventsMeta(end).mstime - blockEventsTrig.mstime;
end


%-------- trigger on sample word events (includes *correct and non-* C and E)
events       = eventsUnfilt;
iSrt         = find( [events.asteriskType]>=1 & [events.asteriskType]<=3  & ~isnan([events.sampleCount])  & [events.responseCorrect]==1 );
events       = events(iSrt);
[aSrt, iSrt] = sort([events.asteriskType]);
events       = events(iSrt);
noAstError   = eventsUnfilt(find( [eventsUnfilt.asteriskType]==3 & ~isnan([eventsUnfilt.sampleCount])  & [eventsUnfilt.responseCorrect]==0 ));
for iEv=1:length(noAstError), noAstError(iEv).asteriskType = 13; end % force a "13" to indicate non-asterisk error... (non ast=3 + 10 for error
events(end+1:end+length(noAstError)) = noAstError;

%- add another copy of the correct asterisk trials, split in half by fast and slow reaction times and labeled as such
noAstCorrSrt = eventsUnfilt(find( [eventsUnfilt.asteriskType]==3 & ~isnan([eventsUnfilt.sampleCount])  & [eventsUnfilt.responseCorrect]==1 ));
[aSrt, iSrt] = sort([noAstCorrSrt.RT]);  
noAstCorrSrt = noAstCorrSrt(iSrt);
iFast = find([noAstCorrSrt.RT] <= median(unique([noAstCorrSrt.RT]))); %- this way allows for different numbers of elements in each meta event (should always be 2 for non-* encoding, but what the hay)
iSlow = find([noAstCorrSrt.RT] >  median(unique([noAstCorrSrt.RT])));  
for iEv=iFast, 
    noAstCorrSrt(iEv).sampleCount = noAstCorrSrt(iEv).sampleCount+max([events.sampleCount])*10;
    noAstCorrSrt(iEv).asteriskType = 23;  % force a "23" to indicate non-asterisk correct that has a fast reaction time... (non ast=3 + 20 for fast correct
    noAstCorrSrt(iEv).asteriskStr = '*nFast';
end 
for iEv=iSlow, 
    noAstCorrSrt(iEv).sampleCount = noAstCorrSrt(iEv).sampleCount+max([events.sampleCount])*10*2;  % multiply by 2 to pass all the samples from iFast
    noAstCorrSrt(iEv).asteriskType = 33;  % force a "33" to indicate non-asterisk correct that has a fast reaction time... (non ast=3 + 30 for slow correct
    noAstCorrSrt(iEv).asteriskStr = '*nSlow';
end 
events(end+1:end+length(noAstCorrSrt)) = noAstCorrSrt;

sampEventsMeta = events;  % includes assocaited + and *
sampEventsTrig = sampEventsMeta(find([sampEventsMeta.isSample]));


%-------- trigger on test word events
events       = eventsUnfilt;
iSrt         = find( [events.asteriskType]>=1 & [events.asteriskType]<=3  & ~isnan([events.testCount])  & [events.responseCorrect]==1 );
events       = events(iSrt);
[aSrt, iSrt] = sort([events.asteriskType]);
events       = events(iSrt);

noAstError   = eventsUnfilt(find( [eventsUnfilt.asteriskType]==3 & ~isnan([eventsUnfilt.testCount])  & [eventsUnfilt.responseCorrect]==0 ));
for iEv=1:length(noAstError), noAstError(iEv).asteriskType = 13; end % force a "13" to indicate non-asterisk error... (non ast=3 + 10 for error
events(end+1:end+length(noAstError)) = noAstError;

foilCorrect  = eventsUnfilt(find( [eventsUnfilt.asteriskType]==4 & ~isnan([eventsUnfilt.testCount])  & [eventsUnfilt.responseCorrect]==1 ));
for iEv=1:length(foilCorrect), foilCorrect(iEv).asteriskType = 14; end % force a "14" to indicate foil error... (foil=4 + 10 for error
events(end+1:end+length(foilCorrect)) = foilCorrect;

allAstCorr   = eventsUnfilt(find( [eventsUnfilt.asteriskType]>=1 & [eventsUnfilt.asteriskType]<=2 & ~isnan([eventsUnfilt.testCount])  & [eventsUnfilt.responseCorrect]==1 ));
for iEv=1:length(allAstCorr), 
    allAstCorr(iEv).asteriskType = 20;
    allAstCorr(iEv).testCount = allAstCorr(iEv).testCount+max([events.testCount])*10;
    allAstCorr(iEv).asteriskStr = '*all';
end % force a "20" to indicate lumped asterisk before and after correct...
events(end+1:end+length(allAstCorr)) = allAstCorr;

% %- add another copy of the correct asterisk trials, split in half by fast and slow reaction times and labeled as such
% noAstCorrSrt = eventsUnfilt(find( [eventsUnfilt.asteriskType]==3 & ~isnan([eventsUnfilt.testCount])  & [eventsUnfilt.responseCorrect]==1 ));
% [aSrt, iSrt] = sort([noAstCorrSrt.RT]);  
% noAstCorrSrt = noAstCorrSrt(iSrt);
% iFast = find([noAstCorrSrt.RT] <= median(unique([noAstCorrSrt.RT]))); %- this way allows for different numbers of elements in each meta event (should always be 2 for non-* encoding, but what the hay)
% iSlow = find([noAstCorrSrt.RT] >  median(unique([noAstCorrSrt.RT])));  
% for iEv=iFast, 
%     noAstCorrSrt(iEv).testCount = noAstCorrSrt(iEv).testCount+max([events.testCount])*10;
%     noAstCorrSrt(iEv).asteriskType = 23;  % force a "23" to indicate non-asterisk correct that has a fast reaction time... (non ast=3 + 20 for fast correct
%     noAstCorrSrt(iEv).asteriskStr = '*nFast';
% end 
% for iEv=iSlow, 
%     noAstCorrSrt(iEv).testCount = noAstCorrSrt(iEv).testCount+max([events.testCount])*10*2;  % multiply by 2 to pass all the samples from iFast
%     noAstCorrSrt(iEv).asteriskType = 33;  % force a "33" to indicate non-asterisk correct that has a fast reaction time... (non ast=3 + 30 for slow correct
%     noAstCorrSrt(iEv).asteriskStr = '*nSlow';
% end 
% events(end+1:end+length(noAstCorrSrt)) = noAstCorrSrt;

%foilError   = eventsUnfilt(find( [eventsUnfilt.asteriskType]==4 & ~isnan([eventsUnfilt.testCount])  & [eventsUnfilt.responseCorrect]==0 ));
%for iEv=1:length(foilError), foilError(iEv).asteriskType = 14; end % force a "14" to indicate foil error... (foil=3 + 10 for error
%events(end+1:end+length(foilError)) = foilError;
testEventsMeta = events;  % includes associated +
testEventsTrig = testEventsMeta(find([testEventsMeta.isTest]));


%-------- trigger on fixation cross events
events       = eventsUnfilt;
iSrt         = find( [events.isCross]>=1 );
events       = events(iSrt);
[aSrt, iSrt] = sort([events.isCross]);
crossEvents  = events(iSrt);


%-------- trigger on asterisk events
events       = eventsUnfilt;
iSrt         = find( [events.isAsterisk]>=1 & [events.responseCorrect]==1 );
events       = events(iSrt);
[aSrt, iSrt] = sort([events.isAsterisk]);
events       = events(iSrt);
crossNoAstCor = eventsUnfilt(find( [eventsUnfilt.responseCorrect]==1 & [eventsUnfilt.asteriskType]==3 & [eventsUnfilt.isCross]==1 ));
events(end+1:end+length(crossNoAstCor)) = crossNoAstCor;
crossNoAstErr = eventsUnfilt(find( [eventsUnfilt.responseCorrect]==0 & [eventsUnfilt.asteriskType]==3 & [eventsUnfilt.isCross]==1 ));
for iEv=1:length(crossNoAstErr), crossNoAstErr(iEv).asteriskType = 13; end % force a "13" to indicate non-asterisk error... (non ast=3 + 10 for error
events(end+1:end+length(crossNoAstErr)) = crossNoAstErr;
astEvents    = events;


%-------- trigger on sample vs asterisk: sample offset time (for *before or *none) AND asterisk offset time (for *after)
events       = eventsUnfilt;
iSrt         = find( [events.responseCorrect]==1 & ([events.isAsterisk]==2 | ([events.asteriskType]==1 & [events.isSample]==1) | ([events.asteriskType]==3 & [events.isSample]==1)) ); %- just the asterisk after events...
events       = events(iSrt);
[aSrt, iSrt] = sort([events.asteriskType]);
events       = events(iSrt);
noAstError   = eventsUnfilt(find( [eventsUnfilt.responseCorrect]==0 & [eventsUnfilt.asteriskType]==3 & [eventsUnfilt.isSample]==1 ));
for iEv=1:length(noAstError), noAstError(iEv).asteriskType = 13; end % force a "13" to indicate non-asterisk error... (non ast=3 + 10 for error
events(end+1:end+length(noAstError)) = noAstError;
for iEv=1:length(events), events(iEv).mstimeEnd = events(iEv).mstime + events(iEv).msDuration; end
sampVSastEvents = events;  % this meta does not include associated stuff... too tricky when comparing events at different time points


switch THIS_TRIGGER,
    case 'blockStart'
        eventsMeta = blockEventsMeta;   metaYval = [eventsMeta.blockCount];    metaZeroMS = [eventsMeta.msBlockStart]; metaYstr = 'Block Start';      metaYtransition = find(diff([eventsMeta.blockCount]));
        eventsTrig = blockEventsTrig;   trigType = [eventsTrig.blockCount];     trigZeroMS = [eventsMeta.msBlockStart]; trigTypeStr = 'blockCount';
        eventsTriggerXlim = [-10 max([eventsMeta.mstime]-metaZeroMS)/1000+ 10]; %10 sec before and after each block
        eventsAveWindowMS = [0 1000]; % list of time windows over which EEG data is averaged for t-tests
    case 'sampleStart'
        disp(['Looking at trigger: ', THIS_TRIGGER]);
        eventsMeta = sampEventsMeta;    
        metaYval = [eventsMeta.sampleCount];    
        metaZeroMS = [eventsMeta.msTextStart];  metaYstr = 'Sample Onset';     metaYtransition = find(diff([eventsMeta.asteriskType]));
        eventsTrig = sampEventsTrig;    trigType = [eventsTrig.asteriskType];   trigZeroMS = [eventsTrig.msTextStart];  trigTypeStr = 'asteriskType';
        eventsTriggerXlim = [-2.25 5.25];
        eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
        %eventsTriggerXlim = [min(min(eventsAveWindowMS)) max(max(eventsAveWindowMS))]/1000;
    case 'testStart'
        eventsMeta = testEventsMeta;    metaYval = [eventsMeta.testCount];      metaZeroMS = [eventsMeta.msTextStart];  metaYstr = 'Test Onset';       metaYtransition = find(diff([eventsMeta.asteriskType]));
        eventsTrig = testEventsTrig;    trigType = [eventsTrig.asteriskType];   trigZeroMS = [eventsTrig.msTextStart];  trigTypeStr = 'asteriskType';
        eventsTriggerXlim = [-1.5 4.5];
        eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
    case 'testEnd'
        eventsMeta = testEventsMeta;    metaYval = [eventsMeta.testCount];      metaZeroMS = [eventsMeta.msTextStart] + double([eventsMeta.RT]);  metaYstr = 'Test Response';       metaYtransition = find(diff([eventsMeta.asteriskType]));
        eventsTrig = testEventsTrig;    trigType = [eventsTrig.asteriskType];   trigZeroMS = [eventsTrig.msTextStart] + double([eventsTrig.RT]);  trigTypeStr = 'asteriskType';
        %eventsMeta = testEventsMeta;    metaYval = [eventsMeta.testCount];      metaZeroMS = [eventsMeta.msTextStop];  metaYstr = 'Test Offset';       metaYtransition = find(diff([eventsMeta.asteriskType]));
        %eventsTrig = testEventsTrig;    trigType = [eventsTrig.asteriskType];   trigZeroMS = [eventsTrig.msTextStop];  trigTypeStr = 'asteriskType';
        eventsTriggerXlim = [-4.5 1.5];
        eventsAveWindowMS = [-3000 -2500; -2500 -2000; -2000 -1500; -1500 1000; -1000 -500; -500 0; 0 500; 500 1000]; % list of time windows over which EEG data is averaged for t-tests
    case 'fixation'
        eventsMeta = crossEvents;       
        metaYval = [eventsMeta.eventCount];     
        metaZeroMS = [eventsMeta.mstime];       
        metaYstr = 'Fixation Onset';   
        metaYtransition = find(diff([eventsMeta.isCross]));
        eventsTrig = crossEvents;       
        trigType = [eventsTrig.isCross];        
        trigZeroMS = [eventsTrig.mstime];       
        trigTypeStr = 'sampleOrTest';
        eventsTriggerXlim = [-.5 7];
        eventsAveWindowMS = [0 1000];                   % list of time windows over which EEG data is averaged for t-tests
    case 'asterisk'
        eventsMeta = astEvents;   
        metaYval = [eventsMeta.eventCount]; 
        metaZeroMS = [eventsMeta.mstime];  
        metaYstr = 'Asterisk Onset';   
        metaYtransition = find(diff([eventsMeta.asteriskType]));
        eventsTrig = astEvents;       
        trigType = [eventsTrig.asteriskType];   
        trigZeroMS = [eventsTrig.mstime];       
        trigTypeStr = 'asteriskType';
        eventsTriggerXlim = [-3 5.5];
        eventsAveWindowMS = [-500 0; 0 500; 500 1000];  % list of time windows over which EEG data is averaged for t-tests
    case 'sampleVSastStart'
        eventsMeta = sampVSastEvents;   metaYval = [eventsMeta.eventCount];     metaZeroMS = [eventsMeta.mstime];       metaYstr = 'Sample vs Asterisk [Onset]'; metaYtransition = find(diff([eventsMeta.asteriskType]));
        eventsTrig = sampVSastEvents;   
        trigType = [eventsTrig.asteriskType];   
        trigZeroMS = [eventsTrig.mstime];       trigTypeStr = 'asteriskType';
        eventsTriggerXlim = [-3 4.5];
        eventsAveWindowMS = [-500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 2500 3000; 3000 3500];  % list of time windows over which EEG data is averaged for t-tests
    case 'sampleVSastEnd'
        eventsMeta = sampVSastEvents;   metaYval = [eventsMeta.eventCount];     metaZeroMS = [eventsMeta.mstimeEnd];    metaYstr = 'Sample vs Asterisk [Offset]'; metaYtransition = find(diff([eventsMeta.asteriskType]));
        eventsTrig = sampVSastEvents;   
        trigType = [eventsTrig.asteriskType];   
        trigZeroMS = [eventsTrig.mstimeEnd];    trigTypeStr = 'asteriskType';
        eventsTriggerXlim = [-3.5 3.5];
        eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 2500 3000];  % list of time windows over which EEG data is averaged for t-tests
    otherwise
        error('no event trigger selected');
end

%- need this fix so gete_ms grabs correctly-offset signals
for iEv=1:length(eventsTrig), eventsTrig(iEv).eegoffset = eventsTrig(iEv).eegoffset + trigZeroMS(iEv) - eventsTrig(iEv).mstime; end

trigType = double(trigType); % so find can be used

% modify meta y-values so they increment without spaces... use unique to maintain meta-event identity
[C,IA,IC] = unique(metaYval, 'stable');
clear metaYval; metaYval(1,1:length(IC)) = IC;
if length(eventsTrig)~=length(C), fprintf('heads up: number of triggered events does not match number of meta events\n'); end


% create sliding windows for anova analysis
winWidth = 100;
winDT    = 50;
slidingWindow100 = [[eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'-winWidth/2 [eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'+winWidth/2];  %-new way... centered for both window sizes
winWidth = 250;
winDT    = 50;
slidingWindow250 = [[eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'-winWidth/2 [eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'+winWidth/2];  %-new way... centered for both window sizes
winWidth = 500;
winDT    = 50;
slidingWindow500 = [[eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'-winWidth/2 [eventsTriggerXlim(1)*1000:winDT:eventsTriggerXlim(2)*1000]'+winWidth/2];  %-new way... centered for both window sizes

%if strcmp(THIS_TRIGGER,'testStart') useSmallWinDT = 1; else useSmallWinDT = 0; end  %- only use the 500 ms window unless "useSmallWinDT" is equal to 1
useSmallWinDT = 1;


% %%- provide option for normalizing number of events in each trigger type here... trim from meta and from trigger
% NORMALIZE_EVENT_COUNTS = 0;
% if NORMALIZE_EVENT_COUNTS==1,
%
%     uniqueTrigType = unique(trigType);
%     numUniqueTrig  = length(uniqueTrigType);
%     for thisTrig=1:numUniqueTrig,  numTrig(thisTrig) = length(find(trigType==uniqueTrigType(thisTrig))); end
%     numTrigSort = sort(numTrig,'descend');
%
%
%     %- option1: just reduce the size of the highest count
%     iTrigTrim = -1;
%     if length(numTrigSort)>1 & numTrigSort(1)>1.5*numTrigSort(2),
%         iKeep=unique(round(1:numTrigSort(1)/numTrigSort(2):numTrigSort(1)));
%         iTrigTrim = find(numTrig==max(numTrig));
%         fprintf(' note: category %d has way more events than the others (%d events), will trim to %d events\n',iTrigTrim,numTrigSort(1),length(iKeep));
%     end
%
%
%     %- option2: make all equal to the smallest count [probably excessive]
%
%
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------------------------------       plot meta event data        -------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_EVENTS=0;
if (PLOT_EVENTS)
    figure(FIG_EVENTS); clf
    set(gcf,'color','w')
    
    
    events  = eventsMeta;
    yValues = metaYval;
    tOffset = metaZeroMS;
    
    
    %-- event start time
    eTimeOnS  = ([events.mstime]-tOffset)/1000;
    plot(eTimeOnS,yValues-.15,'b.'); hold on;               % all events get a blue + to indicate start
    %-- event duration
    eTimeOffS  = ([events.mstime]+[events.msDuration]-tOffset)/1000;
    line([eTimeOnS; eTimeOffS], [yValues-.15; yValues-.15]) % all events get a line indicating duration
    %-- next (meta) event start time
    eTimeEnxtS = ([events.msMetaEventNext]-tOffset)/1000;      % time that next meta event starts
    plot(eTimeEnxtS,yValues-.20,'>');
    
    %-- other event markers that can be used with trigger event or meta event
    eTimeTextOnS  = ([events.msTextStart] -tOffset)/1000;
    h = plot(eTimeTextOnS,yValues-.15,'ks');  % time text is presented
    eTimeCrosOnS  = ([events.msCrossStart]-tOffset)/1000;
    plot(eTimeCrosOnS,yValues-.15,'k+');  % time cross is presented
    eTimeAstrOnS  = ([events.msAsteriskStart]-tOffset)/1000;
    plot(eTimeAstrOnS,yValues-.15,'kp');  % time asterisk is presented
        
    % plot horizontal line indicating meta Y transitions (i.e., change in trigger event
    for iTrans = 1:length(metaYtransition)
        xVals = eventsTriggerXlim;
        yVals = [1 1]*yValues(metaYtransition(iTrans))+.5;
        plot(xVals,yVals,'r--', 'linewidth',2);
    end
    
    % label the asterisk type when triggering on samples (can extend this to test phase too)
    if strcmp(THIS_TRIGGER(1:6),'sample') | strcmp(THIS_TRIGGER(1:4),'test') | strcmp(THIS_TRIGGER,'asterisk') ,
        for iTrans = 1:length(metaYtransition),
            % add text ylabel
            if iTrans==1,
                yVal = yValues(metaYtransition(iTrans))/2;
            else
                yVal = mean( yValues(metaYtransition(iTrans-1:iTrans)) );
            end
            strCorr = 'EC';
            astStr = sprintf('%s (%s)', events(metaYtransition(iTrans)).asteriskStr, strCorr(events(metaYtransition(iTrans)).responseCorrect+1));
            hT = text(eventsTriggerXlim(1)-range(eventsTriggerXlim)*.1, yVal, astStr);
            set(hT,'Rotation',90,'fontsize',figFontAx,'fontweight','bold','horizontalalignment','center')
        end
        if metaYtransition(end)~=length(events),
            yVal = mean( [yValues(metaYtransition(iTrans)) yValues(end)] );
            astStr = sprintf('%s (%s)', events(end).asteriskStr, strCorr(events(end).responseCorrect+1));
            hT = text(eventsTriggerXlim(1)-range(eventsTriggerXlim)*.1, yVal, astStr);
            set(hT,'Rotation',90,'fontsize',figFontAx,'fontweight','bold','horizontalalignment','center')
        end
    else
        ylabel(metaYstr,'fontsize',figFontAx+7)
    end
    
    
    xlabel(sprintf('%s (s)',metaYstr),'fontsize',figFontAx+7)
    
    yTickAuto = get(gca,'ytick');  yTickAuto = unique([1 yTickAuto max(yValues)]);
    set(gca,'fontsize',figFontAx+2, 'ylim',[min(yValues)-.9 max(yValues)+.9], 'ytick',yTickAuto, 'ydir','reverse');
    
    set(gca,'xlim',eventsTriggerXlim)
    box off
    
    %- text labels on each event
    for iE = 1:length(events),
        strTxt   = events(iE).textStr;
        strAbrev = events(iE).textAbrev;
        if (~isempty( strTxt ))
            yOff=0; vAlign = 'baseline';
            
            if (strTxt(1)=='+' | strTxt(1)=='?')    tSize=25; tAlign='center';             strTxt = strTxt(1);
            elseif (strTxt(1)=='*')                 tSize=30; tAlign='left'; yOff=+.05;    strTxt = strTxt(1);
            else                                    tSize=12; tAlign='left';               strTxtOut(1) = {['\fontsize{15}\bf' strAbrev '\rm\fontsize{10}']}; strTxtOut(2) = {strTxt}; strTxtOut(3) = {sprintf('(%s)',events(iE).resultStr)}; strTxt=strTxtOut; end
            
            hT = text(eTimeOnS(iE),yValues(iE)+yOff,strTxt);
            set(hT,'FontSize',tSize, 'HorizontalAlignment',tAlign,'VerticalAlignment',vAlign);
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     plot the simplified (triggered) version too     %%%
    
    
    uniqueTrigType = unique(trigType);
    numUniqueTrig  = length(uniqueTrigType);
    for thisTrig=1:numUniqueTrig,  numTrig(thisTrig) = length(find(trigType==uniqueTrigType(thisTrig))); end
    numTrigSort = sort(numTrig,'descend');
    iTrigTrim = -1;
    if length(numTrigSort)>1 & numTrigSort(1)>1.5*numTrigSort(2),
        iKeep=unique(round(1:numTrigSort(1)/numTrigSort(2):numTrigSort(1)));
        iTrigTrim = find(numTrig==max(numTrig));
        fprintf(' note: category %d has way more events than the others (%d events), will trim to %d events\n',iTrigTrim,numTrigSort(1),length(iKeep));
    end
    
    
    PLOT_TRIGGERED_VERSION = 1;
    if numUniqueTrig>1  &  PLOT_TRIGGERED_VERSION & ~strcmp(THIS_TRIGGER,'blockStart'),
        
        
        mainAx = gca;
        set(mainAx,  'position',[.15 .07 .80 .80]);  %left, bottom, width, height
        subAx = axes('position',[.15 .89 .80 .05]);  %left, bottom, width, height
        
        cTrigAr = 'rkbgmcyr';
        
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
                thisTrigTypeStr = sprintf('%s (%s) [%d ev]',eventsTrig(iTrig(1)).asteriskStr, clnResultStr, length(iTrig));
            else
                thisTrigTypeStr = sprintf('%s %d [%d ev]', trigTypeStr, uniqueTrigType(thisTrig), length(iTrig));
            end
            
            
            events   = eventsTrig(iTrig);
            yValues  = thisTrig; % increment sequentially... don't use value of trigType, which can jump
            tOffset  = trigZeroMS(iTrig);
            
            
            %-- event start and stop time
            eTimeOnS  = ([events.mstime]-tOffset)/1000;
            eTimeOffS = ([events.mstime]+[events.msDuration]-tOffset)/1000;
            hP = plot(eTimeOnS,yValues,'b.'); hold on;               % all events get a blue + to indicate start
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',4); % all events get a line indicating duration
            %-- next (meta) event start time
            eTimeEnxtS = ([events.msMetaEventNext]-tOffset)/1000;   % time that next meta event starts
            plot(eTimeEnxtS,yValues-.1,'k>');
            
            %-- other event markers that can be used with trigger event or meta event
            % time cross is presented
            eTimeOnS  = ([events.msCrossStart]-tOffset)/1000;
            eTimeOffS = ([events.msCrossStop] -tOffset)/1000;
            plot(eTimeOnS,yValues-.15,'k+');  % time cross is presented
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',2); % all events get a line indicating duration
            % time sample/test text is presented
            eTimeOnS  = ([events.msTextStart]-tOffset)/1000;
            eTimeOffS = ([events.msTextStop] -tOffset)/1000;
            hP = plot(eTimeOnS,yValues-.15,'ks');
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',2); % all events get a line indicating duration
            % time asterisk is presented
            eTimeOnS  = ([events.msAsteriskStart]-tOffset)/1000;
            eTimeOffS = ([events.msAsteriskStop] -tOffset)/1000;
            plot(eTimeOnS,yValues-.15,'kp');
            hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',2); % all events get a line indicating duration
            
            %-- label trigger type
            hText = text(min(eventsTriggerXlim)-.001*abs(min(eventsTriggerXlim)), mean(yValues), thisTrigTypeStr);
            set(hText,'HorizontalAlignment','right','FontSize',figFontAx);
            
            set(gca,'ydir','reverse','XAxisLocation','bottom','tickdir','out','fontsize',figFontAx)
            set(gca,'ylim',[1-.2 length(uniqueTrigType)+.2],'ytick',[1:length(uniqueTrigType)],'xlim',[eventsTriggerXlim], 'ydir','reverse')
            set(gca,'xlim',eventsTriggerXlim)
            axis off
            %box off
        end
    end
    
    
    % title applied to the upper-most axis
    titleStr = sprintf('%s : %s',subj,session); titleStr(find(titleStr=='_'))=' ';
    title(titleStr, 'fontsize',figFontAx+7)
    
    if ismac, set(gcf,'position',[ 1000           1        1097        1345]); end
    
    if FIG_OFFSET==0, fprintf('\n\nPRESS ANY KEY TO BEGIN EXTRACTION\n\n');                   pause;
    else              fprintf('\nFIG_OFFSET non-zero... starting extraction in 5 seconds\n'); pause(5); 
    end
    
    if HIDE_FIGURES==1, set(gcf,'visible','off'); fprintf('... ALL FIGURES WILL BE HIDDEN ... \n\n'); pause(0.5); end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------------------------------      Create Channel List          -------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jackSheet = fullfileEEG(docsDir, 'jacksheetMaster.txt');
[chanNums chanTags] = textread(jackSheet,'%d%s%*s');

%%% always look at all electrodes... worry about "good" and "bad" later (bad means inter-ictal activity or seizure activity)
%- three referencing options:  noreref (should manually subtract reference channel), reref bioploar, and reref laplacian
%
switch THIS_REF_TYPE
    case 'noreref'
        fprintf('No rereferencing');
        chanFile      = [talDir '/leads.txt'];
        chanList      = textread(chanFile,'%d');
        chanList      = [chanList;  chanNums(find(strcmp(chanTags,'R1')));  chanNums(find(strcmp(chanTags,'R2'))); chanNums(find(strcmp(chanTags,'EKG1'))); chanNums(find(strcmp(chanTags,'EKG2')))];
        for iChan=1:size(chanList,1),
            %    chanStr{iChan} = sprintf('%d (%s noreref)', chanList(iChan), chanTags{find(chanNums==chanList(iChan))} );
            chanStr{iChan} = sprintf('%s-gnd', chanTags{find(chanNums==chanList(iChan))} );
        end
        chanRef1      = chanNums(find(strcmp(chanTags,'R1')));
        chanRef2      = chanNums(find(strcmp(chanTags,'R2')));
        if length(chanRef1)+length(chanRef2)~=2, fprintf('missing reference channel!'); keyboard; end
        chanRefs      = [chanRef1 chanRef2];
        chanRefs      = [chanRef1 chanRef2 chanNums(find(strcmp(chanTags,'EKG1')))]; %temporary--include EKG to check it out
        eventEEGpath  = '/eeg.noreref/';
        
        iChanListSub  = [1 2 31 45 chanRefs];   %G1, G2, LF1, AST1, R1, R2, EKG1
        
    case 'bipolar'
        fprintf('Bipolar referencing');
        chanFile      = [talDir '/leads_bp.txt'];
        [chan1 chan2] = textread(chanFile,'%d%*c%d');
        chanList      = [chan1 chan2];
        for iChan=1:size(chanList,1),
            %    chanStr{iChan} = sprintf('%d-%d (%s-%s)', chan1(iChan), chan2(iChan), chanTags{find(chanNums==chan1(iChan))}, chanTags{find(chanNums==chan2(iChan))} );
            chanStr{iChan} = sprintf('%s-%s', chanTags{find(chanNums==chan1(iChan))}, chanTags{find(chanNums==chan2(iChan))} );
        end
        chanRefs      = [];
        eventEEGpath  = '/eeg.reref/';
        
        iChanListSub  = [1, 3, 49, 60];                               %G1, G2, LF1, AST1
        iChanListSub  = [1, 3, 49, 52, 55, 60, 63, 66, 69, 72, 77];   %G1, G2, LF1, OF1, TT1, AST1, MST1, PST1, PPST1, TO1, TP1
        
    case 'global'
        fprintf('Global referencing');
        chanFile      = [talDir '/leads.txt'];
        chanList      = textread(chanFile,'%d');
        for iChan=1:size(chanList,1),
            %    chanStr{iChan} = sprintf('%d (%s-global)', chanList(iChan), chanTags{find(chanNums==chanList(iChan))} );
            chanStr{iChan} = sprintf('%s-global', chanTags{find(chanNums==chanList(iChan))} );
        end
        chanRefs      = [];
        eventEEGpath  = '/eeg.reref/';
        
        iChanListSub  = [48 1];            %G1, G2, LF1, AST1,
        try
            iChanListSub = channelIWant;
        catch e
        end
    otherwise
        fprintf('Error, no referencing scheme selected');
end


%%- select all channels, or part of the subset of channels
if USE_CHAN_SUBSET==0,
    iChanList = 1:size(chanList,1);  %all possible channels
else
    iChanList = iChanListSub(1:min([length(iChanListSub) USE_CHAN_SUBSET]));    %select subset of channels
end

%% Printing out how many channels we are processing
chanListUse = [];  chanStrUse = {};
for iChan=iChanList,
    chanListUse(end+1,:) = chanList(iChan,:);
    chanStrUse{end+1}    = chanStr{iChan};
end
chanList = chanListUse;
chanStr  = chanStrUse;
numChannels = size(chanList,1);
fprintf(' -- %d channels to process for %s : %s', numChannels, subj, session);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------------------------------      Create Frequency List and define Bands     -----------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%- gets the range of frequencies using eeganalparams
waveletFreqs = eeganalparams('freqs');
waveletWidth = eeganalparams('width');
for iFB=1:length(waveletFreqs), waveletFreqLabels{iFB} = sprintf('%.0f Hz', waveletFreqs(iFB)); end

%-array of frequency bands
freqBandAr(1).name    = 'delta';
freqBandAr(1).rangeF  = [2 4];          %[2 4]
freqBandAr(2).name    = 'theta';
freqBandAr(2).rangeF  = [4 8];          %[4 8]
freqBandAr(3).name    = 'alpha';
freqBandAr(3).rangeF  = [8 16];         %[8 12]
freqBandAr(4).name    = 'beta';
freqBandAr(4).rangeF  = [16 32];        %[12 30]
freqBandAr(5).name    = 'low gamma';
freqBandAr(5).rangeF  = [32 80];        %[30 70]
freqBandAr(6).name    = 'high gamma';
freqBandAr(6).rangeF  = [80 160];       %[70 150]
freqBandAr(7).name    = 'HFO';
freqBandAr(7).rangeF  = [160 400];      %[150 400]

for iFB=1:length(freqBandAr),
    freqBandAr(iFB).centerF = mean(freqBandAr(iFB).rangeF);
    %freqBandAr(iFB).label   = sprintf('%s-%.0fHz', freqBandAr(iFB).name(1:[ min( [length(freqBandAr(iFB).name), 6] )]), freqBandAr(iFB).centerF);
    freqBandAr(iFB).label   = sprintf('%s [%.0f-%.0f Hz]', freqBandAr(iFB).name, freqBandAr(iFB).rangeF);
end

freqBandYticks  = unique([freqBandAr(1:7).rangeF]);
for iFB=1:length(freqBandYticks), freqBandYtickLabels{iFB} = sprintf('%.0f Hz', freqBandYticks(iFB)); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------------------------------     Extract EEG data (waves and power)    -----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs to gete_ms:
eventTrigger    = eventsTrig;
eventOffsetMS   = eventsTriggerXlim(1)*1000;      % positive = after event time; negative = before event time
eventDurationMS = diff(eventsTriggerXlim)*1000;   % duration includes offset (i.e., if offset -500 and duration 1000, only 500 ms post event will be prsented)

% remap event pointer from default (server) to local copy of the EEG data
for iEvent=1:length(eventTrigger),
    eventTrigger(iEvent).eegfile = regexprep(eventTrigger(iEvent).eegfile,defaultEEGfile,fullfileEEG(subjDir,eventEEGpath));
end

OffsetMS        = eventOffsetMS;     % positive = after event time; negative = before event time
DurationMS      = eventDurationMS;   % duration includes offset (i.e., if offset -500 and duration 1000, only 500 ms post event will be prsented)
BufferMS        = 1000;              % grab excess data before/after event window so filters don't have edge effect
resampledrate   = 1000;              % don't resample... keep the 1kHz native sample rate


% apply a bandpass filter raw data? (i.e. pre-filter the wave?)
if BP_FILTER_RAW==1,
    preFiltFreq      = [1 499];   %[1 499] [2 250]
    preFiltType      = 'bandpass';
    preFiltOrder     = 2;
    preFiltStr       = sprintf('%s filter raw; %.1f - %.1f Hz',preFiltType,preFiltFreq);
    preFiltStrShort  = '_BPfilt';
    FIG_OFFSET       = FIG_OFFSET+100+round(preFiltFreq(2));  %keep this empty to avoid any filtering of the raw data
else
    preFiltFreq      = []; %keep this empty to avoid any filtering of the raw data
    preFiltType      = 'stop';
    preFiltOrder     = 1;
    preFiltStr       = 'Unfiltered raw traces';
    preFiltStrShort  = '_noFilt';
end

%-- pre-allocate memory (to make sure it can be done!)
if PROCESS_CHANNELS_SEQUENTIALLY==0,  numChanPrealloc = numChannels;  else  numChanPrealloc = 1; end
arrayGB = numChanPrealloc * length(eventTrigger) * length(waveletFreqs) * DurationMS * 8 / 2^30; % 8 bytes per double, 2^30 bytes/GB
powerMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
powerMatZ = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
phaseMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);


%-- Loop through the channels: extract, filter, processes, and save...
for iChan = 1:numChannels,
    
    if PROCESS_CHANNELS_SEQUENTIALLY==1,    iChanSave = 1;
    else                                    iChanSave = iChan;  end
    
    thisChan    = chanList(iChan,:);
    thisChanStr = chanStr{iChan};
    strStart    = sprintf('\n Grab %d/%d: %s', iChan, numChannels, thisChanStr );  strStart(end+1:35)=' '; %buffer length so everything lines up
    fprintf('%s', strStart);       tic;
    
    
    % load eeg data: get EEG event data based on ms ranges... manually add the buffer so multiphasevec3 gets the extended waveform
    %                  chan,   event,       signal duration+buffer, t=0 offset,    bufferMS,  bandpass filter info,  resamplingrate
    eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);
    %eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS,OffsetMS,BufferMS,preFiltFreq,preFiltType,preFiltOrder,resampledrate);  % let gete_ms add and cut the buffer... better to do manually
    
    
    % notch filter to eliminate 60 Hz noise
    fprintf(' [%.1f sec] --> notch filt', toc); tic;
    eegWaveV = buttfilt(eegWaveV,[59.5 60.5],resampledrate,'stop',1); %-filter is overkill: order 1 --> 25 dB drop (removing 5-15dB peak)
    
    % get the phase and power
    fprintf(' [%.1f sec] --> freq decomp', toc); tic;
    [rawPhase,rawPow] = multiphasevec3(waveletFreqs,eegWaveV,resampledrate,waveletWidth);
    fprintf(' [%.1f sec] --> save', toc);  tic;
    
    % remove the leading/trailing buffer
    rawPow   = rawPow(:,:,BufferMS+1:end-BufferMS);
    rawPhase = rawPhase(:,:,BufferMS+1:end-BufferMS);
    eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS);
    eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; %in seconds
    if length(eegWaveT)<size(eegWaveV,2), 
        fprintf('wave time length off'); 
        eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
    end
    
    % Do it all in a single function call... not really using the extras in getphasepow anyway (kurtosis check; downsampling; etc)
    %[rawPhase,rawPow] = getphasepow(gl(chan,:),eventList,DurationMS,OffsetMS,BufferMS,'resampledrate',resampledrate);
    
    
    % should never need this... but happens for NIH017 attentionTask
    if length(find(rawPow==0))>0,
        fprintf('\n--------------WARNING:  RAW POWER ZERO AT %d TIME POINTS...  SOMETHING FISHY?-----------------\n', length(find(rawPow==0)))
        keyboard;
        rawPow(find(rawPow==0)) = min(rawPow(find(rawPow>0))); %set to minimum observed value
    end
        
    % normalized to 1 so all can be shifted to fit on a single plot
    eegWaveMeanSub  = eegWaveV-mean(mean(eegWaveV));   %double mean and double max because multiple events from same channel should be normalized together
    eegWaveShift    = (iChanSave-1)*2 + eegWaveMeanSub./max(max(abs(eegWaveMeanSub)));
    eegInstPow      = abs(eegWaveMeanSub).^2;
    eegInstPowShift = (iChanSave-1)*2 + eegInstPow./max(max(abs(eegInstPow)));
    
    
    % temp indicies
    iEv = 1:length(eventTrigger);
    iT  = 1:size(eegWaveV,2);
    iF  = 1:length(waveletFreqs);
    
    
    %       x-axis of time series
    waveT = eegWaveT;
    
    %        chan,event,time
    wavesRaw(iChanSave,iEv,iT)   = eegWaveV;
    wavesMsb(iChanSave,iEv,iT)   = eegWaveMeanSub;
    wavesSft(iChanSave,iEv,iT)   = eegWaveShift;
    wavesSftG(iChanSave,iEv,iT)  = eegWaveShift;  %-set equal to eegWaveShift here... if channels processed sequentially they are the same, else it will be updated after loop
    instPow(iChanSave,iEv,iT)    = eegInstPow;
    instPowSft(iChanSave,iEv,iT) = eegInstPowShift;
    
    %        chan,event,freq,time
    powerMat(iChanSave,iEv,iF,iT) = 10*log10(rawPow);
    phaseMat(iChanSave,iEv,iF,iT) = rawPhase;
    
    
    % for each eegfile stem, z-score each channel and frequency
    fprintf(' [%.1f sec] --> z-score', toc);  tic;
    stemList = unique({eventTrigger.eegfile});
    for iStem=1:length(stemList),
        fprintf('.');
        iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));
        for iF = 1:length(waveletFreqs),
            allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
            mu = mean(allVal); stdev = std(allVal);
            powerMatZ(iChanSave,iEvStem,iF,iT) = (powerMat(iChanSave,iEvStem,iF,iT)-mu)/stdev;
            if sum(isnan(powerMatZ(iChanSave,iEvStem,iF,iT)))>0
                keyboard;
            end
        end
    end
    fprintf(' [%.1f sec]', toc); tic;
    
    
    
    %%%%% create object that stiches together all data for saving and/or analysis  %%%%%%
    %-----------------------------------------------------------------------------------%
    %- figure output parameters
    data.subj           = subj;
    data.session        = session;
    data.figFontAx      = figFontAx;
    data.FIG_OFFSET     = FIG_OFFSET;
    data.HIDE_FIGURES   = HIDE_FIGURES;
    
    %- extraction parameters (filter and wavelet info)
    data.waveletFreqs   = waveletFreqs;
    data.waveletWidth   = waveletWidth;
    data.preFiltFreq    = preFiltFreq;
    data.preFiltType    = preFiltType;
    data.preFiltOrder   = preFiltOrder;
    data.preFiltStr     = preFiltStr;
    data.freqBandAr     = freqBandAr;
    data.freqBandYticks = freqBandYticks;
    data.freqBandYtickLabels = freqBandYtickLabels;
    data.chanList       = chanList;
    data.chanStr        = chanStr;
    data.thisChanStr    = thisChanStr;  %- applies to last channel processed, important when processing is sequential
    
    %- extracted EEG: evoked waveforms
    data.waveT          = waveT;        %- time
    data.wavesRaw       = wavesRaw;     %- raw wave
    data.wavesMsb       = wavesMsb;     %- mean subtracted
    data.wavesSft       = wavesSft;     %- mean shifted so multiple traces do not overlap
    data.wavesSftG      = wavesSftG;    %- mean shifted and scaled to global max to represent absolute voltage
    data.instPow        = instPow;      %- instantaneous power
    data.instPowSft     = instPowSft;   %- power mean shifted for non-overlapping traces
    
    %- extracted EEG: power and phase
    data.powerMat       = powerMat;
    data.phaseMat       = phaseMat;
    data.powerMatZ      = powerMatZ;    %- each wavelet frequency z-scored against all extracted data
    
    %- meta events variables
    data.eventsMeta     = eventsMeta;
    data.metaYval       = metaYval;
    data.metaZeroMS     = metaZeroMS;
    data.metaYstr       = metaYstr;
    data.metaYtransition = metaYtransition;
    
    %- trigger event variables
    data.THIS_TRIGGER   = THIS_TRIGGER;
    data.eventsTrig     = eventsTrig;
    data.trigType       = trigType;
    data.trigZeroMS     = trigZeroMS;
    data.trigTypeStr       = trigTypeStr;
    data.eventsTriggerXlim = eventsTriggerXlim;
    data.eventsAveWindowMS = eventsAveWindowMS;
    %-----------------------------------------------------------------------------------%
    
    
    if  PROCESS_CHANNELS_SEQUENTIALLY==1,

        fprintf(' --> plot & stats');  tic;
        
        [statsStructEvoke, responseStruct, figEvokedSpec] = plotEvoked(data);                			fprintf('.');
        allChanStats(iChan).Evoked = statsStructEvoke;
        allChanResp(iChan)         = responseStruct;
        
        %keyboard
        
        
%         
%         
%         [statsStruct500,   figTimeFreqAve500]  			  = plotTimeFreqAve(data, slidingWindow500);    fprintf('.');
%         allChanStats(iChan).Win500 = statsStruct500;
%         statsStr                   = statsStruct500.statsStr;
%         
%         %%- Epoc analysis was used initially, but now sliding windows are the main approach... comment this out (and the figTimeFreqAveEpoc save line below)
%         %[statsStructEpoc,  figTimeFreqAveEpoc] 		  = plotTimeFreqAve(data, []);                  fprintf('.');
%         %statsStr                   = statsStructEpoc.statsStr;  %- lists the significantly different epochs (frequency x time window)
%         %allChanStats(iChan).Epocs  = statsStructEpoc;
%         
%         %%- For now only extract the short time windows for sampleStart extractions... eventually select a window and use in all extractions
%         if useSmallWinDT,
%             [statsStruct250,   figTimeFreqAve250]  			  = plotTimeFreqAve(data, slidingWindow250);    fprintf('.');
%             allChanStats(iChan).Win250 = statsStruct250;
%             
%            [statsStruct100,   figTimeFreqAve100]              = plotTimeFreqAve(data, slidingWindow100);    fprintf('.');
%             allChanStats(iChan).Win100 = statsStruct100;
%         end
%         
%         
%         fprintf(' [%.1f sec] --> save figs', toc); tic;
%         
%         if iChan==1,
%             if ismac,
%                 %set(figTimeFreqAveEpoc(1), 'position',[ 180           1        1300        1345]); % left, bottom, width, height
%                 set(figEvokedSpec,         'position',[1480           1        1080        1345]); %
%                 [statsStructEvoke, responseStruct, figEvokedSpec] = plotEvoked(data); %re-plot so re-size takes effect
%             else
%                 fprintf('\n\nSIZE FIGURES AS DESIRED AND PRESS ANY KEY TO SAVE\n');
%                 pause;
%             end
%             
%             
%             %- define filenames for ppt and pdf outputs
%             figDirRootRoot  = sprintf('%s/FIGS_temp/%s',   FIG_TEMP_ROOT, THIS_TRIGGER);
%             figDirRoot      = sprintf('%s/%s/%s/%s%s',     figDirRootRoot,subj,session,THIS_REF_TYPE,preFiltStrShort);
%             outputFileSpecs = sprintf('%s/%s_%s_%s%s_%s%s',figDirRoot,THIS_TRIGGER,'Specs',subj,sessStr,THIS_REF_TYPE,preFiltStrShort(2:end));  %- root name of ppt and pdf file (just the file-type suffix will be different)
%             outputFileStats = sprintf('%s/%s_%s_%s%s_%s%s',figDirRoot,THIS_TRIGGER,'Stats',subj,sessStr,THIS_REF_TYPE,preFiltStrShort(2:end));  %- root name of ppt and pdf file (just the file-type suffix will be different)
%             outputFileSwind = sprintf('%s/%s_%s_%s%s_%s%s',figDirRoot,THIS_TRIGGER,'StatsSlideWin',subj,sessStr,THIS_REF_TYPE,preFiltStrShort(2:end));  %- root name of ppt and pdf file (just the file-type suffix will be different)
%             while exist([outputFileSpecs '.pptx'],'file') | exist([outputFileSpecs '.pdf'],'file'), outputFileSpecs = [outputFileSpecs 'X']; end  %-prevent overwritting by creating unique file name if necessary
%             while exist([outputFileStats '.pptx'],'file') | exist([outputFileStats '.pdf'],'file'), outputFileStats = [outputFileStats 'X']; end  %-prevent overwritting by creating unique file name if necessary
%             while exist([outputFileSwind '.pptx'],'file') | exist([outputFileSwind '.pdf'],'file'), outputFileSwind = [outputFileSwind 'X']; end  %-prevent overwritting by creating unique file name if necessary
%             
%             %- save behavioral summary as first figure
%             jwFigToPPTandPDF(sprintf('%s/Specs',figDirRoot), FIG_EVENTS, 'behavior', outputFileSpecs);
%             %jwFigToPPTandPDF(sprintf('%s/Stats',figDirRoot), FIG_EVENTS, 'behavior', outputFileStats);
%         else
%             pause(0.5); %give a moment for figs to resize before exporting to powerpoint
%         end
%         
%         infoStr = sprintf('%s__%s\n%s, %s \nchannel: %s \n%s \n%s', subj, thisChanStr, subj, session, thisChanStr, preFiltStr, statsStr);
%         jwFigToPPTandPDF(sprintf('%s/Specs',figDirRoot),         [figEvokedSpec],                       infoStr, outputFileSpecs);
%         %jwFigToPPTandPDF(sprintf('%s/Stats',figDirRoot),         [figTimeFreqAveEpoc(1)],               infoStr, outputFileStats);
%         %jwFigToPPTandPDF(sprintf('%s/StatsSlideWin',figDirRoot), [figTimeFreqAve250 figTimeFreqAve500], infoStr, outputFileSwind);
%         
%         
%         %- all done... save the mat file and make a copy of the figure files
%         if iChan==numChannels,
%             if ~exist(figDirRoot,'dir'), mkdir(figDirRoot); end;  
%             save(sprintf('%s/allChanStats.mat',figDirRoot),   'allChanStats');  %- save stats for all channels
%             save(sprintf('%s/allChanResp.mat', figDirRoot),   'allChanResp');   %- save responses for all channels
%             
%             %- attempt to create a copy of the powerpoint and pdf in the root-root figure directory
%             FIG_DIR_LIST = {outputFileSpecs, outputFileStats, outputFileSwind};
%             for iFigOut = [1:3],  %[1 2 3]
%                 rootRootCopy = [figDirRootRoot '/' FIG_DIR_LIST{iFigOut}(max(find(FIG_DIR_LIST{iFigOut}=='/'))+1:end)];
%                 [successPPT, message, messageID] = copyfile([FIG_DIR_LIST{iFigOut} '.pptx'],[rootRootCopy '.pptx'],'f');
%                 [successPDF, message, messageID] = copyfile([FIG_DIR_LIST{iFigOut} '.pdf' ],[rootRootCopy '.pdf' ],'f');
%             end
%             
%         end
        
        fprintf(' [%.1f sec]', toc); tic;
        
    end %-- if PROCESS_CHAN_SEP
    
end



if PROCESS_CHANNELS_SEQUENTIALLY==0,
    %- all channels processes, so compute global mean across channels5
    globalScale = max(max(max(abs(wavesMsb))));
    for iChan=1:numChannels
        wavesSftG(iChan,iEv,iT) = (iChan-1)*2 + 2.0*wavesMsb(iChan,iEv,iT)./globalScale;
    end
    data.wavesSftG      = wavesSftG;  %- mean shifted and scaled to global max to represent absolute voltage
end
clear phase pow rawPow rawPhase eegWaveV eegWaveMeanSub eegWaveShift eegInstPow eegInstPowShift allVal; %free up these temps
fprintf('\nAll channels processed\n');



if PROCESS_CHANNELS_SEQUENTIALLY==1,
    fprintf('Channels processed sequentially and plotted one-at-a-time... all done!\n\n');
    if HIDE_FIGURES==1, close all; end
end


fprintf('\n');



