%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
clc;
%% SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/');

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

%%- GET CORRECT EVENTS ONLY
% POST MODIFY EVENTS based on fields we want (e.g. is it correct or not)?
correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);

%% EXTRACT DATA FROM EVENTS
%%%%%%%%%%%%%%%%%%%%%RECENTER THE TIMES WRT MSTIME %%%%%%%%%%%%%%%%
eTimeOn = [events.mstime];
fixationOn = ([events.fixationOnTime] - eTimeOn);     % fixation on time
fixationOff = ([events.fixationOffTime] - eTimeOn);   % fixation off time
matchOnTime = ([events.matchOnTime] - eTimeOn);
probeOffTime = ([events.probeOffTime] - eTimeOn);

% CONVERT ALL TIMES TO SECONDS
eTimeOnS = (eTimeOn-eTimeOn)/1000;           % time on (seconds)
fixationOnS = fixationOn/1000;     % fixation on time
fixationOffS = fixationOff/1000;   % fixation off time
matchOnTimeS = matchOnTime/1000;   
probeOffTimeS = probeOffTime/1000;
responseTimeS = [events.responseTime]/1000;

sessionNumbers = [events.sessionNum];
sessionsHad = unique([events.sessionNum]);

for i=1:length(unique([events.sessionNum])),
    index = sessionsHad(i);
    
    %%- Get events for this session
    session_events = events(sessionNumbers==index);
    
    %%- get response times
    session_responses = [session_events.responseTime]/1000;
    
    %%- plot
    figure;
    hist(session_responses);
    set(gca, 'XLim', [0, 4])
    title(sprintf('Response Time for session %s with %s events', num2str(index), num2str(length(session_events))))
    xlabel('Response Time (seconds)')
    ylabel('Frequency')
end
    
    
    
    