% return for a subject, every session-block's begin and end time

subj = 'NIH034';
%% LOAD EVENTS STRUCT AND SET DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
% eegRootDirHome = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

% Either go through all the sessions, or a specific session
if sessNum == -1 | length(sessNum)>1, % all sessions
    disp('STEP 1: Going through all sessions')
    session = 'Meta Session [all]';
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap');
    sessStr = '[all]';
else                                  % one session
    disp('STEP 1: Going through one session')
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap/', session);
    sessStr = sprintf('[%d]',sessNum);
end

subjDir = fullfileEEG(eegRootDir,subj); % directory to subject (e.g. NIH034)
docsDir = fullfileEEG(subjDir,'docs');  % directory to the docs (electordes.m, tagNames.txt, etc.)
talDir  = fullfileEEG(subjDir,'tal');
defaultEEGfile = fullfileEEG('/Volumes/Shares/FRNU/data/eeg/',subj,'/eeg.reref/');  % default event eegfile fields point here... switch to local before loading

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

% SPLIT INTO SESSIONS AND BLOCKS
subjSessions = unique({events.sessionName}); % e.g. sessions 0, 1, 2
subjBlocks = unique({events.blocknumber});   % e.g. blocks 0,1,2,3,4,5

if strcmp(subj, 'NIH039')
    subjSessions = subjSessions([1,2,4]);
elseif strcmp(subj, 'NIH034')
    subjSessions = subjSessions([3, 4]);
end

for iSesh=1:length(subjSessions),
    sessionFieldName = strcat('session', num2str(iSesh));
    blockTimes = zeros(length(subjBlocks),2);
    for iBlock=1:length(subjBlocks),
        sessionBlockIndices = strcmp({events.sessionName}, subjSessions(iSesh)) & ...
                                strcmp({events.blocknumber}, subjBlocks(iBlock));
        sessionBlockEvents = events(sessionBlockIndices);
        beginTime = round(min([sessionBlockEvents.mstime]));
        endTime = round(max([sessionBlockEvents.mstime]));
        
        blockFieldName = strcat('block', num2str(iBlock));
        blockTimes(iBlock,1) = beginTime;
        blockTimes(iBlock,2) = endTime;
    end
    
    subject.(sessionFieldName) = reshape(blockTimes(:) - min(blockTimes(:)), 6, 2)./1000./60;
end

sessions = fields(subject);

fig = figure;
legStr = {};
colors = {'r', 'm', 'c', 'b', 'k', 'g'};
for s=1:length(sessions)
    for i=1:6
        plot(subject.(sessions{1})(i,:), [s, s], colors{i});
        hold on
        if s==1,
            legStr{i} = strcat('block_', num2str(i-1));
        end
    end
end
legend(legStr, 'Location', 'best');
title(['Subject ', subj, ' Length of each Session-Block'])
xlabel('Time in seconds');
ylabel('Session');
ylim([0, s+1]);


