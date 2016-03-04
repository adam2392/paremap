%%- change for different type of expeirments
% tailored for checking paremap events

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

% Either go through all the sessions, or a specific session
if sessNum == -1 | length(sessNum)>1, % all sessions
    disp('STEP 1: Going through all sessions')
    session = 'Meta Session [all]';
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/');
    sessStr = '[all]';
else                                  % one session
    disp('STEP 1: Going through one session')
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/', session);
    sessStr = sprintf('[%d]',sessNum);
end

subjDir = fullfileEEG(eegRootDir,subj); % directory to subject (e.g. NIH034)

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------------------------------       plot meta event data        -------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_EVENTS=1;
if ~exist('FIG_OFFSET','var'), FIG_OFFSET = 0; end %- default to 0, but if called as function then allow that value to persist
FIG_EVENTS = 100 + FIG_OFFSET;

% Plot histogram of each field
if (PLOT_EVENTS)
    figure(FIG_EVENTS); clf % open figure for plotting
    set(gcf, 'color', 'w')  % set color background to white
   
    %% EXTRACT DATA FROM EVENTS
    eTimeOn = [events.mstime];
    fixationOn = ([events.fixationOnTime] - eTimeOn);     % fixation on time
    fixationOff = ([events.fixationOffTime] - eTimeOn);   % fixation off time
    matchOnTime = ([events.matchOnTime] - eTimeOn);
    probeOffTime = ([events.probeOffTime] - eTimeOn);
%     responseTime = ([events.responseTime] - eTimeOn);
    
    eTimeOnS = (eTimeOn-eTimeOn)/1000;           % time on (seconds)
    fixationOnS = fixationOn/1000;     % fixation on time
    fixationOffS = fixationOff/1000;   % fixation off time
    matchOnTimeS = matchOnTime/1000;   
    probeOffTimeS = probeOffTime/1000;
    responseTimeS = [events.responseTime]/1000;
    
    %%- Create Histograms of events metadata
    labeloffset = 50;
    
    % extract the information from a histogram function
    [n1, xout1] = hist(fixationOnS);
    [n2, xout2] = hist(eTimeOnS);
    [n3, xout3] = hist(fixationOffS);
    [n4, xout4] = hist(matchOnTimeS);
    [n5, xout5] = hist(probeOffTimeS);
    [n6, xout6] = hist(responseTimeS);
    
    %% BAR CHARTS
    %%-01 plot fixation on/off times - should be < 0
    subplot(3,1,1)
%     axSpec(1)=gca;
    bar(xout1,n1,'r'); 
    grid; hold on;
    bar(xout3, n3, 'b');
    
    set(gca, 'XLim', [min(xout1)-0.05, max(xout3)+0.05])
    % add text
    text(mean(fixationOnS),max(n1)+labeloffset, 'fixation On')
    text(mean(xout3), max(n3)+labeloffset, 'fixation off')
    
    % add legend
    legend('fixation On', 'fixation off', 'Location', 'northeast')
    
    %%-02 plot the rest on/off times - should be < 0
    subplot(3,1,2)
%     axSpec(2)=gca;
    bar(xout2,n2,'g'); grid; hold on;
    bar(xout4, n4, 'w');
    bar(xout5, n5, 'r');

    set(gca, 'XLim', [0-0.05, max([xout4, xout5])+0.05])

    text(0, max(n2)+2*labeloffset,'probeWordOn')
    text(mean(xout4), max(n4)+labeloffset, 'match came on')
    text(mean(xout5), max(n5)+labeloffset*3, 'match/probe off')
    
    legend('probewordOn', 'match came on', 'match/probe off', 'response Time')
    
    %%-03 response times
    subplot(3,1,3)
%     axSpec(3)=gca;
    bar(xout6, n6, 'b'); grid;
    text(mean(xout6), max(n6)+labeloffset*6, 'response Time')
    
    %%% add a title
    subplot(3,1,1)
    title(sprintf('Plotting Events meta data for %s with %s events', subj, num2str(length(events))))
end