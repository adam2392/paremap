clear all
clc
%% SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];

 % set this to the session's of each subject
 % Should be 'session_<sessions>'
sessions = {'0a', '0b', '1', '2'};

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
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap/');
    sessStr = '[all]';
else                                  % one session
    disp('STEP 1: Going through one session')
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap/', session);
    sessStr = sprintf('[%d]',sessNum);
end

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load Annotations and Read Data  ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n');
disp('STEP 2: Getting Data From Annotated Files');

%%- Reading the annotations for each session
for i=1:length(sessions)
    % get list of annotated files
    sessionDir = sprintf('%s/session_%s', behDir, sessions{i}); 
    list_of_annotations = dir(fullfile(sessionDir, '*.ann'));
    list_of_annotations = {list_of_annotations.name}';
    
    % debug print statements
    fprintf('On this session: %s', sessions{i});
    fprintf('\n');
    
    % initialize 1) session data array to hold data for each session's
    % annotations and 2) length of each annotation to loop cell array
    sessionDataArray = {};
    sessionDataArray{1} = [];
    sessionDataArray{2} = [];
    sessionDataArray{3} = [];
    lengthOfAnnotations = [];
    for ann_index=1:length(list_of_annotations) % loop through each annotated file
        % debug print statements
%         fprintf('\n');
%         fprintf('On this annotated file: %s', list_of_annotations{ann_index});
        
        % get the full file location
        ann_file = fullfile(sessionDir, list_of_annotations{ann_index});
        
        % Open the text file.
        fileID = fopen(ann_file, 'r');     
        tline = fgetl(fileID);
        msoffset = [];
        typeResponse = [];
        responseWord = {};
        
        %%- loop through entire annotation file
        while ischar(tline) 
            tline = fgetl(fileID); % get new line of the file
            
            if(tline ~= -1)        % if new line is valid
                if(tline(1) ~= '#')% beginning of data lines
%                     disp(tline)
                    % create annotated cell data split by delimiter
                    annotation_data = strsplit(tline);
                    
                    % create array of the annotated data
                    msoffset = [msoffset; str2num(annotation_data{1})];            % get the time (ms) after recording starts
                    typeResponse = [typeResponse; str2num(annotation_data{2})];        % type of response (-1, 1, 2, 3, 4, 5)
                    responseWord = [responseWord; annotation_data{3}];        % the actual word response
                end
            end
        end
        
        % each column is another annotation file's results
        sessionDataArray{1} = [sessionDataArray{1}; msoffset];
        sessionDataArray{2} = [sessionDataArray{2}; typeResponse];
        sessionDataArray{3} = [sessionDataArray{3}; responseWord];
        lengthOfAnnotations = [lengthOfAnnotations; length(msoffset)];
        
        % Close the text file.
        fclose(fileID);
    end
    break
end

fprintf('\n')
disp('sessionDataArray and lengthOfAnnotations are the resulting data structs');
fprintf('For session # %s, Length of annotations were: %d', sessions{i}, sum(lengthOfAnnotations));
fprintf('\n')

%%- read in session.log file
delimiter = '\t';
startRow = 4;
formatSpec = '%s%s%s%s%s%[^\n\r]';
filename = fullfile(sessionDir, 'session.log')

%%- Open the sessionlog and read it in
fileID = fopen(filename,'r');
sessionLogArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

