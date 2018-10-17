%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written By: Adam Li
% Contact: ali39@jhu.edu
% 
% Inputs:   patseizevent - the struct with all the seizure events for this
%                          patient.
%           patinfo      - the struct with all the patient info, including
%                          filename and time frame of each recording
%           varargin     - the optional argument to pass in (offset in 
%                          SECONDS). Is used for generating an array of 
%                          time frames +/- offset of the seizure.
%
% Output:   A - struct with two fields:
%                   1. seizefiles - a list of recording files with seizures
%                   2. seiztimes - time frame of those seizures (optional:
%                   with offset)
%                   3. rectimes - an array of the start/end times of each
%                   of the seizure files
% Ver.: 1.0 - Date: 08/13/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [seizfiles, seiztimes, rectimes] = isolate_seizures(patseizevent, patinfo, varargin)
    nVarargs = length(varargin);

    %% Check to see if any seizure events take place over 2 separate files
    numseiz = length(patseizevent.dur);     % the number of seizures
    recindices = zeros(numseiz,1);          % the row indices in patientinfo where seizures occur
    seizfiles = cell(numseiz,1);            % the file names of the recordings with seizures
    rectimes = zeros(numseiz,2);            % the recording times of the seizure files
    seiztimes = patseizevent.time;
    seiztimes = round(seiztimes);    % round to whole numbers
    % loop through each seizure event
    for i=1:numseiz
        time = patinfo.time;
        % get the row where the seizure starts in patinfo
        start = patseizevent.time(i,1);
        endseize = patseizevent.time(i,2);

        %%%% If there are no empty matrices, then all seizures happen in
        %%%% separate recording files.
        
        % find where start time is less then seizure time AND the end time
        % for that recording is greater then seizure time -> rows in
        % patinfo struct
        recindices(i) = find(start > time(1:end, 1) & endseize < time(1:end, 2));
        seizfiles{i} = patinfo.filename{recindices(i)}; % return as is
        rectimes(i,1) = patinfo.time(recindices(i),1);
        rectimes(i,2) = patinfo.time(recindices(i),2);
    end
    if (ismember(1, cellfun(@isempty,seizfiles)))
        disp('Seizures are located in separate files...');
    else
        disp('Seizure events located in same rec. files');
    end
    
    % set times with seizure inside +/ some offset
    if (nVarargs == 0)
        seiztimes = patseizevent.time
    elseif (nVarargs == 1) % need to use offset
        seiztimes(:,1) = patseizevent.time(:,1) - varargin{1};
        seiztimes(:,2) = patseizevent.time(:,2) + varargin{1};
    else
        seiztimes = patseizevent.time;
        disp('There were too many inputs into isolate_seizures.m!');
    end
    seiztimes = round(seiztimes);
%     seiztimes(:,1) = floor(seiztimes(:,1));
%     seiztimes(:,2) = ceil(seiztimes(:,2));
end

