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
%                   1. iitimes - time frame of interictal time periods
%                   (optional: with offset)
%                   2. rectimes - an array of the start/end times of each
%                   of the seizure files
% Ver.: 1.0 - Date: 08/13/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iitimes] = isolate_interictal(patseizevent, patinfo, varargin)
    nVarargs = length(varargin);
    %% Check to see if any seizure events take place over 2 separate files
    numseiz = length(patseizevent.dur);     % the number of seizures
%     recindices = zeros(numseiz,1);          % the row indices in patientinfo where seizures occur
%     rectimes = zeros(numseiz,2);            % the recording times of the interictal files
    iitimes = zeros(numseiz+1, 2);          % the interictal time frames
    
    time = patinfo.time;
    startrec = time(1,1);           % store beginning of the recording session
    endrec = time(end,2);           % store ending of the recording session
    
    % loop through each seizure event and locate it
    for i=1:numseiz
        % get the row where the seizure starts in patinfo
        startseize = patseizevent.time(i,1);
        endseize = patseizevent.time(i,2);

        % store the interictal start, end times 
        iitimes(i,1) = startrec;
        iitimes(i,2) = startseize;
        
        % set times with seizure inside +/ some offset
        if (nVarargs == 0)
            % do nothing
        elseif (nVarargs == 1) % need to use offset
            if i==1                       % if first rec session
               iitimes(i,2) = iitimes(i,2) - varargin{1};   % set ending of first ii period back
            elseif i ~= 1 && i ~= numseiz % if not first or last recording session
               iitimes(i,1) = iitimes(i,1) - varargin{1}; 
               iitimes(i,2) = iitimes(i,2) + varargin{1};
            elseif i == numseiz           % if last rec session
               iitimes(i,1) = iitimes(i,1) + varargin{1};
            end
        else % too many inputs
            % 1) Construct an MException object to represent the error.
            msgID = 'MYFUN:BadIndex';
            msg = 'Too many inputs.';
            baseException = MException(msgID,msg);
            throw(baseException);
        end
        
        % reset startrec and endrec for a loop
        startrec = endseize;
    end
    iitimes(end, 1) = startrec;
    iitimes(end, 2) = endrec;
    
    iitimes = round(iitimes);
end