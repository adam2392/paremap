function [events] = paRemap_ExtractEvents(sessLogFile, subject, sessionName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from paRemap %%%%
%
%
%   extraction designed for paRemap v4.2, which is implemented in pyEPL and was used for NIH030 and beyond
%           (earlier versions were implemented in psychoPy) and were used for NIH028 and 029... those session logs will need some tweaking to use with this extraction
%
%   training section NOT saved to session log... for earlier versions this MUST BE DELETED from the session log
%
%
%
%%%%%   create an event for every presented word and every response (no words from the training sections should be included)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% uncomment following lines to directly run script
% clear all
%
% rootEEGdir  = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg';
% rootEEGdir  = '/Volumes/Shares/FRNU/dataWorking/eeg';
% subject     = 'NIH031';   % EEG002  NIH016
% sessionName = 'session_1';
%
% sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/paRemap',sessionName);
% sessLogFile = fullfileEEG(sessionDir,'session.log');
% eventFile   = fullfileEEG(sessionDir,'events.mat');
% priorEvents = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['\nOn session: ' sessionName])

%- For NIH028 and NIH029, copy the session log and add the correct msoffset
sessFolderPath  = sessLogFile(1:strfind(sessLogFile,sessionName)+length(sessionName));
paRemap2Sesslog = [sessFolderPath 'paRemap2_session.log'];
    
if exist(paRemap2Sesslog,'file'), % ~exist(sessLogFile,'file')         
    
    dateStrPsycho = sessionName(strfind(sessionName,'2015'):end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert psycho date/time into pyEPL date/time which comes from javaSDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %- imports required to convert "epoch time" saved by pyepl into real date
    import java.lang.System;
    import java.text.SimpleDateFormat;
    import java.util.Date;
    javaSDF = SimpleDateFormat('MM/dd/yyyy HH:mm:ss.SS');  %this java object is used for mstime to date conversion
    
    %- grab the date from the excell file
    dateNumPsycho  = datenum(dateStrPsycho, 'yyyy_mmm_dd_HHMM');
    dateStrPsycho2 = datestr(dateNumPsycho, 'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
    
    %- convert matlab datenum into milisecond number used by javaSDF
    % javaSDF.format(Date(0))              --> '12/31/1969 19:00:00.00'
    % javaSDF.format(Date(60000))          --> '12/31/1969 19:01:00.00'     % 60000 is increment of 1 minute in javatime (javatime is in miliseconds, starting at 12/31/1969  1 min = 60000 microsec)
    % javaSDF.format(Date(1424183683378))  --> '02/17/2015 09:34:43.378'    % example mstime from pyepl session.log
    dateNum0java   = datenum(char(cell(javaSDF.format(Date(0)))));     % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
    dayJava        = 24 * 60 * 60 * 1000;                                 % number of miliseconds in a day
    dayMatlab      = datenum('01/02/01 1:00')-datenum('01/01/01 1:00');   % number of days in a matlab date num (should be exactly 1)
    daysToAdd      = (dateNumPsycho-dateNum0java)/dayMatlab;
    msStartPyEPL   = round(dayJava*daysToAdd);
    dateNumMSstart = datenum(char(cell(javaSDF.format(Date(msStartPyEPL)))));  % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
    dateStrMSstart = datestr(dateNumMSstart, 'mm/dd/yy HH:MM PM');  % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
    %fprintf('\n >> converting to pyepl time reference: msoffset %d = %s  (should match psychoPy xls date = %s) << ', msStartPyEPL, dateStrMSstart, dateStrPsycho2);  %- uncomment to confirm date conversion is working
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %- copy the psychopy session log with msoffset shifted based on folder name/date
    fidRead  = fopen(paRemap2Sesslog,'r');
    fidWrite = fopen(sessLogFile,'w+'); 
    while true
        thisLine = fgetl(fidRead);
        if ~ischar(thisLine); break; end
        [msTime,pos] = textscan(thisLine,'%f',1);
        fprintf(fidWrite,'%d%s\n',msTime{1}+msStartPyEPL,thisLine(pos:end));
    end
    fclose(fidWrite);    
    fclose(fidRead);  
      
    %- copy the psychopy session log
    paRemap2eeglog = [sessFolderPath 'paRemap2_eeg.eeglog'];
    standardeeglog = [sessFolderPath 'eeg.eeglog'];

    fidRead  = fopen(paRemap2eeglog,'r'); 
    fidWrite = fopen(standardeeglog,'w+'); 
    while true
        thisLine = fgetl(fidRead);
        if ~ischar(thisLine); break; end
        [msTime,pos] = textscan(thisLine,'%f',1);
        fprintf(fidWrite,'%d%s\n',msTime{1}+msStartPyEPL,thisLine(pos:end));
    end
    fclose(fidWrite);    
    fclose(fidRead);  
    fprintf('\n New copies of session.log and eeg.eeglog were created in %s ',  sessFolderPath);
     
    %- Old version, just make the sure copy the first number copy the psychopy session log
    %[SUCCESS,MESSAGE,MESSAGEID] = copyfile(paRemap2Sesslog,sessLogFile);
    
    %[SUCCESS,MESSAGE,MESSAGEID] = copyfile(paRemap2eeglog,standardeeglog);
    
    %fprintf('\n\n Update the new copy of session.log and eeg.eeglog in %s \n --->  REPLACE FIRST TIME POINT WITH % s\n HIT return when done or break with shift F5.',  sessFolderPath, num2str(msStartPyEPL));
    %keyboard
end

fid    = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric))>1, iKeep=[1:find(diff(strNumeric)>1,1,'first')]; fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep))); strNumeric=strNumeric(iKeep); end;
sessionNum = str2num( sessionName(strNumeric) );               if isempty(sessionNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...



%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 0;
while true
    thisLine            = fgetl(fid); % get new line
    if ~ischar(thisLine); break; end  % reached EOF
    
    
    %- Generic text scan to get time, offset, and type
    xTOT                = textscan(thisLine,'%f%d%s');
    msoffset            = xTOT{2}(1);
    type                = xTOT{3}{1};
    
    
    %- default Parameters (details will be filled out/altered based on type)
    experiment          = 'paRemap';
    subject             = subject   ;
    sessionName         = sessionName ;
    sessionNum          = sessionNum  ;  %- store in state var so all events are assigned a sessionNum
    msoffset            = msoffset  ;
%     isCorrect           = nan;
    
    switch type
        case {'REC_START'} % get name of .ann file and extract the contents
            annindex = 0;   % initialize going through annotation file index
            recstartTOT=textscan(thisLine,'%f%d%s%s'); % grab the block number
            annfilename = fullfile(sessionDir, strcat(recstartTOT{4}{1}, '.ann'));
            
            disp(['On this file: ' annfilename])
            annfid = fopen(annfilename, 'r');    % open file to read
            
            recStart = recstartTOT{1}(1);  % get this block's recording start time
            
            if ~strcmp(sessionName, 'session_0b')
                %%- Error check on opening up annotation .ann file
                if (annfid==-1), error('%s.ann not found: \n Exiting.',xTOT{4}(1)); end
            end
            
            % get lines until reached annotated data
            while true  
                tempLine = fgetl(annfid); % get new line
                if isempty(tempLine), break; end  % reached EOF
            end
            clear tempTot tempLine
            
            % read rest of annotation
            while true
                line = fgetl(annfid);
                if ~ischar(line); break; end  % reached EOF
                %- Generic text scan to get time, offset, and type
                annTOT                = textscan(line,'%f%d%s');
                
                % add to list these annotation data
                if annindex == 0
                    timeAfterRec          = annTOT{1}(1) + recStart;   %- must add (1) because numbers after the string in the above line cause overflow to first %f
                    wordType              = annTOT{2}(1);
                    vocalizedWord         = {annTOT{3}{1}};
                    annindex = annindex + 1;
                else
                    timeAfterRec          = [timeAfterRec; annTOT{1}(1) + recStart];   %- must add (1) because numbers after the string in the above line cause overflow to first %f
                    wordType              = [wordType; annTOT{2}(1)];
                    vocalizedWord         = [vocalizedWord; annTOT{3}{1}];
                    annindex = annindex + 1;
                end
            end
            
            % store total number of annotations in .ann file
            TOTAL_ANN = annindex;
            
            % variable to help in finding response words/times
            annwordindex = 1; % loop through the annotated words 1:TOTAL_ANN
            probeFound = 0;
            fclose(annfid); % close annotation file
        case {'BLOCK_0', 'BLOCK_1', 'BLOCK_2', 'BLOCK_3', 'BLOCK_4', 'BLOCK_5'}
            blockTOT = textscan(thisLine, '%f%d%s%s%s'); 
            if strcmp(blockTOT{5}(1), 'TEST')
                blocknumber = type;             % store the block number
                miniblocknumber = blockTOT{4}(1);   % store miniblocknumber
            end  
        case {'FIXATION_ON'}
            fixationOnTime = xTOT{1}(1);
        case {'FIXATION_OFF'}
            fixationOffTime = xTOT{1}(1);
        case {'PROBEWORD_ON'}
            isProbe    = 1;
            probeTOT=textscan(thisLine,'%f%d%s%s%s%s%s'); % grab the block number
            probeWord   = probeTOT{4}{1};
            targetWord  = probeTOT{6}{1};
            mstime = probeTOT{1}(1);
            type = type;
            %thisAnnFile = xTOT{7}{1};  %- one version has an annotation file associated with each word, eventually that was jettisoned.  
            
            probeFound = 1;
        case {'MATCHWORD_ON'}
            matchOnTime = xTOT{1}(1); % get the mstime of this line
            
            %%- When matchword comes on, should have vocalized...
            % determine timerange words can occur from 
            % probewordon -> matchword on (0-timeRange)
            timeRange = matchOnTime - mstime;                          % time Range the word can come on
            timeVocalizations = timeAfterRec(:) - mstime;              % convert all vocalizations wrt to the mstime
            validIndices = find(timeVocalizations > 0 & timeVocalizations < timeRange); % find valid indices of word by response times

            %%- found one word that was correct, no other vocalizations
            if length(validIndices) == 1, 
                if strcmp(targetWord, vocalizedWord{validIndices})
                    isCorrect = 1;

                    % LOG THE EVENT FIELDS and increment index through ann file
                    responseTime = timeVocalizations(validIndices); % responseTime
                    responseWord = vocalizedWord{validIndices};
                    annwordindex = validIndices(end);
                else % incorrect word vocalized
                    isCorrect = 0;
                    responseTime = timeVocalizations(validIndices); % responseTime
                    responseWord = vocalizedWord{validIndices};
                    
                    disp('Error in strcmp first if...');
                end
            elseif length(validIndices) == 0, %%- no word response in this frame period
                isCorrect = 0;
                responseTime = 0;
                responseWord = 'none';

                annwordindex = annwordindex + 1;
            else %%- more then 1 word vocalization found within time frame
                isCorrect = 0; % defined since they vocalized more then 1 word
                if strcmp(targetWord, vocalizedWord{validIndices(1)}) % first try was correct
                    responseTime = timeVocalizations(validIndices(1));
                    responseWord = vocalizedWord{validIndices(1)};
                elseif strcmp(targetWord, vocalizedWord{validIndices(end)}) % last vocalization was correct
                    responseTime = timeVocalizations(validIndices(end));
                    responseWord = vocalizedWord{validIndices(end)};
                else % no vocalization was correct, or it was in the middle
                    responseTime = timeVocalizations(validIndices(1)); % get the first response word
                    responseWord = vocalizedWord{validIndices(1)}; % get the first vocalized word
                end
            end
        case {'PROBEWORD_OFF'} % SAME TIME MATCHWORD TURNS OFF
            probeOffTime = xTOT{1}(1); % get the mstime of this line
                        
            if(probeFound == 1), index = index+1; end
    end
    
    %%%%%% CURRENTLY ONLY APPENDING EVENTS WHEN INDEX ++ed by PROBEWORD ON
    %- asign values to events array
    if index>length(events),
        
        % just making sure all the times are in order
        if matchOnTime < mstime || probeOffTime < matchOnTime ...
            || fixationOnTime > probeOffTime ...
            || fixationOffTime < fixationOnTime
            disp('error in times')
        end
        
        
        %- create dummy event structure that is upddated below based on type
        clear thisEvent
        thisEvent.experiment        = experiment  ;
        thisEvent.subject           = subject     ;
        thisEvent.sessionName       = sessionName ;
        thisEvent.sessionNum        = sessionNum  ;   % store in state var so all events are assigned a sessionNum  %% JW updated 2/2015
        thisEvent.type              = type        ;
        thisEvent.msoffset          = msoffset    ;
        thisEvent.mstime            = mstime      ;
        
        %- event identity
%         thisEvent.isProbe           = isProbe        ;   %- 1 or 0
%         thisEvent.isResponce        = isResponse     ;   %- 1 or 0
        
        thisEvent.probeWord         = probeWord     ;
        thisEvent.targetWord        = targetWord    ;
        thisEvent.isCorrect         = isCorrect     ;
        thisEvent.blocknumber       = blocknumber   ;
        thisEvent.miniblocknumber   = miniblocknumber;
        thisEvent.matchOnTime       = matchOnTime   ;
        thisEvent.probeOffTime      = probeOffTime  ;
        thisEvent.fixationOnTime    = fixationOnTime;
        thisEvent.fixationOffTime   = fixationOffTime;
        thisEvent.responseTime      = responseTime;
        thisEvent.vocalizedWord     = responseWord;
        
        if (index==1)
            events        = thisEvent; %- before events defined must convert to structure
        else
            events(index) = thisEvent;
        end
    end
    
end
fclose(fid);  % close session.log