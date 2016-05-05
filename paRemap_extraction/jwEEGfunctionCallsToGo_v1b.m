%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   JW's calls to eeg_toolbox alignment functions
%   modified: AL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%-----------------------------------------------------------------------------------------------------%
%----------------------  Choose the subject of interest and root EEG directory -----------------------%
<<<<<<< HEAD
subj       = 'NIH037';
=======
subj       = 'NIH039';
>>>>>>> dbf44b78f069cad9758117ab7a41149dced9f341


rootEEGdirWork  = '/Volumes/JW24TB/data24TB/eeg';                        %office-local
%rootEEGdirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg';   %office-local
rootEEGdirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';   %home
<<<<<<< HEAD
rootEEGdirHome = '/home/adamli/paremap';
=======
% rootEEGdirHome = '/Volumes/NIL_PASS';
>>>>>>> dbf44b78f069cad9758117ab7a41149dced9f341
if      exist(rootEEGdirWork, 'dir'), rootEEGdir = rootEEGdirWork; 
elseif  exist(rootEEGdirHome, 'dir'), rootEEGdir = rootEEGdirHome; 
else    fprintf('WARNING: cant find root directory (not home or work'); end

% rootEEGdir = '/Volumes/Shares/FRNU/data/eeg';                      %office-server
%rootEEGdir = '/Volumes/Shares/FRNU/dataWorking/eeg';                %office-server

%-----------------------------------------------------------------------------------------------------%
%----------------------  Extract Behavioral Events for the Attention Task ----------------------------%
taskList = {'attentionTask' 'palRam' 'palRamStim' 'paRemap' 'languageTask' 'moveTask' 'pa3' 'paRepeat'  'playPass' 'stimMapping'  'auditoryLexicalDecision', 'auditoryVowels', 'goAntiGo', 'SequenceMem'};

% call behavioral processing for certain tasks within taskList
for iTask = [4],
    %behavioralProcessing(subj, rootEEGdir, taskList{iTask}, 'behavioral_preOp');   % extract preOp tasks
    %behavioralProcessing(subj, rootEEGdir, taskList{iTask}, 'behavioral_postOp');  % extract postOp tasks
    behavioralProcessing(subj, rootEEGdir, taskList{iTask}, 'behavioral');         % extract iEEG tasks
end%

%-----------------------------------------------------------------------------------------------------%
%---------------------- Prep for and Align; Create alignmentSummary.txt  -----------------------------%
eegPrepAndAlign(subj, rootEEGdir);%
