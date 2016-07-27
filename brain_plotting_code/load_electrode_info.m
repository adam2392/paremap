% FUNCTION: LOAD_ELECTRODE_INFO
% Description: To load in the electrode information for either monopolar,
% or bipolar electrodes
% Input:
% - pat_s = patient ID (e.g. NIH039)
% - bp_flag = 0, or 1 for monopolar, or bipolar respectively
% Output:
% - els = struct containing metadata about each electrode
function els = load_electrode_info(pat_s,bp_flag)
    % check each possible directory for the subject information
    eegRootDirWork = fullfile('/Users/liaj/Documents/MATLAB/paremap/', pat_s);     % work
    eegRootDirHome = fullfile('/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/', pat_s);  % home
    eegRootDirVolume = fullfile('/Volumes/NIL_PASS/', pat_s);
    eegRootDirJhu = fullfile('/home/adamli/paremap/', pat_s);

    % Determine which directory we're working with automatically
    if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
    elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
    elseif ~isempty(dir(eegRootDirVolume)), eegRootDir = eegRootDirVolume;
    elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
    else   error('Neither Work nor Home EEG directories exist! Exiting'); end

    %Load a struct containing the electrode information for the bipolar
    %electrodes
    if exist('bp_flag','var')
        if bp_flag==0
            %e = load(['/Volumes/Shares/FRNU/data/eeg/tal/monopolar/talOrig/allTalLocs_GM_' pat_s '.mat']);
            e = load([eegRootDir, '/tal/talSurf_monopolar.mat']);
            els = e.events;
        else
            e = load([eegRootDir, '/tal/talSurf_bipolar.mat']);
            els = e.events;
        end
    end
end