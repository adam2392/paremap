function [chanList, chanStr, numChannels, eventEEGpath] = loadChannels(docsDir, talDir, ref_type, USE_CHAN_SUBSET)
    %%- Get 1. # of channels to use, 
    %%-     2. list of channels
    %%-     3. list of channel names
    %%- select all channels, or part of the subset of channels
    jackSheet = fullfileEEG(docsDir, 'jacksheetMaster.txt');
    [chanNums chanTags] = textread(jackSheet,'%d%s%*s');  

    %%% always look at all electrodes... worry about "good" and "bad" later (bad means inter-ictal activity or seizure activity)
    %- three referencing options:  noreref (should manually subtract reference channel), reref bioploar, and reref laplacian
    chanStr = {};   % cell for all the channel names
    chanFile = 0;   % file for the channels (e.g. ~/NIH034/tal/leads.txt) 
    chanList = [];  % list of the channels (e.g. 1-96)
    iChanListSub = []; % list of the subset of channels we want to analyze (e.g. [48 1])

    switch ref_type
        case 'noreref'  
        case 'bipolar'
            fprintf('Bipolar referencing');
            chanFile      = [talDir '/leads_bp.txt'];
            [chan1 chan2] = textread(chanFile,'%d%*c%d');
            chanList      = [chan1 chan2];
            for iChan=1:size(chanList,1),
                chanStr{iChan} = sprintf('%s-%s', chanTags{find(chanNums==chan1(iChan))}, chanTags{find(chanNums==chan2(iChan))} );
            end
            chanRefs      = [];
            eventEEGpath  = '/eeg.reref/';
        case 'global' % look at global electrodes / monopolar
            fprintf('STEP 1: Using Global referencing\n');
            chanFile      = [talDir '/leads.txt'];
            chanList      = textread(chanFile,'%d'); % read in the list of channels nums

            % set the names for each channel
            for iChan=1:size(chanList,1),
                chanStr{iChan} = sprintf('%s-global', chanTags{find(chanNums==chanList(iChan))} );
            end
            eventEEGpath  = '/eeg.reref/';
        otherwise
            fprintf('Error, no referencing scheme selected');
    end
    fprintf('\n');
    
    iChanListSub  = 58:104;            %G1, G2, LF1, AST1,
    %%- select all channels, or part of the subset of channels
    if USE_CHAN_SUBSET==0,
        iChanList = 1:size(chanList,1);  %all possible channels
    else
        iChanList = iChanListSub;
    end

    % what is this doing here?
    chanListUse = [];  chanStrUse = {};
    for iChan=iChanList,
        chanListUse(end+1,:) = chanList(iChan,:);
        chanStrUse{end+1}    = chanStr{iChan};
    end

    % reset variables and create list of channels and their corresponding names
    chanList = chanListUse;
    chanStr  = chanStrUse;
    numChannels = size(chanList,1);
end