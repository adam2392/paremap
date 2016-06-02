% function createAcrossBlocksVocalizedGroupReinstatement(subj)
    close all;
    
    subj = 'NIH034';
    sessNum = [0, 1, 2];
    addpath('./m_reinstatement/');
    
    %% INITIALIZATION OF SESSION AND BLOCKS TO LOOK AT
    TYPE_TRANSF = 'morlet_spec_vocalization';
    disp('WITHIN BLOCKS');
    disp(TYPE_TRANSF);
    
    dataDir = strcat('./condensed_data_', subj);
    dataDir = fullfile(dataDir, TYPE_TRANSF);
    sessions = dir(dataDir);
    sessions = {sessions(3:end).name};
    % sessions = sessions(3:end);

    if strcmp(subj, 'NIH039')
        sessions = sessions([1,2,4]);
    elseif strcmp(subj, 'NIH034')
        sessions = sessions([3, 4]);
    end
    sessions
    blocks = dir(fullfile(dataDir, sessions{1}));
    blocks = {blocks(3:end).name};
    
    %% RUN ANALYSIS
    %%- LOOP THROUGH SESSIONS AND BLOCKS
    for iSesh=1:length(sessions),
        for iBlock=1:length(blocks),
            %%- 01: BUILD WORD PAIRS
            % get word pairs in this session-block
            targetWordsOne = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
            targetWordsOne = {targetWordsOne(3:end).name};
            targetWordsTwo = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock+1}));
            targetWordsTwo = {targetWordsTwo(3:end).name};
            
            test = {};
            for i=1:4
                test{i,1} = targetWordsOne{i};
                test{i,2} = targetWordsTwo{i};
            end
            test
            combnk(test, 2)
            wordPairs = createAcrossVocalizedWordGroups(targetWordsOne, targetWordsTwo)
            
        end
    end
% end