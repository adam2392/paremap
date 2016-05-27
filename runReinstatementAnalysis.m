clear
close all
clc

subj_list = {'NIH034', 'NIH039'};
OPTIONS = [0, 1]; % PROBEWORDON, VOCALIZATION, MATCHWORD

%%- run within blocks analysis
% createWithinBlocksReinstatementMat('NIH034', 0, 1); % sync to probewordon
% createWithinBlocksReinstatementMat('NIH039', 0, 1); % sync to probewordon
% createWithinBlocksReinstatementMat('NIH034', 0, 0); % sync to probewordon
% createWithinBlocksReinstatementMat('NIH039', 0, 0); % sync to probewordon
% createWithinBlocksReinstatementMat('NIH034', 1, 0); % sync to vocalization
% createWithinBlocksReinstatementMat('NIH039', 1, 0); 

% %%- run across blocks analysis
% createAcrossBlocksReinstatementMat('NIH034', 0, 1); % sync probewordon
% createAcrossBlocksReinstatementMat('NIH039', 0, 1); % sync probewordon
% createAcrossBlocksReinstatementMat('NIH034', 0, 0); % sync to probewordon
% createAcrossBlocksReinstatementMat('NIH039', 0, 0); % sync to probewordon
% createAcrossBlocksReinstatementMat('NIH034', 1, 0);
% createAcrossBlocksReinstatementMat('NIH039', 1, 0);