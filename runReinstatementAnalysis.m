clear
close all
clc

subj_list = {'NIH034', 'NIH039'};
vocalization = [0, 1];
%%- run within blocks analysis
% createWithinBlocksReinstatementMat('NIH034', 0); % sync to probewordon
% createWithinBlocksReinstatementMat('NIH039', 0); % sync to probewordon
% createWithinBlocksReinstatementMat('NIH034', 1); % sync to vocalization
% createWithinBlocksReinstatementMat('NIH039', 1); 

%%- run across blocks analysis
createAcrossBlocksReinstatementMat('NIH034', 0); % sync probewordon
createAcrossBlocksReinstatementMat('NIH039', 0); % sync probewordon
createAcrossBlocksReinstatementMat('NIH034', 1);
createAcrossBlocksReinstatementMat('NIH039', 1);

% parfor i=1:length(subj_list)
%     createWithinBlocksReinstatementMat(subj_list{i}, 0);
%     createAcrossBlocksReinstatementMat(subj_list{i}, 0);
%     createWithinBlocksReinstatementMat(subj_list{i}, 1);
%     createAcrossBlocksReinstatementMat(subj_list{i}, 1);
% end