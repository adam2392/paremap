clear
close all
clc
addpath('./m_reinstatement/');
subj_list = {'NIH034', 'NIH039'};
OPTIONS = [0, 1]; % PROBEWORDON, VOCALIZATION, MATCHWORD

%% Looking at Pairs of Words
%%- MORLET BIPOLAR VOCALIZATION within/across blocks
% createWithinBlocksReinstatementMat('NIH034', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH034', 'morlet', 'vocalization', 'bipolar');


% VOCALIZED WORD GROUPS
% createWithinBlocksVocalizedGroupReinstatement('NIH034', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH034', 'morlet', 'vocalization', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH039', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH039', 'morlet', 'vocalization', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH039', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH039', 'multitaper', 'vocalization', 'bipolar');



%%- MULTITAPER GLOBAL VOCALIZATION
% createWithinBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createWithinBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'global');
% % createWithinBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'global');
% createAcrossBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'global');
% createAcrossBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'global');

% NIH039
createWithinBlocksReinstatementMat('NIH039', 'morlet', 'vocalization', 'bipolar');
createAcrossBlocksReinstatementMat('NIH039', 'morlet', 'vocalization', 'bipolar');
createWithinBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'bipolar');
createAcrossBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'bipolar');
