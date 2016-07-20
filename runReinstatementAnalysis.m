clear
close all
clc
addpath('./m_reinstatement/');
subj_list = {'NIH034', 'NIH039'};
OPTIONS = [0, 1]; % PROBEWORDON, VOCALIZATION, MATCHWORD

%% Looking at Pairs of Words
% NIH034
%%- MORLET BIPOLAR VOCALIZATION within/across blocks
% createWithinBlocksReinstatementMat('NIH034', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH034', 'morlet', 'vocalization', 'bipolar');
%%- MULTITAPER GLOBAL VOCALIZATION
% createWithinBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'bipolar');


% VOCALIZED WORD GROUPS
% createWithinBlocksVocalizedGroupReinstatement('NIH034', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH034', 'morlet', 'vocalization', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH039', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH039', 'morlet', 'vocalization', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH039', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH039', 'multitaper', 'vocalization', 'bipolar');

% NIH039
createWithinBlocksReinstatementMat('NIH039', 'morlet', 'vocalization', 'bipolar');
createAcrossBlocksReinstatementMat('NIH039', 'morlet', 'vocalization', 'bipolar');
createWithinBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'bipolar');
createAcrossBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'bipolar');
