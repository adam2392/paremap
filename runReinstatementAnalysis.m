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
% %%- MULTITAPER GLOBAL VOCALIZATION
% createWithinBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH034', 'multitaper', 'vocalization', 'bipolar');


% VOCALIZED WORD GROUPS
% createWithinBlocksVocalizedGroupReinstatement('NIH034', 'morlet', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH034', 'morlet', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH034', 'multitaper', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH034', 'multitaper', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH039', 'morlet', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH039', 'morlet', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH039', 'multitaper', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH039', 'multitaper', 'bipolar');
% 
% % NIH039
% createWithinBlocksReinstatementMat('NIH039', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH039', 'morlet', 'vocalization', 'bipolar');
% createWithinBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH039', 'multitaper', 'vocalization', 'bipolar');

plotVocalizedReinstatement('NIH034', 'morlet', 'bipolar', 'within_blocks_vocalizationWord');
plotVocalizedReinstatement('NIH034', 'morlet', 'bipolar', 'across_blocks_vocalizationWord');
plotVocalizedReinstatement('NIH034', 'multitaper', 'bipolar', 'within_blocks_vocalizationWord');
plotVocalizedReinstatement('NIH034', 'multitaper', 'bipolar', 'across_blocks_vocalizationWord');
plotVocalizedReinstatement('NIH039', 'morlet', 'bipolar', 'within_blocks_vocalizationWord');
plotVocalizedReinstatement('NIH039', 'morlet', 'bipolar', 'across_blocks_vocalizationWord');
plotVocalizedReinstatement('NIH039', 'multitaper', 'bipolar', 'within_blocks_vocalizationWord');
plotVocalizedReinstatement('NIH039', 'multitaper', 'bipolar', 'across_blocks_vocalizationWord');
% plotVocalizedReinstatement('NIH037', 'morlet', 'bipolar', 'within_blocks_vocalizationWord');
% plotVocalizedReinstatement('NIH037', 'morlet', 'bipolar', 'across_blocks_vocalizationWord');
% plotVocalizedReinstatement('NIH037', 'multitaper', 'bipolar', 'within_blocks_vocalizationWord');
% plotVocalizedReinstatement('NIH037', 'multitaper', 'bipolar', 'across_blocks_vocalizationWord');

% setupSessionBlocksSpect('NIH034', 'morlet', 'bipolar', 'within_blocks');
% setupSessionBlocksSpect('NIH034', 'multitaper', 'bipolar', 'within_blocks');
% 
% setupSessionBlocksSpect('NIH034', 'morlet', 'bipolar', 'across_blocks');
% setupSessionBlocksSpect('NIH034', 'multitaper', 'bipolar', 'across_blocks');
% 
% setupSessionBlocksSpect('NIH039', 'morlet', 'bipolar', 'within_blocks');
% setupSessionBlocksSpect('NIH039', 'multitaper', 'bipolar', 'within_blocks');
% 
% setupSessionBlocksSpect('NIH039', 'morlet', 'bipolar', 'across_blocks');
% setupSessionBlocksSpect('NIH039', 'multitaper', 'bipolar', 'across_blocks');
% 
% % NIH037
% createWithinBlocksReinstatementMat('NIH037', 'morlet', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH037', 'morlet', 'vocalization', 'bipolar');
% createWithinBlocksReinstatementMat('NIH037', 'multitaper', 'vocalization', 'bipolar');
% createAcrossBlocksReinstatementMat('NIH037', 'multitaper', 'vocalization', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH037', 'morlet', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH037', 'morlet', 'bipolar');
% createWithinBlocksVocalizedGroupReinstatement('NIH037', 'multitaper', 'bipolar');
% createAcrossBlocksVocalizedGroupReinstatement('NIH037', 'multitaper', 'bipolar');