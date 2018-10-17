% subj = 'NIH034';
% typeTransform = 'morlet';
% timeLock = 'vocalization';
% referenceType = 'bipolar';
% typeReinstatement = 'across_blocks';

subj = 'NIH034';
% rankCosineSimilarity(subj, 'morlet', 'vocalization', 'bipolar', 'across_blocks');
% rankCosineSimilarity(subj, 'morlet', 'vocalization', 'bipolar', 'within_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization', 'bipolar', 'across_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization', 'bipolar', 'within_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization', 'global', 'across_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization', 'global', 'within_blocks');
rankVocalizationCosineSim(subj, 'morlet', 'vocalizationWord', 'bipolar', 'within_blocks');
% rankVocalizationCosineSim(subj, 'morlet', 'vocalizationWord', 'bipolar', 'across_blocks');
rankVocalizationCosineSim(subj, 'multitaper', 'vocalizationWord', 'bipolar', 'within_blocks');
% rankVocalizationCosineSim(subj, 'multitaper', 'vocalizationWord', 'bipolar', 'across_blocks');

subj = 'NIH039';
% rankCosineSimilarity(subj, 'morlet', 'vocalization', 'bipolar', 'across_blocks');
% rankCosineSimilarity(subj, 'morlet', 'vocalization_allPairs', 'bipolar', 'within_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization', 'bipolar', 'across_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization_allPairs', 'bipolar', 'within_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization', 'global', 'across_blocks');
% rankCosineSimilarity(subj, 'multitaper', 'vocalization', 'global', 'within_blocks');
rankVocalizationCosineSim(subj, 'morlet', 'vocalizationWord', 'bipolar', 'within_blocks');
% rankVocalizationCosineSim(subj, 'morlet', 'vocalizationWord', 'bipolar', 'across_blocks');
rankVocalizationCosineSim(subj, 'multitaper', 'vocalizationWord', 'bipolar', 'within_blocks');
% rankVocalizationCosineSim(subj, 'multitaper', 'vocalizationWord', 'bipolar', 'across_blocks');

subj = 'NIH037';
rankVocalizationCosineSim(subj, 'morlet', 'vocalizationWord', 'bipolar', 'within_blocks');
% rankVocalizationCosineSim(subj, 'morlet', 'vocalizationWord', 'bipolar', 'across_blocks');
rankVocalizationCosineSim(subj, 'multitaper', 'vocalizationWord', 'bipolar', 'within_blocks');
% rankVocalizationCosineSim(subj, 'multitaper', 'vocalizationWord', 'bipolar', 'across_blocks');
