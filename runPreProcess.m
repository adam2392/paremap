%%- vocalization
% morlet wavelet preprocessed
ALParemap_eCog_PreProcess('NIH034', 'morlet', 'vocalization', 'bipolar', 500, 100);
ALParemap_eCog_PreProcess('NIH034', 'morlet', 'vocalization', 'global', 500, 100);
ALParemap_eCog_PreProcess('NIH039', 'morlet', 'vocalization', 'global', 500, 100);
ALParemap_eCog_PreProcess('NIH039', 'morlet', 'vocalization', 'bipolar', 500, 100);

% multitaper FFT preprocessed
ALParemap_eCog_PreProcess('NIH034', 'multitaper', 'vocalization', 'global', 500, 100);
ALParemap_eCog_PreProcess('NIH034', 'multitaper', 'vocalization', 'bipolar', 500, 100);
% ALParemap_eCog_PreProcess('NIH039', 'multitaper', 'vocalization', 'global', 500, 100);
% ALParemap_eCog_PreProcess('NIH039', 'multitaper', 'vocalization', 'bipolar', 500, 100);
