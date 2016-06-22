function [fx, PS, PSD, freqs_FFT, all_phase, t_sec] = fft_spectral(X,Fs,T,overlap,varargin)
%
%FFT_SPECTRAL
%Written by: Julio I. Chapeton
%
%Calculates the power spectrum (PS) and power spectral density (PSD) of a signal X.
%By default the PSD is calculated by using STFTs with each time segment being tapered with a rectangular window (Welch's method). Alternatively different
%window types can be specified (e.g. Hann, Bartlett,etc...), or the multitaper option can be specified to use Thomson's multitaper PSD estimate.
%The PS is calculated by averaging the STFTs over time.
%For most cases,PS is scaled such that the height of the PS corresponds to the amplitude of the signal in the time domain, and the PSD is scaled to satisfy Parseval's theorem.
%When multitapers are used the scale is somewhat arbitrary. If called with no arguments it will plot the spectrum and a spectrogram if there are multiple windows.
%
%Usage:
%FFT_SPECTRAL(X,Fs,T,overlap,TW,varargin)
%FFT_SPECTRAL(...,'PropertyName',PropertyValue,...)
%
%Inputs:
%    Required:
%        eeg        = vector (contains time series)
%        Fs         = scalar (sampling frequency in Hz)
%        T          = scalar (size of the time window in seconds)
%        overlap    = scalar (size of overlap for adjacent windows, as a fraction of T)
%    Optional:
%        (name-value pairs, not case sensitve)
%        window     = handle (defines the type of window that will be used for Welch's method, the default is rectangular. Input a handle (@blah), not a string ('@blah') )
%                   Allowed property values:
%                   @bartlett       - Bartlett window.
%                   @barthannwin    - Modified Bartlett-Hanning window.
%                   @blackman       - Blackman window.
%                   @blackmanharris - Minimum 4-term Blackman-Harris window.
%                   @bohmanwin      - Bohman window.
%                   @chebwin        - Chebyshev window.
%                   @flattopwin     - Flat Top window.
%                   @gausswin       - Gaussian window.
%                   @hamming        - Hamming window.
%                   @hann           - Hann window.
%                   @kaiser         - Kaiser window.
%                   @nuttallwin     - Nuttall defined minimum 4-term Blackman-Harris window.
%                   @parzenwin      - Parzen (de la Valle-Poussin) window.
%                   @rectwin        - Rectangular window.
%                   @taylorwin      - Taylor window.
%                   @tukeywin       - Tukey window.
%                   @triang         - Triangular window.
%       multitapers = scalar (option to use multitapers, property value should be the time-bandwidth product TW, 2*TW-1 will be the number of tapers,...
%                            and must be multiple of 1/2. If this option is used the tapers (windows) will be discrete prolate spheroidal (Slepian) sequences)
%       weights     = string (the property value should be 'eigen' to weight the result of each taper by its eigenvalue, must be using the 'multitaper' option)
%
%Outputs:
%    Required:
%        fx         = array  (complex array with the one-sided fft of the form (frequency x window x taper))
%        PS         = vector (spectrum up to Fnyquist, these are normalized to give the correct amplitudes in Volts)
%        PSD        = matrix (freq x timestep, power has been normalized to conserve the energy in the signal)
%        freqs_FFT  = vector (vector of frequencies in Hz where the power was calculated)
%        all_phase  = array  (phase angles (in rads) of the form (frequency x window x taper))
%        t_sec      = vector (time (in seconds) at the center of each STFT window)
%
%Pending:
%       1.Update so that the phase is being calculated correctly, since this can be used for coherence methods when the 'multitapers' option is on.

%% parse variables
inp_pars = inputParser;
def_window = @rectwin;
expected_windows={'bartlett','barthannwin','blackman','blackmanharris','bohmanwin','chebwin',...
    'flattopwin','gausswin','hamming','hann','kaiser','nuttallwin','parzenwin','rectwin','taylorwin','tukeywin','triang'};
def_multitapers=0;
def_weights=0;

addParamValue(inp_pars,'window',def_window,@(x) (isa(x,'function_handle') && any(validatestring(func2str(x),expected_windows))));
addParamValue(inp_pars,'multitapers',def_multitapers,@(x) (isnumeric(x) && ~mod(x,0.5)));
addParamValue(inp_pars,'weights',def_weights,@(x) strcmp(x,'eigen'));

parse(inp_pars,varargin{:})
%inp_pars.Results
assert(~(inp_pars.Results.multitapers && ~any(strcmp('window',inp_pars.UsingDefaults))),'The ''window'' and ''multitaper'' options cannot be used together since for multitapers the DPSSs are used')
assert(~(any(strcmp('multitapers',inp_pars.UsingDefaults)) && ~any(strcmp('weights',inp_pars.UsingDefaults))),'The ''weights'' option is only valid when using the ''multitapers'' option')

%% data and parameters
X=X(:);
X(isnan(X))=[];
N = length(X);
Nsamples=round(T*Fs);
nfft=ceil((Nsamples+1)/2);
freqs_FFT=(0:nfft-1).*Fs/Nsamples;
%freqs_FFT=linspace(0,Fs/2,Nsamples/2+1);
assert(Nsamples<=N,'The lenght of the window (T) exceeds the duration of the clip (Fs*length(X))')
ovrlp_samples=round(overlap*Nsamples);
df=1/T;
dt=T/Nsamples;

%% set up tapers
if inp_pars.Results.multitapers
    taper_ind=1;
    V_weights=1;
    no_tapers = 2*inp_pars.Results.multitapers-1;
    if no_tapers<1 || no_tapers>Nsamples
        str=sprintf('\nthere are m DPSSs of length Nsamples, so require #tapers <= Nsamples per window and # tapers >= 1')
        error(str)
    end
    [E,V]=dpss(Nsamples,no_tapers/2);
    
    if strcmp(inp_pars.Results.weights,'eigen')
        clear V_weights
        V_weights=ones(1,ovrlp_samples>0,length(V));
        V_weights(1,1,:)=V/sum(V);
    end
else
    no_tapers=1;
    E=window(inp_pars.Results.window,Nsamples);
    taper_ind=0;
    V_weights=1;
end

% tic
% [Pxx F]=pmtm(X,3,[],Fs);
% toc

% pmtm is too slow

%% split signal into windows and apply tapers to each one
[eeg_windowed, Z] = buffer(X,Nsamples,ovrlp_samples,'nodelay'); %%NOTE: the second output argumen 'Z' is necessary to enforce full windows
% if ~isempty(Z)
%     length(Z)
% end
eeg_windowed = repmat(eeg_windowed,[1,1,no_tapers]);
windows=size(eeg_windowed,2);
E_windowed = permute(repmat(E,[1,1,windows]),[1 3 2]);
t_sec=nan(windows,2);
t_sec(:) = [((1:windows)-1)*(Nsamples-ovrlp_samples),((1:windows)*(Nsamples-ovrlp_samples)+ovrlp_samples)]/Fs;

%% get coefficients, power and phase
% freq x window x taper
fx = fft(E_windowed.*eeg_windowed);
fx = fx(1:nfft,:,:)/sqrt(Nsamples); %%% unitary scaling
all_pow = (fx.*conj(fx));%*(dt.^2);

all_pow = [all_pow(1,:,:);2*all_pow(2:nfft,:,:)];
if ~rem(Nsamples,2) % if Nsamples is odd, Nyquist component is not evaluated
    all_pow(end,:,:) = all_pow(end,:,:)/2;
end
all_phase=angle(fx);

% test_windowed=(sum(all_pow(:))-sum((E_windowed(:).*eeg_windowed(:)).^2))/(sum((E_windowed(:).*eeg_windowed(:).^2))); % Parseval check, for now only works if NFFT=Nsamples
% test_dim = sum((eeg(:).^2))-sum((eeg_windowed(:)).^2)/(Nsamples*windows/N); % not useful for tapers

%% mean phase (in rad), amplitude (in Volts) and power (in Volts^2)

% tempx=squeeze(sum(cos(all_phase),3)); %%%% this is definitely not correct for multitapers, haven't really checked for other cases.
% tempy=squeeze(sum(sin(all_phase),3));
% phase = atan2(tempy,tempx);

% coherent gain factors for amplitude and energy
taper_amp_norm=(1-(taper_ind && no_tapers==1))+(taper_ind && no_tapers==1)/mean(E(:));
taper_pow_norm=(1-(taper_ind && no_tapers==1))+(taper_ind && no_tapers==1)/rms(E(:));

% average over tapers, weighted by eigen values 'weights' is being used
PSD=mean(bsxfun(@times,all_pow,V_weights),3); % will average only if there is a third dim, i.e tapers.
% average over windows and scale amplitude
PS=mean(PSD,2);

PS = [sqrt(PS(1));sqrt(2*PS(2:end))]/sqrt(Nsamples)*(taper_amp_norm);
if ~rem(Nsamples,2)% if Nsamples is odd, Nyquist component is not evaluated
    PS(end)=PS(end)/sqrt(2);
end
% the division by (Nsamples*windows/N) in PSD is there to account for windowing, it is not used in PS since PS is an average over all windows, i.e amplitudes for windows aren't added. This is an approximation and may need to be fixed
PSD=(PSD*taper_pow_norm.^2);%/(Nsamples*windows/N);

% test_orig=(sum(PSD(:))-sum((eeg(:).^2)))/sum(eeg(:).^2);
if ~nargout
    figure;plot(freqs_FFT,PS)
    title('Spectrum')
    h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    set(gcf,'units','normalized','outerposition',[1 1 1 1])
    
    
    if ~isvector(PSD)
        figure;imagesc(mean(t_sec,2),freqs_FFT,10*log10(PSD))
        axis xy
        ylim([0 Fs/4])
        title('Spectrogram')
        h=gcf;
        set(h,'PaperOrientation','landscape');
        set(h,'PaperUnits','normalized');
        set(h,'PaperPosition', [0 0 1 1]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
    end
end
end