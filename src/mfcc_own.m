function cn = mfcc_own(s, fs, N, p, M, speaker, plot)
%{
 mfcc: Calculate the Mel-Frequency Ceptstral Coefficients
 In:
 s - Audio vector 
 fs - Sampling Frequency
 N - Hamming window size for stft
 p - # of filter banks
 M - stft overlap, in general M = N/3
 speaker - text, to label the graph in the mfcc plot. 
 plot - boolean, to decide if plotting the mfcc or not. 

 Out:
 cn - Mel-Frequency Ceptstral Coefficients
 y_t - Time-domain for STFT(Short-Time Fourier Transform)
%}

% A> STFT
% S = amplitude of stft, M x N
% F = frequency domain stft, M x 1
% T = time domain stft, N x 1
[S,F,T] = stft(s,fs,'Window',hamming(N),'OverlapLength',M,'FFTLength',N);


% B> Mel filter bank
% p: number of filters in filterbank
% n: length of fft
% fs: sample rate in Hz
m = melfb(p, N, fs);

% C> periodgram = abs(fft)^2
Sc = S((N/2):end, :);  % positive half of the frequency
melxstft = m * abs(Sc).^2; % matrix multiply mel with stft

% D> apply log and DCT
sk = log10(melxstft);           % mel spectrum coefficients
cn = dct(sk);                   % Discrete Cosine Transform

% E> Normalize mfcc
cn = cn ./ max(max(abs(cn))); % Using L-inf normalization

% F>  Plot
%check if the label for speaker exists, otherwise set it to a blank space
if ~exist('speaker', 'var') || isempty(speaker)
    speaker = ' ';
end

if plot
    figure;
    surf(T, 1:p, cn, 'EdgeColor','none'); 
    view(0, 90); 
    colorbar;
    caxis([-1 1]);
    xlim([min(T), max(T)]); 
    ylim([1 p]);
    xlabel('Time (s)'); 
    ylabel('mfc coefficients');
    title(speaker);

end