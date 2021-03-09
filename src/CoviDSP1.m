 %include: wavread, hamming, fft, dct, and own function melfb_own.m
clear;
clc;

addpath(genpath(pwd));

%% 1. Read signal
[s1, Fss1] = audioread('s1.wav');
[s2, Fss2] = audioread('s2.wav');
[s3, Fss3] = audioread('s3.wav');
[s4, Fss4] = audioread('s4.wav');
[s5, Fss5] = audioread('s5.wav');
[s6, Fss6] = audioread('s6.wav');
[s7, Fss7] = audioread('s7.wav');
[s8, Fss8] = audioread('s8.wav');
[s9, Fss9] = audioread('s9.wav');
[s10, Fss10] = audioread('s10.wav');
[s11, Fss11] = audioread('s11.wav');

%% 2. eliminate quiet regions
s1=s1(2201:10200, 1);
s2=s2(2201:10200, 1);
s3=s3(2201:10200, 1);
s4=s4(2501:11500, 1);
s5=s5(4801:14800, 1);
s6=s6(3001:12000, 1);
s7=s7(2501:11500, 1);
s8=s8(2601:11600, 1);
s9=s9(4701:14700, 1);
s10=s10(6501:14500, 1);
s11=s11(8001:16000, 1);

% Test2: plot the signal to view it in the time domain
%{

 figure(1)
 plot(s1)
 figure(2)
 plot(s2)
 figure(3)
 plot(s3)
 figure(4)
 plot(s4)
 figure(5)
 plot(s5)
 figure(6)
 plot(s6)
 figure(7)
 plot(s7)
 figure(8)
 plot(s8)
 figure(9)
 plot(s9)
 figure(10)
 plot(s10)
 figure(11)
 plot(s11)
%}


%% 3. Time domain
%{
t = (0:length(s1)-1)/Fss1;
figure; 
plot(t, s1); 
xlim([min(t), max(t)]);
xlabel('Time (s)'); 
ylabel('Amplitude'); 
title('Time domain');
%}

%% 4. obtain mel coefficients
N = 512; % window size
M = 200; % overlap
p = 20;  % number of filters in filterbank
%replace previous code with function
%{
%% 1 STFT
%s2 = s2./ max(s2);
%N = 256; %window size
%M = 100; %overlap
[S,F,T] = stft(s2,Fss2,'Window',hamming(N),'OverlapLength',M,'FFTLength',N);
% S = amplitude output of stft
% F = frequency
% T = time

% 3D plot
%{
 %surf(F,T,abs(S(:,:,1))')
 %%mesh(F,T,abs(S(:,:,1))');
 shading interp;
 xlabel('Frequency (Hz)');
 ylabel('Time (seconds)');
 zlabel('Amplitude');
%}

%% 2 mel-frequency filter bank
p = 20;                         % number of filters in filterbank
n = 256;                        % n length of fft
fs = Fss1;                      % fs sample rate in Hz
m = melfb(p, n, fs);

%% 3 periodgram = abs(fft)^2
Sc = S((n/2):end, :);           % positive half of the frequency
melxstft = m * abs(Sc).^2;      % matrix multiply mel with stft

%% 4
sk = log10(melxstft);           % mel spectrum coefficients
cn = dct(sk);                   % Discrete Cosine Transform

%% 5
figure; 
surf(T, 1:p, cn./ max(max(abs(cn))),'EdgeColor','none'); 
view(0, 90); 
colorbar;
caxis([-1 1]);
xlim([min(T), max(T)]); 
ylim([1 p]);
xlabel('Time (s)'); 
ylabel('mfc coefficients');
title('s1');

%}

cn1 = mfcc_own(s1, Fss1, N, p, M, 's1', true);
cn2 = mfcc_own(s2, Fss2, N, p, M, 's2', true);
cn3 = mfcc_own(s3, Fss3, N, p, M, 's3', true);
cn4 = mfcc_own(s4, Fss4, N, p, M, 's4', true);
cn10 = mfcc_own(s10, Fss10, N, p, M, 's10', true);

%% 5. plot 2D with any two speakers ,2D

%TODO: how to choose the index of the mel fb to plot.
%Also, why only some indices seem to be in cluster, while others not. 

figure;
%plot(cn1(2,:)', cn1(3,:)', 'x');
hold on;
%plot(cn2(2,:)', cn2(3,:)', '*');
plot(cn3(2,:)', cn3(3,:)', '+');
%plot(cn4(2,:)', cn4(3,:)', '.');
plot(cn10(2,:)', cn10(3,:)', 'o');
xlabel('mfcc-2'); ylabel('mfcc-3');
%legend("Speaker 1", "Speaker 2", "Speaker 3", "Speaker 4", "Speaker 10");
legend("Speaker 3", "Speaker 10");
grid on;
title("mfcc space");
hold off


