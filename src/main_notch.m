%include: wavread, hamming, fft, dct, and own function melfb_own.m
clc;
close all;
N = 128; % window size
M = 100; % overlap
p = 16;  % number of filters in filterbank

lbg_p = 15; % length of the column vector for the lbg clustering. 
K = 16; % number of clusters
error_thresh = 0.05;

type_signal = 'edit'; %can be 'edit' or 'raw'. if not specified, 'edit' by default

%notch filter parameters
f_notch = 500; % frequency in hertz 
f_notch_bw = 500;

dim1_signal = 1;    %dimensions to plot in the mfcc
dim2_signal = 3;

plot_all = false; % boolean to plot the graphs of all the speakers

angle1 = 135; %angles to rotate the 3d plots
angle2 = 60;
%% define counters
numFiles = 11; %number of files
fig_count = 1; % initialize figure counter

%% define directory - subfolder
directory = './Train/';

%% 1. Read signals
% First, create the array with the names of the files. Assume that all the
% files follow the standard 's<i>.wav', where <i> is the identifier of the
% speaker. 
files = cell(1,numFiles);
for i = 1:numFiles
    files{i} = ['s',num2str(i),'.wav'];
end

% Then, load the information of the signal and the sampling rate
s = cell(1,numFiles);
Fss = cell(1,numFiles);


% save raw data for generating new testing sample using filters

for i = 1:numFiles
    [s{i},Fss{i}]=audioread([directory, files{i}]);

end

%% 2. eliminate quiet regions
% normalize and remove quiet regions at the beginning and in the end. 
s_n = cell(1,numFiles);
for i = 1:numFiles
    s_n{i}=normAudio(s{i});
end




%% 4. obtain mel coefficients

%initialize the cells to store the data
cn_raw_signal = cell(1,numFiles);
cn_edit_signal = cell(1,numFiles);
T_raw = cell(1,numFiles); 
T_edit = cell(1,numFiles); 

%obtain the mfcc 
for i = 1:numFiles
    [cn_raw_signal{i},T_raw{i}]=mfcc_own(s{i}(:,1) - mean(s{i}(:,1)), Fss{i}, N, p, M);
    [cn_edit_signal{i},T_edit{i}]=mfcc_own(s_n{i}(:,1), Fss{i}, N, p, M);
end


%% 5. Clustering via LBG algorithm & k-means

% pick the MFCC for speakers the selected speakers.
% select the vectors to plot, based on the original signal 'raw' or the
% cropped 'edit'. 
if strcmp(type_signal, 'raw')
    cn_signal = cn_raw_signal;
else 
    cn_signal = cn_edit_signal;
end

centroids_codebook = zeros(11, K, 15);
for i = 1:11
    S_N = cn_signal{i}(1:lbg_p,:)';
    centroids_N = lbg(S_N, K, 0.01, 0.001);
    centroids_codebook(i, :, :) = centroids_N;
end




%% Testing
recognition_rate=zeros(numFiles, 1);
%% Read, load, normalize, mfcc for test signals
for i = 1:numFiles
    %get the signal and apply the normalization step. 
    gen_sig = s{i};
    %gen_sig = gen_sig + 0 * randn(size(gen_sig));
    gen_sig = normAudio(gen_sig);
    
    %notch fiter 
    %https://github.com/heguanda2/EEC201_genliu/blob/master/main_script_notch.m
    df = Fss{i}/(length(gen_sig)); %get the delta in frequency
    notch_n = floor(f_notch_bw/df); %get the index of the notch frequency
    step = floor(f_notch/df)+1; %get the step
    gen_sig_spec = fft(gen_sig); %apply fft to the signal
    gen_sig_spec(step:step+notch_n) = db2mag(-100);
    gen_sig_spec((length(gen_sig) -(step+notch_n):(length(gen_sig) - step))) = db2mag(-100); %flat the frequency response around f_notch
    gen_sig = ifft(gen_sig_spec);
    
    
    figure()
    %plot the frequency response
    plot(-1:2/(length(gen_sig)-1):1 , 20*log10(abs(fftshift(gen_sig_spec))),'k');
    xlabel('Normalized Frequency \omega/\pi'); ylabel('Magnitude [dB]');
    grid on;
    title(['Frequency response ', files{i}]);
    xlim([0, 1]);
    %ylim([0 1]);
    hold off;
    
    figure()
    [Y, F, T, P] = spectrogram (gen_sig, hamming(N,'periodic'), M, N, Fss{i});
    surf(T_edit{i},F,20*log10(abs(P)),'EdgeColor','none');
    %axis tight;
    view(angle1, angle2); colorbar; %caxis([-60 0]);
    xlabel('Time[s]'); ylabel('Frequency [Hz]');zlabel('Amplitude [dB]')
    title(['Spectrogram ', files{i}, ' edited']);
    ylim([0,Fss{i}/2]);


    if strcmp(type_signal, 'raw')
        [gen_mfcc, ~] = mfcc_own(gen_sig(:,1) - mean(gen_sig(:,1)), Fss{i}, N, p, M);
    else 
        [gen_mfcc, ~] = mfcc_own(gen_sig(:, 1), Fss{i}, N, p, M);
    end
    
    %% TEST classification performance using trained codebook
    test_edit_signal = gen_mfcc(1:lbg_p,:)';                            % turn cell into matrix  
    err_vec = zeros(1, numFiles);
    for j = 1: numFiles
        codebook = squeeze(centroids_codebook(j,:,:));
        % centroids_codebook is the array of training-signals centroids line 426
        % codebook is the matrix of each training signal
        t_num = length(test_edit_signal(:,1));                                 % time of the test signal
        for c = 1: t_num
            dist = zeros(1, K);
            for k = 1: K
                dist(1, k) = norm(test_edit_signal(c, :) - codebook(k, :), 2); % Euclidean distance between test signals and codebooks
            end
            [val, ~] = min(dist);                                            % minimum mean distance with each centroids
            err_vec(j) = err_vec(j) + val;                                     % accumulate error for the whole time
        end
    end
    [val, ind] = min(err_vec);                                                 % find the closest codebook index
        
    recognition_rate(i) = ind == i; 
end
display(recognition_rate');
