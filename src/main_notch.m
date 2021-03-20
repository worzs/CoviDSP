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

f_notch = 3000; % 0 Hz, 500 Hz, 1000 Hz, 3000 Hz
f_notch_bw = 500;


%% Testing
recognition_rate=zeros(numFiles, 1);
%% Read, load, normalize, mfcc for test signals
for i = 1:numFiles
    gen_sig = s{i};
    gen_sig = gen_sig + 0 * randn(size(gen_sig));
    gen_sig = normAudio(gen_sig);
    
    %notch fiter 
    df = 12.5*1000/(length(gen_sig));
    notch_n = floor(f_notch_bw/df);
    stp = floor(f_notch/df)+1;
    gen_sig_spec = fft(gen_sig);
    gen_sig_spec(stp:stp+notch_n) = 1e-20;
    gen_sig_spec((length(gen_sig) -(stp+notch_n):(length(gen_sig) - stp))) = 1e-20;
    gen_sig = ifft(gen_sig_spec);
    
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
