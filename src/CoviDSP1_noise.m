 %include: wavread, hamming, fft, dct, and own function melfb_own.m
clear;
clc;
close all;
N = 256; % window size
M = 128; % overlap
p = 16;  % number of filters in filterbank

lbg_p = 16; % length of the column vector for the lbg clustering. 
K = 4; % number of clusters
error_thresh = 0.001;
a = 0; % the variance of the noise

type_signal = 'edit'; %can be 'edit' or 'raw'. if not specified, 'edit' by default
signal_indexA = 5; %signal to plot their mfcc and the lbg clustering
signal_indexB = 8;
signal_indexC = 2;

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

SNR= [];
for i = 1:numFiles
    [sK,FssK]=audioread([directory, files{i}]);
    N_S = a*randn(size(sK)); % randn noise added 
    sK = sK + N_S;
    SNR = snr(sK, N_S); % signal/noise
    s{i} = sK;
    Fss{i} = FssK;
end
% for i = 1:numFiles
%     [s{i},Fss{i}]=audioread([directory, files{i}]);
% end

%% 2. eliminate quiet regions
% normalize and remove quiet regions at the beginning and in the end. 
s_n = cell(1,numFiles);
for i = 1:numFiles
    s_n{i}=normAudio(s{i});
end


% %% 3. Time domain plots
% %plot the time domain graphs for the selected speakers
% 
% %generate the time domain graph
% 
% 
% %speaker A
% figure(fig_count);
% fig_count = fig_count+1;
% subplot(1,2,1)
% t = (0:length(s{signal_indexA})-1)/Fss{signal_indexA};
% plot(t, s{signal_indexA}); 
% xlim([min(t), max(t)]);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title(['Time domain (original) ', files{signal_indexA}]);
% subplot(1,2,2)
% t = (0:length(s_n{signal_indexA})-1)/Fss{signal_indexA};
% plot(t, s_n{signal_indexA}); 
% xlim([min(t), max(t)]);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title(['Time domain (edited) ', files{signal_indexA}]);
% 
% %speaker B
% figure(fig_count);
% fig_count = fig_count+1;
% subplot(1,2,1)
% t = (0:length(s{signal_indexB})-1)/Fss{signal_indexB};
% plot(t, s{signal_indexB}); 
% xlim([min(t), max(t)]);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title(['Time domain (original) ', files{signal_indexB}]);
% subplot(1,2,2)
% t = (0:length(s_n{signal_indexB})-1)/Fss{signal_indexB};
% plot(t, s_n{signal_indexB}); 
% xlim([min(t), max(t)]);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title(['Time domain (edited) ', files{signal_indexB}]);
% 
% if plot_all
%     % plot in time domain, in sets of 4, of all the speakers.  
%     %compare the time domain graphs, before and after cropping the audio .
%     subplots = ceil(numFiles/4); % # of subplots
%     countFigs = numFiles; % # of figures to plot
%     for i = 1:subplots
%         figure(fig_count);
%         fig_count = fig_count+1;
% 
%         for j = 1:4
%             if countFigs >0
%                 countFigs = countFigs-1;
%                 index = 4*(i-1)+j; 
%                 subplot(2,4,j);
%                 plot(s{index});
%                 title([files{index} , ' original']);
% 
%                 subplot(2,4,j+4);
%                 plot(s_n{index});
%                 title([files{index} , ' edited']);
%             end
% 
% 
%         end 
%     end
% end


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


% %plot mfcc
% if plot_all
%     % plot the mfcc in sets of 8, of all the speakers. 
%     subplots = ceil(numFiles/8); % # of subplots
%     countFigs = numFiles; % # of figures to plot
%     for i = 1:subplots
%         figure(fig_count);
%         fig_count = fig_count+1;
% 
%         for j = 1:8
%             if countFigs >0
%                 countFigs = countFigs-1;
%                 index = 8*(i-1)+j; 
%                 subplot(2,4,j);
% 
%                 if strcmp(type_signal, 'raw')
%                     title([files{index}]);
%                     surf(T_raw{index}, 1:p, cn_raw_signal{index}, 'EdgeColor','none'); view(0, 45); colorbar;caxis([-1 1]);
%                     xlim([min(T_raw{index}), max(T_raw{index})]); ylim([1 p]);
%                 else 
%                     title([files{index} , ' edited']);
%                     surf(T_edit{index}, 1:p, cn_edit_signal{index}, 'EdgeColor','none'); view(0, 45); colorbar;caxis([-1 1]);
%                     xlim([min(T_edit{index}), max(T_edit{index})]); ylim([1 p]);
%                 end
%                 
%                 xlabel('Time(s)'); ylabel('mfcc');
% 
%                 %subplot(2,4,j+4);
%                 %plot(s_n{index});
%                 %title([files{index} , ' edited']);
%             end
% 
% 
%         end 
%     end
% end
% 
% 
% %Plot the spectrograms and the mfcc for the selected speakers
% figure(fig_count);
% fig_count = fig_count+1;
% subplot(2,2,1)
% [Y, F, T, P] = spectrogram (s{signal_indexA}(:,1), hamming(N,'periodic'), M, N, Fss{signal_indexA});
% surf(T_raw{signal_indexA},F,20*log10(abs(P)),'EdgeColor','none');
% %axis tight;
% view(angle1, angle2); colorbar; %caxis([-60 0]);
% xlabel('Time[s]'); ylabel('Frequency [Hz]');zlabel('Amplitude [dB]')
% title(['Spectrogram ', files{signal_indexA}, ' original']);
% 
% subplot(2,2,2)
% %title([files{signal_indexA}]);
% surf(T_raw{signal_indexA}, 1:p, 20*log10(abs(cn_raw_signal{signal_indexA})), 'EdgeColor','none'); view(angle1, angle2); colorbar;caxis([-60 0]);
% xlim([min(T_raw{signal_indexA}), max(T_raw{signal_indexA})]); ylim([1 p]);
% xlabel('Time[s]'); ylabel('MFCC');zlabel('Amplitude [dB]')
% title(['MFCC ', files{signal_indexA}, ' original']);
% ylim([2,p])
% 
% subplot(2,2,3);
% %spectrogram(s_n{signal_indexA}(:,1));
% [Y, F, T, P] = spectrogram (s_n{signal_indexA}(:,1), hamming(N,'periodic'), M, N, Fss{signal_indexA});
% surf(T_edit{signal_indexA},F,20*log10(abs(P)),'EdgeColor','none');
% %axis tight;
% view(angle1, angle2); colorbar; %caxis([-60 0]);
% xlabel('Time[s]'); ylabel('Frequency [Hz]');zlabel('Amplitude [dB]')
% title(['Spectrogram ', files{signal_indexA}, ' edited']);
% 
% subplot(2,2,4)
% %title([files{signal_indexA}]);
% surf(T_edit{signal_indexA}, 1:p, 20*log10(abs(cn_edit_signal{signal_indexA})), 'EdgeColor','none'); view(angle1, angle2); colorbar;caxis([-60 0]);
% xlim([min(T_edit{signal_indexA}), max(T_edit{signal_indexA})]); ylim([1 p]);
% xlabel('Time[s]'); ylabel('MFCC');zlabel('Amplitude [dB]')
% title(['MFCC ', files{signal_indexA}, ' edited']);
% ylim([2,p])
% 
% % plot the filter bank
% figure(fig_count);
% fig_count = fig_count+1;
% plot(linspace(0,(Fss{signal_indexA}/2), N/2+1), melfb(p, N, Fss{signal_indexA})');
% title('Mel filterbank');
% xlabel('Frequency (Hz)');

%% 5. Clustering via LBG algorithm & k-means

% pick the MFCC for speakers the selected speakers.
% select the vectors to plot, based on the original signal 'raw' or the
% cropped 'edit'. 
if strcmp(type_signal, 'raw')
    cn_signal = cn_raw_signal;
else 
    cn_signal = cn_edit_signal;
end

% lbg algorithm 
% new_centroids = lbg(samples, M_max, step_size, error_threshold)

centroids_codebook = zeros(numFiles, K, lbg_p);
for i = 1: numFiles
    S_N = cn_signal{i}(1:lbg_p, :)';
    centroids_N = lbg(S_N, K, 0.01, error_thresh);
    centroids_codebook(i, :, :) = centroids_N;
end

%% 6. plot 2D with any two speakers, the MFCC and the centroids. 
% S_A = cn_signal{signal_indexA}(1:lbg_p,:)';
% S_B = cn_signal{signal_indexB}(1:lbg_p,:)';
% S_C = cn_signal{signal_indexC}(1:lbg_p,:)';
% 
% % lbg algorithm 
% % new_centroids = lbg(samples, M_max, step_size, error_threshold)
% centroids_A = lbg(S_A, K, 0.01, error_thresh);
% centroids_B = lbg(S_B, K, 0.01, error_thresh);
% 
% %plot MFCC and centroids for the first speaker
% figure(fig_count);
% fig_count = fig_count+1;
% plot(cn_signal{signal_indexA}(dim1_signal,:)', cn_signal{signal_indexA}(dim2_signal,:)','ro');
% hold on;
% plot(centroids_A(:,dim1_signal)', centroids_A(:,dim2_signal)','k*');
% xlabel(['mfcc-',num2str(dim1_signal)]); ylabel(['mfcc-',num2str(dim2_signal)]);
% legend(files{signal_indexA}, 'centroids');
% grid on;
% title(["MFCC ", files{signal_indexA}]);
% xlim([-1 1]);
% ylim([-1 1]);
% hold off;

% %plot MFCC and centroids for the second speaker
% figure(fig_count);
% fig_count = fig_count+1;
% plot(cn_signal{signal_indexB}(dim1_signal,:)', cn_signal{signal_indexB}(dim2_signal,:)','ro');
% hold on;
% plot(centroids_B(:,dim1_signal)', centroids_B(:,dim2_signal)','b*');
% xlabel(['mfcc-',num2str(dim1_signal)]); ylabel(['mfcc-',num2str(dim2_signal)]);
% legend(files{signal_indexB}, 'centroids');
% grid on;
% title(["MFCC ", files{signal_indexB}]);
% xlim([-1 1]);
% ylim([-1 1]);
% hold off;



for i = 7
    centroids_N=[];
    S_N = cn_signal{i}(1:lbg_p,:)';
    centroids_N(:, :) = lbg(S_N, K, 0.01, error_thresh);
    figure(fig_count);
    fig_count = fig_count+1;
    plot(cn_signal{i}(dim1_signal,:)', cn_signal{i}(dim2_signal,:)','ro');
    hold on;
    plot(centroids_N(:,dim1_signal)', centroids_N(:,dim2_signal)','b*');
    xlabel(['mfcc-',num2str(dim1_signal)]); ylabel(['mfcc-',num2str(dim2_signal)]);
    legend(files{i}, 'centroids');
    grid on;
    title(["TEST MFCC ", files{i}]);
    xlim([-1 1]);
    ylim([-1 1]);
    hold off;
end

%% Testing
% define counters
numFiles = 8; %number of test files
% define directory - subfolder
directory = './Test/';
% 1. Read test signals
% First, create the array with the names of the files. Assume that all the
% files follow the standard 's<i>.wav', where <i> is the identifier of the
% speaker. 
%% TEST Buffer cells
t_files = cell(1,numFiles);
t_s = cell(1,numFiles);
t_Fss = cell(1,numFiles);
t_s_n = cell(1,numFiles);
t_cn_raw_signal = cell(1,numFiles);
t_cn_edit_signal = cell(1,numFiles);
t_T_raw = cell(1,numFiles); 
t_T_edit = cell(1,numFiles);


recognition_rate=zeros(1, numFiles);
%% Read, load, normalize, mfcc for test signals
for i = 1:numFiles
    t_files{i} = ['s',num2str(i),'.wav'];
    [t_s{i}, t_Fss{i}] = audioread([directory, t_files{i}]);  
    t_s_n{i} = normAudio(t_s{i}); 
    [t_cn_raw_signal{i}, t_T_raw{i}] = mfcc_own(t_s{i}(:,1) - mean(t_s{i}(:,1)), t_Fss{i}, N, p, M);
    [t_cn_edit_signal{i}, t_T_edit{i}] = mfcc_own(t_s_n{i}(:,1), t_Fss{i}, N, p, M);

    if strcmp(type_signal, 'raw')
        t_cn_signal = t_cn_raw_signal{i};
    else 
        t_cn_signal = t_cn_edit_signal{i};
    end
    

    %% plot test signals
    t_centroids_N=[];
    t_S_N = t_cn_signal(1:lbg_p,:)';
    t_centroids_N(:, :) = lbg(t_S_N, K, 0.01, error_thresh);   %t_new_centroids = lbg(samples, M_max, step_size, error_threshold);
    figure(fig_count);
    fig_count = fig_count+1;
    plot(t_cn_signal(dim1_signal,:)', t_cn_signal(dim2_signal,:)','ro');
    hold on;
    plot(t_centroids_N(:,dim1_signal)', t_centroids_N(:,dim2_signal)','b*');
    xlabel(['mfcc-',num2str(dim1_signal)]); ylabel(['mfcc-',num2str(dim2_signal)]);
    legend(t_files{i}, 'centroids');
    grid on;
    title(["TEST MFCC ", t_files{i}]);
    xlim([-1 1]);
    ylim([-1 1]);
    hold off;
     
    test_edit_signal = t_cn_signal(1:lbg_p,:)';
    err_vec = zeros(1, numFiles);
    for j = 1: numFiles
        codebook = squeeze(centroids_codebook(j, :, :));
        % centroids_codebook is the array of training-signals centroids line 426
        % codebook is the matrix of each training signal
        t_num = length(test_edit_signal(:, 1));                                % time of the test signal
        for t = 1: t_num
            dist = zeros(1, K);
            for k = 1: K
                dist(1, k) = norm(test_edit_signal(t, :) - codebook(k, :), 2); % Euclidean distance between test signals and codebooks
            end
            [val, ~] = min(dist);                                              % minimum mean distance with each centroids
            err_vec(j) = err_vec(j) + val;                                     % accumulate error for the whole time
        end
    end
    [val, ind] = min(err_vec);                                                 % find the closest codebook index
    recognition_rate(i) = ind; 
end
disp(recognition_rate);




