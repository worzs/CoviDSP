 %include: wavread, hamming, fft, dct, and own function melfb_own.m
clc;
close all;
N = 256; % window size
M = 100; % overlap
p = 20;  % number of filters in filterbank

lbg_p = 15; % length of the column vector for the lbg clustering. 
K = 32; % number of clusters
error_thresh = 0.05;
start_index_p = 2; %index to remove the first MFCC

type_signal = 'edit'; %can be 'edit' or 'raw'. if not specified, 'edit' by default
signal_indexA = 3; %signal to plot their mfcc and the lbg clustering
signal_indexB = 10;
signal_indexC = 2;

dim1_signal = 3;    %dimensions to plot in the mfcc
dim2_signal = 5;


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
for i = 1:numFiles
    [s{i},Fss{i}]=audioread([directory, files{i}]);
end
 
%% 2. eliminate quiet regions
% normalize and remove quiet regions at the beginning and in the end. 
s_n = cell(1,numFiles);
for i = 1:numFiles
    s_n{i}=normAudio(s{i});
end

%{
s1_n = normAudio(s1);
s2_n = normAudio(s2);
s3_n = normAudio(s3);
s4_n = normAudio(s4);
s5_n = normAudio(s5);
s6_n = normAudio(s6);
s7_n = normAudio(s7);
s8_n = normAudio(s8);
s9_n = normAudio(s9);
s10_n = normAudio(s10);
s11_n = normAudio(s11);


% then remove quiet regions
% s1_n = quitSilence(s1_n);


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


%% 3. Time domain plots
%{
t = (0:length(s1)-1)/Fss1;
figure; 
plot(t, s1); 
xlim([min(t), max(t)]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Time domain');
%}

%plot the time domain graphs for the selected speakers

%generate the time domain graph


%speaker A
figure(fig_count);
fig_count = fig_count+1;
subplot(1,2,1)
t = (0:length(s{signal_indexA})-1)/Fss{signal_indexA};
plot(t, s{signal_indexA}); 
xlim([min(t), max(t)]);
xlabel('Time (s)');
ylabel('Amplitude');
title(['Time domain (original) ', files{signal_indexA}]);
subplot(1,2,2)
t = (0:length(s_n{signal_indexA})-1)/Fss{signal_indexA};
plot(t, s_n{signal_indexA}); 
xlim([min(t), max(t)]);
xlabel('Time (s)');
ylabel('Amplitude');
title(['Time domain (edited) ', files{signal_indexA}]);

%speaker B
figure(fig_count);
fig_count = fig_count+1;
subplot(1,2,1)
t = (0:length(s{signal_indexB})-1)/Fss{signal_indexB};
plot(t, s{signal_indexB}); 
xlim([min(t), max(t)]);
xlabel('Time (s)');
ylabel('Amplitude');
title(['Time domain (original) ', files{signal_indexB}]);
subplot(1,2,2)
t = (0:length(s_n{signal_indexB})-1)/Fss{signal_indexB};
plot(t, s_n{signal_indexB}); 
xlim([min(t), max(t)]);
xlabel('Time (s)');
ylabel('Amplitude');
title(['Time domain (edited) ', files{signal_indexB}]);

if plot_all
    % plot in time domain, in sets of 4, of all the speakers.  
    %compare the time domain graphs, before and after cropping the audio .
    subplots = ceil(numFiles/4); % # of subplots
    countFigs = numFiles; % # of figures to plot
    for i = 1:subplots
        figure(fig_count);
        fig_count = fig_count+1;

        for j = 1:4
            if countFigs >0
                countFigs = countFigs-1;
                index = 4*(i-1)+j; 
                subplot(2,4,j);
                plot(s{index});
                title([files{index} , ' original']);

                subplot(2,4,j+4);
                plot(s_n{index});
                title([files{index} , ' edited']);
            end


        end 
    end
end

%{
figure(fig_count);
fig_count = fig_count+1;

subplot(2,5,1);

subplot(2,5,2);
plot(s2);
title('s2');
subplot(2,5,3);
plot(s3);
title('s3');
subplot(2,5,4);
plot(s4);
title('s4');
subplot(2,5,5);
plot(s5);
title('s5');

subplot(2,5,6);
plot(s1_n);
title('s1_n');
subplot(2,5,7);
plot(s2_n);
title('s2_n');
subplot(2,5,8);
plot(s3_n);
title('s3_n');
subplot(2,5,9);
plot(s4_n);
title('s4_n');
subplot(2,5,10);
plot(s5_n);
title('s5_n');

figure(fig_count);
fig_count = fig_count+1;

subplot(2,5,1);
plot(s6);
title('s6');
subplot(2,5,2);
plot(s7);
title('s7');
subplot(2,5,3);
plot(s8);
title('s8');
subplot(2,5,4);
plot(s9);
title('s9');
subplot(2,5,5);
plot(s11);
title('s11');

subplot(2,5,6);
plot(s6_n);
title('s6_n');
subplot(2,5,7);
plot(s7_n);
title('s7_n');
subplot(2,5,8);
plot(s8_n);
title('s8_n');
subplot(2,5,9);
plot(s9_n);
title('s9_n');
subplot(2,5,10);
plot(s11_n);
title('s11_n');
%}


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
 
%{
[cn1, T1] = mfcc_own(s1, Fss1, N, p, M);
[cn1_n,T1_n] = mfcc_own(s1_n, Fss1, N, p, M);
[cn2_n,T2_n] = mfcc_own(s2_n, Fss2, N, p, M);
[cn3_n,T3_n] = mfcc_own(s3_n, Fss3, N, p, M);
[cn4_n,T4_n] = mfcc_own(s4_n, Fss4, N, p, M);
[cn5_n,T5_n] = mfcc_own(s5_n, Fss5, N, p, M);
[cn6_n,T6_n] = mfcc_own(s6_n, Fss6, N, p, M);
[cn7_n,T7_n] = mfcc_own(s7_n, Fss7, N, p, M);
[cn8_n,T8_n] = mfcc_own(s8_n, Fss8, N, p, M);
[cn9_n,T9_n] = mfcc_own(s9_n, Fss9, N, p, M);
[cn10_n,T10_n] = mfcc_own(s10_n, Fss10, N, p, M);
[cn11_n,T11_n] = mfcc_own(s11_n, Fss11, N, p, M);
%}

%plot mfcc
if plot_all
    % plot the mfcc in sets of 8, of all the speakers. 
    subplots = ceil(numFiles/8); % # of subplots
    countFigs = numFiles; % # of figures to plot
    for i = 1:subplots
        figure(fig_count);
        fig_count = fig_count+1;

        for j = 1:8
            if countFigs >0
                countFigs = countFigs-1;
                index = 8*(i-1)+j; 
                subplot(2,4,j);

                
                

                if strcmp(type_signal, 'raw')
                    title([files{index}]);
                    surf(T_raw{index}, 1:p, cn_raw_signal{index}, 'EdgeColor','none'); view(0, 45); colorbar;caxis([-1 1]);
                    xlim([min(T_raw{index}), max(T_raw{index})]); ylim([1 p]);
                else 
                    title([files{index} , ' edited']);
                    surf(T_edit{index}, 1:p, cn_edit_signal{index}, 'EdgeColor','none'); view(0, 45); colorbar;caxis([-1 1]);
                    xlim([min(T_edit{index}), max(T_edit{index})]); ylim([1 p]);
                end
                
                xlabel('Time(s)'); ylabel('mfcc');

                %subplot(2,4,j+4);
                %plot(s_n{index});
                %title([files{index} , ' edited']);
            end


        end 
    end
end




%Plot the spectrograms and the mfcc for the selected speakers
figure(fig_count);
fig_count = fig_count+1;
subplot(2,2,1)
[Y, F, T, P] = spectrogram (s{signal_indexA}(:,1), hamming(N,'periodic'), M, N, Fss{signal_indexA});
surf(T_raw{signal_indexA},F,20*log10(abs(P)),'EdgeColor','none');
%axis tight;
view(angle1, angle2); colorbar; %caxis([-60 0]);
xlabel('Time[s]'); ylabel('Frequency [Hz]');zlabel('Amplitude [dB]')
title(['Spectrogram ', files{signal_indexA}, ' original']);

subplot(2,2,2)
%title([files{signal_indexA}]);
surf(T_raw{signal_indexA}, 1:p, 20*log10(abs(cn_raw_signal{signal_indexA})), 'EdgeColor','none'); view(angle1, angle2); colorbar;caxis([-60 0]);
xlim([min(T_raw{signal_indexA}), max(T_raw{signal_indexA})]); ylim([1 p]);
xlabel('Time[s]'); ylabel('MFCC');zlabel('Amplitude [dB]')
title(['MFCC ', files{signal_indexA}, ' original']);
ylim([start_index_p,p])

subplot(2,2,3);
%spectrogram(s_n{signal_indexA}(:,1));
[Y, F, T, P] = spectrogram (s_n{signal_indexA}(:,1), hamming(N,'periodic'), M, N, Fss{signal_indexA});
surf(T_edit{signal_indexA},F,20*log10(abs(P)),'EdgeColor','none');
%axis tight;
view(angle1, angle2); colorbar; %caxis([-60 0]);
xlabel('Time[s]'); ylabel('Frequency [Hz]');zlabel('Amplitude [dB]')
title(['Spectrogram ', files{signal_indexA}, ' edited']);

subplot(2,2,4)
%title([files{signal_indexA}]);
surf(T_edit{signal_indexA}, 1:p, 20*log10(abs(cn_edit_signal{signal_indexA})), 'EdgeColor','none'); view(angle1, angle2); colorbar;caxis([-60 0]);
xlim([min(T_edit{signal_indexA}), max(T_edit{signal_indexA})]); ylim([1 p]);
xlabel('Time[s]'); ylabel('MFCC');zlabel('Amplitude [dB]')
title(['MFCC ', files{signal_indexA}, ' edited']);
ylim([start_index_p,p])

% plot the filter bank
figure(fig_count);
fig_count = fig_count+1;
plot(linspace(0,(Fss{signal_indexA}/2), N/2+1), melfb(p, N, Fss{signal_indexA})');
title('Mel filterbank');
xlabel('Frequency (Hz)');

%{
%for i = 1:numFiles
%    subplot(2,5,i);
%    surf(T1_n, 1:p, cn1_n, 'EdgeColor','none'); view(0, 90); colorbar;caxis([-1 1]);
%    xlim([min(T1_n), max(T1_n)]); ylim([1 p]);
%    xlabel('Time (s)'); ylabel('mfcc');
%end

%subplot(2,5,2);
%surf(T2_n, 1:p, cn1_n, 'EdgeColor','none'); view(0, 90); colorbar;caxis([-1 1]);
%xlim([min(T2_n), max(T1_n)]); ylim([1 p]);
%xlabel('Time (s)'); ylabel('mfcc');
%}



%% 5. Clustering via LBG algorithm & k-means

%pick the MFCC for speakers the selected speakers.
% select the vectors to plot, based on the original signal 'raw' or the
% cropped 'edit'. 
if strcmp(type_signal, 'raw')
    cn_signal = cn_raw_signal;
else 
    cn_signal = cn_edit_signal;
end

S_A = cn_signal{signal_indexA}(start_index_p:lbg_p,:)';
S_B = cn_signal{signal_indexB}(start_index_p:lbg_p,:)';
S_C = cn_signal{signal_indexC}(start_index_p:lbg_p,:)';

% lbg algorithm 
% new_centroids = lbg(samples, M_max, step_size, error_threshold)
centroids_A = lbg(S_A, K, 0.01, 0.001);
centroids_B = lbg(S_B, K, 0.01, 0.001);

centroids_codebook = zeros(11, K, 15);
for i = 1:8
    S_N = cn_signal{i}(1:lbg_p,:)';
    centroids_N = lbg(S_N, K, 0.01, 0.001);
    centroids_codebook(i, :, :) = centroids_N;
end
%{  
for i = 1:11
    centroids_N=[];
    S_N = cn_signal{i}(1:lbg_p,:)';
    centroids_N(:, :) = lbg(S_N, K, 0.01, 0.001);
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
%}
%% 6. plot 2D with any two speakers, the MFCC and the centroids. 

%plot MFCC and centroids for the first speaker
figure(fig_count);
fig_count = fig_count+1;
plot(cn_signal{signal_indexA}(dim1_signal,:)', cn_signal{signal_indexA}(dim2_signal,:)','ro');
hold on;
plot(centroids_A(:,dim1_signal-1)', centroids_A(:,dim2_signal-1)','k*');
xlabel(['mfcc-',num2str(dim1_signal)]); ylabel(['mfcc-',num2str(dim2_signal)]);
legend(files{signal_indexA}, 'centroids');
grid on;
title(["MFCC ", files{signal_indexA}]);
xlim([-1 1]);
ylim([-1 1]);
hold off;

%plot MFCC and centroids for the second speaker
figure(fig_count);
fig_count = fig_count+1;
plot(cn_signal{signal_indexB}(dim1_signal,:)', cn_signal{signal_indexB}(dim2_signal,:)','ro');
hold on;
plot(centroids_B(:,dim1_signal-1)', centroids_B(:,dim2_signal-1)','b*');
xlabel(['mfcc-',num2str(dim1_signal)]); ylabel(['mfcc-',num2str(dim2_signal)]);
legend(files{signal_indexB}, 'centroids');
grid on;
title(["MFCC ", files{signal_indexB}]);
xlim([-1 1]);
ylim([-1 1]);
hold off;


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

recognition_rate=zeros(11, 1);
test_num_max = 1;
test_num = 0;
a=0;
%% Read, load, normalize, mfcc for test signals
for i = 1:numFiles
    t_files{i} = ['s',num2str(i),'.wav'];
    [t_s{i},t_Fss{i}]=audioread([directory, t_files{i}]);
    for j_1=1:test_num_max    
        test_num = test_num + 1; 
        err_vec = zeros(11,1); %matrix use for error in testing
        t_s_n{i}=normAudio(t_s{i}); 
        [t_cn_raw_signal{i},t_T_raw{i}]=mfcc_own(t_s{i}(:,1) - mean(t_s{i}(:,1)), t_Fss{i}, N, p, M);
        [t_cn_edit_signal{i}, t_T_edit{i}]=mfcc_own(t_s_n{i}(:,1), t_Fss{i}, N, p, M);

        if strcmp(type_signal, 'raw')
            t_cn_signal = t_cn_raw_signal;
        else 
            t_cn_signal = t_cn_edit_signal;
        end
        %t_new_centroids = lbg(samples, M_max, step_size, error_threshold);
        %% plot test signals
        t_centroids_N=[];
        t_S_N = t_cn_signal{i}(1:lbg_p,:)';
        t_centroids_N(:, :) = lbg(t_S_N, K, 0.01, 0.001);
        figure(fig_count);
        fig_count = fig_count+1;
        plot(t_cn_signal{i}(dim1_signal,:)', t_cn_signal{i}(dim2_signal,:)','ro');
        hold on;
        plot(t_centroids_N(:,dim1_signal)', t_centroids_N(:,dim2_signal)','b*');
        xlabel(['mfcc-',num2str(dim1_signal)]); ylabel(['mfcc-',num2str(dim2_signal)]);
        legend(t_files{i}, 'centroids');
        grid on;
        title(["TEST MFCC ", t_files{i}]);
        xlim([-1 1]);
        ylim([-1 1]);
        hold off;
        
        %%
        test_edit_signal = cn_signal{i}(1:lbg_p,:)';     % turn cell into matrix   
        match_num = 0;
        for k=1:11
            codebook=squeeze(centroids_codebook(k,:,:)); 
            % centroids_codebook is the array of training-signals centroids line 414
            % codebook is the matrix of each training signal
                                                        
            err = 0;
            t_num = length(test_edit_signal(:,1));  % time of the test signal
            codebook_num = length(codebook(:, 1));  % time of the codebook
            for c=1:t_num
                cur_err = [];
                for j=1:codebook_num
                    cur_err = [cur_err; norm(test_edit_signal(c,:) - codebook(j,:),2)]; % Euclidean distance between test signals and codebooks
                end
                [val, ind] = min(cur_err); % minimum error with each time slot
                err = err + val;           % error for the whole time
            end
            err = err/t_num;               % mean error with each of the codebook
            err_vec(k) = err;              % save each mean error in a vector
        end
        [val, ind] = min(err_vec);         % minimum mean error with each signal
        if ind==i
           match_num = match_num + 1;
%         elseif val > 0.2
%            print('no match found');
        end
        
    end
    recognition_rate(i) = match_num/test_num; 
    display(match_num);
    test_num = 0;
    match_num = 0;
    display(t_files{i});
    display(recognition_rate');
end
%}












