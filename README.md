# CoviDSP 🔊🎙️

## Final Project: Speaker recognition system

Elaborated by: 

Howard Kao - hkao [at] ucdavis [dot] edu

William Orozco - worozco [at] ucdavis [dot] edu

EEC201 - University of California, Davis. Winter Quarter 2021


---

### A. Introduction 

The purpose of the project is to build an automatic speaker recognition system. Features were extracted from input speech [files](https://github.com/worzs/CoviDSP/tree/main/src/Train), by applying th Fourier Transform to the signal, and then obtaining the Mel frequeny cepstrum coefficients (MFCC). The characteristics of an audio signal vary over time. Therefore, applying Windowing and the Short Time Fourier Transform is convenient to locate the regions with useful information, and isolate the useless sectors. Signals left are the features that are used to train and to evaluate each test speaker.

After the feature extraction, we are ready to calculate the centroids using the LBG algorithm. They are the codewords of the codebook for each speaker. Finally, we test the system by identifying the speaker in a different dataset. 

TODO: add flow diagrams
---
### B. Data preprocessing
The input signal contains 11 different individual speaking the word "Zero." Each of the sampling rate is 12.5 KHz.

 Before the feature extraction process, the signals were analyzed to obtain general characteristics like shape, amplitude, mean, noise and quiet regions. In the following figures we compare the raw  and edited signals of speaker 3. The original samples show quiet regions. They were cropped (removing all samples with amplitude lower than -10 dB at the beginning and at the end), and we also performed an amplitude normalization to be between -1 and 1. Only monophonic sound was considered. 
 
![s3?time](https://github.com/worzs/CoviDSP/blob/main/src/img/time_original_edited_s3.png)

Speakers 9, 10 and 11 have stereo recordings. In the next figure, we compare the raw data and the edited signal for speaker 10. 

![s10_time](https://github.com/worzs/CoviDSP/blob/main/src/img/time_original_edited_s10.png)

Data visualization and preprocessing is critical to improve the feature extraction stage. If we train our system with useless datasets, our results will be affected and the most important features for each speaker will be hidden. 

---
### C. Feature extraction

After the data preprocessing, the Short-Time Fourier Transform was applied to the blocks of the original signal. The size of the frame is N=256 with a Hamming window applied, and an overlap of M = 100 samples between. 
The spectrograms of the speakers were plotted. The spectral content of the signals shows higher manigtudes near low frequencies. Then, the Mel filterbank with size 20 wrapped the fourier transform of the signal blocks and we also obtained the cepstrum coefficients by applying the Discrete Cosine Transform.  The first component of the DCT was dropped, as it is not relevant to the features of the speaker. 

#### Mel filter bank frequency response:

![melfb](https://github.com/worzs/CoviDSP/blob/main/src/img/melfb.png)

Below we show the relevance of data preprocessing for this project. If we take the audio signal as it is, we will include useless data to obtain the centroids. In general, the frequency response of the speakers tends to be in low ranges. Quiet regions have an almost flat response. The Mel frequency warping compress the information, providing representative features of the speaker. The graphics to the right in the next plots are a reduced-resolution version of the ones on the left,  confirming that the extraction of the Mel cepstrum coefficients is working. 

#### Speaker 3
The next figure shows the spectrogram and MFCCs for Speaker 3. 

Top-Left: Spectrogram of the raw signal.

Top-Right: Frequency response mapped to the accoustic vectors (MFCCs) of the raw signal.

Bottom-Left: Spectrogram of the cropped,normalized signal.

Bottom-Right: Frequency response mapped to the accoustic vectors (MFCCs) of the cropped,normalized signal.

![mfcc_s3](https://github.com/worzs/CoviDSP/blob/main/src/img/spectrogram_mfcc_s3.png)


#### Speaker 10
The next figure shows the spectrogram and MFCCs for Speaker 10. 

Top-Left: Spectrogram of the raw signal.

Top-Right: Frequency response mapped to the accoustic vectors (MFCCs) of the raw signal.

Bottom-Left: Spectrogram of the cropped,normalized signal.

Bottom-Right: Frequency response mapped to the accoustic vectors (MFCCs) of the cropped,normalized signal.


![mfcc_s10](https://github.com/worzs/CoviDSP/blob/main/src/img/spectrogram_mfcc_s10.png)

---
### D. Clustering

From the previous step, we take the MFCCs of the normalized, cropped signal because they provide useful information from the active regions of the audio. Then we calculate the centroids using the LGB algorithm with K centroids, epsilon (splitting parameter) 0.01, and an error threshold (distortion) 0.001. The distance criteria is the euclidean norm (L2). 

In the next graphs, the 3rd and 5th dimensions of the centroids and the accoustic vectors of two speakers are plotted, K = 8. The codewords converge among the clusters of the speakers. 

#### Speaker 3

![codebook_features_s3](https://github.com/worzs/CoviDSP/blob/main/src/img/centroids_mfcc3_mfcc5_s3.png)

#### Speaker 10

![codebook_features_s10](https://github.com/worzs/CoviDSP/blob/main/src/img/centroids_mfcc3_mfcc5_s10.png)

---
### E. Results

####      1. Original samples. 

In the following table we summarize the results of the system when it matches speakers from the original training and test sets. From the spectrograms shown previously, the frequency content is mainly concentrated on the indexes 15 and below. Between the indexes 15 and 20, the amplitudes of the spectrum are lower, and they can be discarded to find the codewords. Parameters:
```matlab
%windowing parameters
N = 256; % window size
M = 100; % overlap
p = 20;  % number of filters in filterbank

%lbg parameters
lbg_p = 15; % length of the column vector for the lbg clustering. 
K = 32; % number of clusters

%normalization
type_signal = 'edit'; %'edit': normalized signal; 'raw': original signal
```

The accuracy is 100% or 0%, and it represents if the system recognizes the speaker from the test dataset based on the codebook created with the training dataset. It was able to recognize the 8 speakers. 

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% 

Next, we show the codebook error (based on the distances from the codebook) of each speaker. The values are not higher than 5%. 

![error_original_parameters](https://user-images.githubusercontent.com/33579806/111860293-db2d0500-8903-11eb-8f23-a2624caf9137.png)

A comparison between the centroids obtained with the training data (dimensions 3 and 5), and the test dataset fit. 

Speaker 8 Training | Speaker 8 Test
--- | --- 
![centroids_mfcc3_mfcc5_s8_k32](https://user-images.githubusercontent.com/33579806/111860649-29430800-8906-11eb-88a2-2508bb971ed9.png) | ![centroids_mfcc3_mfcc5_s8_k32_test](https://user-images.githubusercontent.com/33579806/111860669-48da3080-8906-11eb-9e3a-5509818252d5.png)



-----------------------------
We modified the parameters to determine their relevance. 
First, we tried with the original signal without the normalization step. 

```matlab
%windowing parameters
N = 256; % window size
M = 100; % overlap
p = 20;  % number of filters in filterbank

%lbg parameters
lbg_p = 15; % length of the column vector for the lbg clustering. 
K = 32; % number of clusters

%normalization
type_signal = 'raw'; %'edit': normalized signal; 'raw': original signal
```

The system is still able to recognize the all the speakers. 

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% 

The error is almost constant, between 4% and 5%, for each speaker. 

![error_original_parameters_raw](https://user-images.githubusercontent.com/33579806/111860861-8be8d380-8907-11eb-89f9-dc8741b9b531.png)

Again, we made a comparison between the centroids obtained with the training data (dimensions 3 and 5), and the test dataset fit. It is clear that the clusters are closer when the normalization is not performed.

Speaker 8 Training | Speaker 8 Test
--- | --- 
![centroids_mfcc3_mfcc5_s8_k32_raw](https://user-images.githubusercontent.com/33579806/111860948-eda93d80-8907-11eb-9b15-2729fa677b38.png) | ![centroids_mfcc3_mfcc5_s8_k32_raw_test](https://user-images.githubusercontent.com/33579806/111860959-f8fc6900-8907-11eb-885c-43849897de55.png)



-----------------------------
Now we tried with a lower window size, K= 2 clusters, using the first 5 elements in a 10-element mel filterbank, and including the normalization step for the signal. 

```matlab
%windowing parameters
N = 100; % window size
M = 30; % overlap
p = 10;  % number of filters in filterbank

%lbg parameters
lbg_p = 5; % length of the column vector for the lbg clustering. 
K = 2; % number of clusters
```

This time, we compared the behavior of the system with the parameters listed above, when the input signal is normalized and when it is not. 
```matlab
%normalization
type_signal = 'edit'; %'edit': normalized signal; 'raw': original signal
```


Without any upper threshold for the tolerable error (maximum distance from the codebook), our system is still able to identify the speakers. The codebook error increased to be around 30% for the speakers. As we said before, the accuracy is still 100% without distance limitations. However, if we set our system to tolerate a maximum of 10% of error, the accuracy goes to 0% for all cases. 

![error_edit_30_percent](https://user-images.githubusercontent.com/33579806/111861344-b8521f00-890a-11eb-8d96-b92576a30d37.png)
 

The next graphs show the relevance of the normalization step when the system has reduced characteristics. In other words, a few filters in the filter bank, a smaller window size and a couple of clusters. 
```matlab
%normalization
type_signal = 'raw'; %'edit': normalized signal; 'raw': original signal
```

Speaker 4 Test edited signal | Speaker 4 original signal
--- | --- 
![centroids_mfcc3_mfcc5_s4_k2](https://user-images.githubusercontent.com/33579806/111861321-6e693900-890a-11eb-8691-efe5a462aaf7.png)  |  ![centroids_mfcc3_mfcc5_s4_k2_raw](https://user-images.githubusercontent.com/33579806/111861449-842b2e00-890b-11eb-9edc-9f30b3ac0e23.png)

What is Happening? With a "poor" system, spacing the clusters critical. If they are too close, the speaker identification fails because the codewords are close enough to create  false positives. The normalization of the signal helps here. 

The following table summarizes the accuracy of the system. Speakers 1 and 6 are not matched correctly. 

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 0% | 100% | 100% | 100% | 100% | 0% | 100% | 100% 



---
#### 2. Original samples with noise added. 
```matlab
%windowing parameters
N = 256; % window size
M = 100; % overlap
p = 20;  % number of filters in filterbank

%lbg parameters
lbg_p = 15; % length of the column vector for the lbg clustering. 
K = 32; % number of clusters

%normalization
type_signal = 'edit'; %'edit': normalized signal; 'raw': original signal
```
We use the function randn in MATLAB to generate the noise to the signals. The randn function generates pseudo-random numbers whose elements are normally distributed with mean 0 and variance 1 (standard normal). Which is also known as the function to generate Gaussian distributed variables. We added 4 different variances to generate 4 degree of noise which is 0.001, 0.0035, 0.006, and 0.013 corresponds to the signal-to-noise ratio (SNR) 25 dB, 15dB, 10dB, and 5dB, respectively. The figure is an example of the signals of speaker 4 and speaker 6 with and without noise added. The figure shows that the noise is added with a variance of 0.013 (~SNR = 5dB), which is a level that our system starts to lose accuracy from recognition.

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% 

With 0.0035 as the variance of the noise (~SNR = 15 dB)

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% 

With 0.0060 as the variance of the noise (~SNR = 10 dB)

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 100% | 100% |  0%  | 100% |  0%  | 100% | 0% 

With 0.013 as the variance of the noise (~SNR = 5 dB)

![S4_orignial]https://github.com/worzs/CoviDSP/blob/main/src/img/S4_original.jpg
![S4_with_noise]
![S6_orignial]
![S6_with_noise]

Table X: Results 

---
#### 3. Original samples with notch filters at different frequencies. 


We modified the parameters of the filter bank size, clusters and we also set the signal type to be normalized. 
```matlab
%windowing parameters
N = 256; % window size
M = 100; % overlap
p = 16;  % number of filters in filterbank

%lbg parameters
lbg_p = 15; % length of the column vector for the lbg clustering. 
K = 16; % number of clusters

%normalization
type_signal = 'edit'; %'edit': normalized signal; 'raw': original signal
```

Our first notch filter was centered at 1 KHz, with a bandwidth of 400 Hz. 
Speaker 4 frequency response | Speaker 4 spectrogram
--- | --- 
![s4_notch_1k_bw_400](https://user-images.githubusercontent.com/33579806/111867789-8dc88c00-8933-11eb-851c-2f65b4f286f4.png) | ![s4_notch_1k_bw_400_spectrogram](https://user-images.githubusercontent.com/33579806/111867790-902ae600-8933-11eb-89e7-b02c76fa968f.png)


We got an accuracy of 100% for all the speakers. If we see in the spectrogram, the majority of spectral content is located in lower frequencies. 

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% 


Now we tried with a notch filter centered at 500 Hz, with a bandwidth of 500 Hz. The accuracy decreased to 6/8, two speakers were not matched. 

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 0% | 100% | 100% | 100% | 0% | 100% | 100% 

If we examine the frequency response and the spectrogram of speaker 6, we notice that an important slice of data was truncated by the notch filter. In low frequencies there is essential spectral content. If we remove it, the system will not recognize the speaker. 

Speaker 6 frequency response | Speaker 6 spectrogram
--- | --- 
![s6_notch_500_bw_500](https://user-images.githubusercontent.com/33579806/111868032-ea787680-8934-11eb-96fb-50e772c9ed3e.png) | ![s6_notch_500_bw_500_spectrogram](https://user-images.githubusercontent.com/33579806/111868033-ec423a00-8934-11eb-98e6-9f9da8b66f29.png)

A similar behavior is observed in the spectrogram and frequency response of speaker 2, the second that was not recognized in the test set. 

Speaker 2 frequency response | Speaker 2 spectrogram
--- | --- 
![s2_notch_500_bw_500](https://user-images.githubusercontent.com/33579806/111868124-5e1a8380-8935-11eb-8eef-53b8726ad751.png) | ![s2_notch_500_bw_500_spectrogram](https://user-images.githubusercontent.com/33579806/111868127-61157400-8935-11eb-9c53-9ac70ba5e0bc.png)

In conclusion, different speakers have their own features in the frequency domain. It is possible to fool a speaker recognition system by removing the main frequencies of the speaker, applying notch filters. For the provided samples, the most relevant spectral content was located in lower frequencies (around 500 Hz). 

---
#### 4. Testing with different audio signals 

TODO

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Accuracy | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

### F. Instructions to run the code

Start by cloning the repository.

CoviDSP1.m generates the codebooks and test the recognition based on the samples from all t.

You can modify the parameters at the top of CoviDSP1.m before running it, to explore the clustering and accuracy of the system.  Next, main_notch.m trains and tests the system with notch filters applied. You can modify the same parameters, including the notch frequency and the bandwidth of the filter. 
 
 The most relevant:
```matlab
%STFT parameters
N = 256; % window size
M = 100; % overlap

%Mel Filter Bank
p = 20;  % number of filters in filterbank

%LBG parameters
lbg_p = 15; % the first lbg_p filters from the original p will be used for the clustering
K = 8; % number of clusters
error_thresh = 0.05; 
```



### G. References

Y. Linde, A. Buzo & R. Gray, “An algorithm for vector quantizer design”, IEEE Transactions on Communications, Vol.
28, pp.84-95, 1980

Matlab documentation [Online]. In: https://www.mathworks.com/help/matlab/ Access: March, 2021. 



