# CoviDSP 🔊🎙️

## Final Project: Speaker recognition system

Elaborated by: 

Howard Kao - hkao [at] ucdavis [dot] edu

William Orozco - worozco [at] ucdavis [dot] edu

EEC201 - University of California, Davis. Winter Quarter 2021


---

### A. Introduction 

The purpose of the project is to build an automatic speaker recognition system. Features were extracted from input speech [files](https://github.com/worzs/CoviDSP/tree/main/src/Train), by applying th Fourier Transform to the signal, and then obtaining the Mel frequeny cepstrum coefficients (MFCC). The characteristics of an audio signal vary over time. Therefore, applying Windowing and the Short Time Fourier Transform is convenient to locate the regions with useful information, and isolate the useless sectors. 

After the feature extraction, we are ready to calculate the centroids using the LBG algorithm. They are the codewords of the codebook for each speaker. Finally, we test the system by identifying the speaker in a different dataset. 

TODO: add flow diagrams
---
### B. Data preprocessing

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

From the previous step, we take the MFCCs of the normalized, cropped signal because they provide useful information from the active regions of the audio. Then we calculate the centroids using the LGB algorithm with K = 8 centroids, epsilon (splitting parameter) 0.01, and an error threshold (distortion) 0.001. The distance criteria is the euclidean norm (L2). 

In the next graphs, the 3rd and 5th dimensions of the centroids and the accoustic vectors of two speakers are plotted. The codewords converge among the clusters of the speakers. 

#### Speaker 3

![codebook_features_s3](https://github.com/worzs/CoviDSP/blob/main/src/img/centroids_mfcc3_mfcc5_s3.png)

#### Speaker 10

![codebook_features_s10](https://github.com/worzs/CoviDSP/blob/main/src/img/centroids_mfcc3_mfcc5_s10.png)

---
### E. Results

#####      1. Original samples. 

Accuracy results with the Test dataset. 

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 
--- | --- | --- | --- |--- |--- |--- |--- |--- 
Accuracy | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% 

Our system is able to match the 8 speakers from the training and test datasets. 

---
##### 2. Original samples with noise added. 

TODO

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Accuracy | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

Table X: Results 

---
##### 3. Original samples with notch filters at different frequencies. 

TODO

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Accuracy | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

---
##### 4. Testing with different audio signals 

TODO

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Accuracy | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

### F. Instructions to run the code

Start by cloning the repository.

CoviDSP1.m generates the codebooks and test the recognition based on the samples from all t. 

You can modify the parameters at the top of CoviDSP1.m before running it, to explore the clustering and accuracy of the system. 
 
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
