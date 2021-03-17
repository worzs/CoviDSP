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
 
![image](https://github.com/worzs/CoviDSP/blob/main/src/img/time_original_edited_s3.png)

Speakers 9, 10 and 11 have stereo recordings. In the next figure, we compare the raw data and the edited signal for speaker 10. 

![image](https://github.com/worzs/CoviDSP/blob/main/src/img/time_original_edited_s10.png)

Data visualization and preprocessing is critical to improve the feature extraction stage. If we train our system with useless datasets, our results will be affected and the most important features for each speaker will be hidden.

---
### C. Results

---
#####      1. Original samples. 


Figure X: VQ for Speaker 3. 

Speaker | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Accuracy | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

Table X: Results

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

---

### D. Instructions to run the code

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
