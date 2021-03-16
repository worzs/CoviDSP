# CoviDSP

## Final Project: Speaker recognition system

Elaborated by: 

Howard Kao - hkao [at] ucdavis [dot] edu

William Orozco - worozco [at] ucdavis [dot] edu

EEC201 - University of California, Davis. Winter Quarter 2021


---
---
### A. Introduction 

TODO

---
---
### B. Results

---
##### 1. Original samples. 


![image](https://user-images.githubusercontent.com/33579806/111231274-3ed2cd80-85a6-11eb-8292-277afc8f1f41.png)

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
---
### C. Instructions to run the code

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
