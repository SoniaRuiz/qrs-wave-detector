
# qrs-wave-detector
Shiny App to detect the QRS complex.

## Introduction

This Shiny App has been developed to detect QRS complexes, and P&T waves in electrocardiogram signaling samples.

## Methodology

The algorithm proposed by Pan and Tompkins has been implemented to detect QRS complexes. It has also been implemented the algorithm proposed by Mohamed Elgendi, Bjoern Eskofier, and Derek Abbott for the detection of P and T waves.

This app uses 20 electrocardiogram signaling samples randomly downloaded from the [ECG-ID Database (ecgiddb)](https://physionet.org/content/ecgiddb/1.0.0/) available on the Physionet platform, Massachusetts Institute of Technology (MIT).

## Results

### Person1

![Person1](https://github.com/SoniaRuiz/qrs-wave-detector/blob/main/results/person1.PNG "Person1")

### Person2

![Person1](https://github.com/SoniaRuiz/qrs-wave-detector/blob/main/results/person2.PNG "Person2")

See 'results' folder for more images.

# Frameworks used

This code has been partially addapted from Matlab to R using the functions below:
* [Pan Tomkin](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/45840/versions/11/previews/pan_tompkin.m/index.html)
* [Modified form of Pan Tomkin](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/50953/versions/2/previews/pt.m/index.html)

