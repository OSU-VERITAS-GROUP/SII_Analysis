# Noise Analysis 

As you may or may not be aware of there is a 79.5 MHz Noise that frequently interferes with observations at the VERITAS site. It appears sinusoidally in the correlation functions and there are a number of things we can do in order to remove this noise. 

$$A*sin(\omega\tau +\phi)$$

We use sine and cosine fourier components now to remove the noise but it still takes on the sinusoidal form in any given correlation. We need to understand these parameters to be confident in our noise removal.  

## The Frequency!!

The frequency of ***the*** noise is very apparent at 79.5 MHz. You can see this in both the correlated and uncorrelated data. However, we can get the most accurate reading of the frequency using the raw data since the precision of the FFt is a function of the number of points fed into the FFT. We can take at FFT of the Raw data with ***GoodFFT.cpp***. 

## Geometric Noise Analysis (Phase)
Using the Phase Parameter we can extract important information about the location of our potential beacon. ***GeometericNoiseAnalysis.C*** contains callibation data of incramenetenting the telescopes at fitting with a point source noise model to find the location of the noise beacon. 

*More of this data lives on the CHPC but I have not had a chance to analyze it*

## Comments on the Amplitude

The amplitude (somtimes) changes smoothly for high noise old runs, but it is unlikely that there is any reason for this. 
