from __future__ import division
import sys, os, warnings

from pylab import rfft, find
from numpy import argmax, log, diff
from scipy.io import wavfile
from scipy.signal import fftconvolve
from scipy.signal.windows import blackmanharris

def parabolic(f, x):
    xv = 1/2 * (f[x-1] - f[x+1]) / (f[x-1] - 2 * f[x] + f[x+1]) + x
    yv = f[x] - 1/4 * (f[x-1] - f[x+1]) * (xv - x)
    return (xv, yv)



#This function preprocesses the given signal to be a 1-dimensional array of floats, normalized
def prepare(fs, signal):
    if(len(signal.shape) > 1):
        signal = signal[:, 0]
    signal = signal.astype(float) / 2**16
    fs = float(fs)
    return fs, signal


def freq_from_fft(sig, fs):
    windowed = sig * blackmanharris(len(sig))
    f = rfft(windowed)
    
    i = argmax(abs(f)) # Just use this for less-accurate, naive version
    true_i = parabolic(log(abs(f)), i)[0]
    
    return fs * true_i / len(windowed)


def autocorrelation(sig, fs):

    corr = fftconvolve(sig, sig[::-1], mode='full')
    corr = corr[len(corr)//2:]
    
    d = diff(corr)
    start = find(d > 0)[0]
    peak = argmax(corr[start:]) + start
    px, _ = parabolic(corr, peak)
    
    return fs / px


def main():
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    
    ok_counter = 0
    for path in sys.argv[1:]:
        try:
            fs , wf = wavfile.read(path)
            fs, signal = prepare(fs, wf)
            autocorr = autocorrelation(signal, fs)
            ffs = freq_from_fft(signal, fs)

            if(autocorr > 180 and ffs > 180):
                out = "K"
            elif(autocorr > 180 and ffs < 10):
                out = "M"
            elif(autocorr <= 180):
                out = "M"
            else:
                out = "K"    
        except:
            pass

if __name__ == "__main__":
    main()