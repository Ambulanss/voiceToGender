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


def prepare(fs, signal):
    if(len(signal.shape) > 1):
        signal = signal[:, 0]
    signal = signal.astype(float) / 2**16
    fs = float(fs)
    return fs, signal


def autocorrelation(signal, fs):

    correlation = fftconvolve(signal, signal[::-1], mode='full')
    correlation = correlation[len(correlation)//2:]
    
    d = diff(correlation)
    start = find(d > 0)[0]
    peak = argmax(correlation[start:]) + start
    px, _ = parabolic(correlation, peak)
    
    return fs / px


def main():
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    
    for path in sys.argv[1:]:
        try:
            fs , wf = wavfile.read(path)
            fs, signal = prepare(fs, wf)
            autocorr = autocorrelation(signal, fs)
            if(autocorr > 180):
                out = "K"
            else:
                out = "M"  
            print(out)
        except:
            pass

if __name__ == "__main__":
    main()