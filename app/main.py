
import sys, os, wave

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.fftpack import fft
from scipy import signal
from scipy.signal import argrelextrema, fftconvolve

   

def main():
    for path in sys.argv[1:]:
        wf = wave.open(path, 'r')
        signal = wf.readframes(-1)
        signal = np.fromstring(signal, 'Int16')
        rate = wf.getframerate()
        swidth = wf.getsampwidth()
        signal = signal - np.mean(signal)
        corr = fftconvolve(signal, signal[::-1], mode="full")
        corr = corr[len(corr)//2:]
        diff = np.diff(corr)
        n = [i for i in range(0,len(diff)) if diff[i]>0][0]
        peak = np.argmax(corr[n:]) + n
        if(rate/peak < 200):
            print('M')
        else:
            print('K')
        print(rate/peak)
    




if __name__ == "__main__":
    main()