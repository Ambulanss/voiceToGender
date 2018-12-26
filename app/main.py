from __future__ import division
import sys, os, warnings

from pylab import *
from numpy import *
from scipy.io import wavfile
from scipy.signal import fftconvolve
from scipy.signal.windows import blackmanharris

def parabolic(f, x):
    """Quadratic interpolation for estimating the true position of an
    inter-sample maximum when nearby samples are known.
   
    f is a vector and x is an index for that vector.
   
    Returns (vx, vy), the coordinates of the vertex of a parabola that goes
    through point x and its two neighbors.
   
    Example:
    Defining a vector f with a local maximum at index 3 (= 6), find local
    maximum if points 2, 3, and 4 actually defined a parabola.
   
    In [3]: f = [2, 3, 1, 6, 4, 2, 3, 1]
   
    In [4]: parabolic(f, argmax(f))
    Out[4]: (3.2142857142857144, 6.1607142857142856)
   
    """
    # Requires real division.  Insert float() somewhere to force it?
    xv = 1/2 * (f[x-1] - f[x+1]) / (f[x-1] - 2 * f[x] + f[x+1]) + x
    yv = f[x] - 1/4 * (f[x-1] - f[x+1]) * (xv - x)
    return (xv, yv)

if not sys.warnoptions:
    warnings.simplefilter("ignore")

#This function preprocesses the given signal to be a 1-dimensional array of floats, normalized
def prepare(fs, signal):
    if(len(signal.shape) > 1):
        signal = signal[:, 0]
    signal = signal.astype(float) / 2**16
    fs = float(fs)
    return fs, signal


def freq_from_fft(sig, fs):
    """Estimate frequency from peak of FFT
    
    Pros: Accurate, usually even more so than zero crossing counter 
    (1000.000004 Hz for 1000 Hz, for instance).  Due to parabolic interpolation 
    being a very good fit for windowed log FFT peaks?
    https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
    Accuracy also increases with signal length
    
    Cons: Doesn't find the right value if harmonics are stronger than 
    fundamental, which is common.
    
    """
    # Compute Fourier transform of windowed signal
    windowed = sig * blackmanharris(len(sig))
    f = rfft(windowed)
    
    # Find the peak and interpolate to get a more accurate peak
    i = argmax(abs(f)) # Just use this for less-accurate, naive version
    true_i = parabolic(log(abs(f)), i)[0]
    
    # Convert to equivalent frequency
    return fs * true_i / len(windowed)


def freq_from_autocorr(sig, fs):
    """Estimate frequency using autocorrelation
    
    Pros: Best method for finding the true fundamental of any repeating wave, 
    even with strong harmonics or completely missing fundamental
    
    Cons: Not as accurate, currently has trouble with finding the true peak
    
    """
    # Calculate circular autocorrelation (same thing as convolution, but with 
    # one input reversed in time), and throw away the negative lags
    corr = fftconvolve(sig, sig[::-1], mode='full')
    corr = corr[len(corr)//2:]
    
    # Find the first low point
    d = diff(corr)
    start = find(d > 0)[0]
    
    # Find the next peak after the low point (other than 0 lag).  This bit is 
    # not reliable for long signals, due to the desired peak occurring between 
    # samples, and other peaks appearing higher.
    # Should use a weighting function to de-emphasize the peaks at longer lags.
    # Also could zero-pad before doing circular autocorrelation.
    peak = argmax(corr[start:]) + start
    px, py = parabolic(corr, peak)
    
    return fs / px

def freq_from_crossings(sig, fs):
    """Estimate frequency by counting zero crossings
    
    Pros: Fast, accurate (increasing with signal length).  Works well for long 
    low-noise sines, square, triangle, etc.
    
    Cons: Doesn't work if there are multiple zero crossings per cycle, 
    low-frequency baseline shift, noise, etc.
    
    """
    # Find all indices right before a rising-edge zero crossing
    indices = find((sig[1:] >= 0) & (sig[:-1] < 0))
    
    # Naive (Measures 1000.185 Hz for 1000 Hz, for instance)
    #crossings = indices
    
    # More accurate, using linear interpolation to find intersample 
    # zero-crossings (Measures 1000.000129 Hz for 1000 Hz, for instance)
    crossings = [i - sig[i] / (sig[i+1] - sig[i]) for i in indices]
    
    # Some other interpolation based on neighboring points might be better. Spline, cubic, whatever
    
    return fs / mean(diff(crossings))

def main():
    for path in sys.argv[1:]:
        fs , wf = wavfile.read(path)
        fs, signal = prepare(fs, wf)
        autocorr = freq_from_autocorr(signal, fs)
        ffs = freq_from_fft(signal, fs)
        print("corr: ", freq_from_autocorr(signal, fs))
        print("ffs: ", freq_from_fft(signal, fs))
        

        
    




if __name__ == "__main__":
    main()