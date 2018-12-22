
import sys
import scipy.io.wavfile as wvf

def main():
    path = sys.argv[1]
    sampling_rate, data = wvf.read(path)
    
    return




if __name__ == "__main__":
    main()