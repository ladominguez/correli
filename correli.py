import obspy
import numpy as np
import sys
import argparse

def four1(data, nn, isign): 
    N=2**nn
    j=1
    for i in range(1, N+1, 2):
        if j > i:
            tempr = data[j]
            tempi = data[j+1]
            data[j], data[j+1] = data[i], data[i+1]
            data[i], data[i+1] = tempr, tempi
        m = N // 2
        while m >= 2 and j > m:
            j -= m
            m = m // 2
        j += m
    mmax = 2
    while N > mmax:
        istep = 2 * mmax
        theta = 6.28318530717959 / (isign * mmax)
        wtemp = np.sin(0.5 * theta)
        wpr = -2.0 * wtemp * wtemp
        wpi = np.sin(theta)
        wr = 1.0
        wi = 0.0
        for m in range(1, mmax, 2):
            for i in range(m, N+1, istep):
                j = i + mmax
                tempr = wr * data[j] - wi * data[j+1]
                tempi = wr * data[j+1] + wi * data[j]
                data[j] = data[i] - tempr
                data[j+1] = data[i+1] - tempi
                data[i] += tempr
                data[i+1] += tempi
            wr = (wtemp := wr) * wpr - wi * wpi + wr
            wi = wi * wpr + wtemp * wpi + wi
        mmax = istep

def twofft(data1, data2, fft1, fft2, N):
    c1 = 0.5  + 0.0j
    c2 = -0.0 + 0.5j

    fft1 = data1 + data2*1j
    fft1 = four1(fft1, N, 1)
    fft2(0) = np.imag(fft1(0)) + 0j
    fft1(0) = np.real(fft1(0)) + 0j

    N2 = N + 2
    for j in range(2, N+1, 1):
        h1 = c1 * (fft1(j) + np.conj(fft1(N2-j)))
        h2 = c2 * (fft1(j) - np.conj(fft1(N2-j)))
        fft1(j) = h1
        fft1(N2-j) = np.conj(h1)
        fft2(j) = h2
        fft2(N2-j) = np.conj(h2)

    return fft1, fft2

if __name__ == "__main__":
    if len(sys.argv) != 4:
        parser = argparse.ArgumentParser(description="Process SAC files.")
        parser.add_argument("sac_file1", type=str, help="First SAC file")
        parser.add_argument("sac_file2", type=str, help="Second SAC file")
        parser.add_argument("output_file", type=str, help="Output file")
        parser.add_argument("-w", type=float, default=0.3, help="Window (default: 0.3s)")

        args = parser.parse_args()

        sac_file1 = args.sac_file1
        sac_file2 = args.sac_file2
        output_file = args.output_file
        window = args.w
        print("Usage: python correli.py <sac_file1> <sac_file2> <output_file>")
        sys.exit(1)


    try:
        st1 = obspy.read(sac_file1)
        st2 = obspy.read(sac_file2)
    except Exception as e:
        print(f"Error reading SAC files: {e}")
        sys.exit(1)

    delta1 = st1[0].stats.sac.delta
    delta2 = st2[0].stats.sac.delta
    npts1 = st1[0].stats.sac.npts
    npts2 = st2[0].stats.sac.npts

    if delta1 != delta2:
        print("Error: delta values are different")
        sys.exit(1)
    delta = delta1

    if npts1 != npts2:
        npts = min(npts1, npts2)

    lim_1 = 1
    lim_2 = np.rpund(window / delta)

    ancdelta = 0.0
    aintervalo = 0.0

    while True:
        aintervalo = aintervalo + window
        lim_2 = lim_2 + 1

        vec1 = st1[0].data[lim_1:lim_2]
        vec2 = st2[0].data[lim_1:lim_2]

        aint1 = np.sum(vec1**2)
        aint2 = np.sum(vec2**2)

        lim_1 = lim_1 + 1







    

print("SAC files loaded successfully")