import obspy
import numpy as np
import sys
import argparse

alpha = 8313.0
beta = 4800.0

def correl(data1, data2, N):
    fft1, fft2 = twofft(data1, data2, N)
    NO2 = N // 2
    fft2 - fft1*np.conj(fft2)/NO2
    fft2(0) = np.real(fft2(0)) + np.real(fft2(N/2+1))
    return realft(fft2, N/2, -1)
    

def realft(data, N, isign):
    theta = 3.141592653589793238 / N
    c1 = 0.5
    c2 = -0.5 * isign
    if isign == 1:
        four1(data, N, 1)
    else:
        theta = -theta
    wpr = -2.0 * np.sin(0.5 * theta)**2
    wpi = np.sin(theta)
    wr = 1.0 + wpr
    wi = wpi
    n2p3 = 2*N + 3
    for i in range(1, N//2+1, 1):
        i1 = 2*i - 1
        i2 = i1 + 1
        i3 = n2p3 - i2
        i4 = i3 + 1
        h1r = c1 * (data[i1] + data[i3])
        h1i = c1 * (data[i2] - data[i4])
        h2r = -c2 * (data[i2] + data[i4])
        h2i = c2 * (data[i1] - data[i3])
        data[i1] = h1r + wr*h2r - wi*h2i
        data[i2] = h1i + wr*h2i + wi*h2r
        data[i3] = h1r - wr*h2r + wi*h2i
        data[i4] = -h1i + wr*h2i + wi*h2r
        wtemp = wr
        wr = wr*wpr - wi*wpi + wr
        wi = wi*wpr + wtemp*wpi + wi
    if isign == 1:
        data[0] = (h1r := data[0]) + data[1]
        data[1] = h1r - data[1]

    else:
        data[0] = c1 * ((h1r := data[0]) + data[1])
        data[1] = c1 * (h1r - data[1])
        four1(data, N, -1)
    return data

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
        theta = 2*np.pi / (isign * mmax)
        
        wpr = -2.0 * np.sin(0.5 * theta)
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
    lim_2 = np.round(window / delta)

    ancdelta = 0.0
    aintervalo = 0.0

    noise = True
    while True:
        aintervalo = aintervalo + window
        lim_2 = lim_2 + 1

        vec1 = st1[0].data[lim_1:lim_2]
        vec2 = st2[0].data[lim_1:lim_2]

        aint1 = np.sum(vec1**2)
        aint2 = np.sum(vec2**2)

        lim_1 = lim_1 + 1

        cor = correl(vec1, vec2, len(vec1))

        k = (npts/4) + 1
        res = np.zeros(npts/2)
        for i in range(0, npts/4, 1):
            res[k] = np.real(cor[i])
            k += 1
        
        for i in range(0, npts/4, 1):
            res[i] = np.real(cor[k])
            k += 1

        corrim = -(delta * npts/2)
        delta1 = delta*2
        fac = (aint1*aint2)**0.5

        amax = np.max(res)
        amaxc = amax/fac

        u1 = np.sum(vec1[lim_1:lim_2+1]**2)
        u2 = np.sum(vec2[lim_1:lim_2+1]**2)
        u1 = u1/(lim_2 - lim_1 + 1)
        u2 = u2/(lim_2 - lim_1 + 1)

        if noise:
            an1 = np.sum(vec1[lim_1:lim_2+1]**2)
            an2 = np.sum(vec2[lim_1:lim_2+1]**2)
            an1 = an1/(lim_2 - lim_1 + 1)
            an2 = an2/(lim_2 - lim_1 + 1)
            noise = False

        anoisec = np.sqrt(1 - (an1/u1))*np.sqrt(1 - (an2/u2))
        amaxc = amaxc/anoisec

        amax = np.clip(amax, None, 1)

        sigma = 2*(1 - amaxc)/355.3
        anum = 7*((2/alpha**6) + (3/beta**6))
        aden = (6/alpha**8) + (7/beta**8)
        dist = np.sqrt(sigma*anum/aden)
        ancdelta = ancdelta + delta

        













    

print("SAC files loaded successfully")
