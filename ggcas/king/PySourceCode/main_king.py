import numpy as np

def main():
    temp1 = 0.0
    temp2 = 0.0
    second(temp1)

    print('IW1,IW2,IWSTEP  1,200,1  W0 = PARW0(IW)')
    IW1 = 1
    IW2 = 6
    IWSTEP = 1

    indgau = 4
    relerr = 1e-9
    estr1 = 0.0
    pai = 3.141592654
    w0 = 6.19
    w = np.zeros(8000)
    dw = np.zeros(8000)
    w[0] = w0
    dw[0] = 0.0
    wwww = w0
    estr2 = w0
    qd = prova(funzione, "PROVA", 2.0)
    rho0 = qd
    rapp = np.zeros(8000)
    rapp[0] = qd / rho0

    AAA = -1.5
    XMIN = (-relerr * w0 / AAA) ** 0.5
    IIII = np.log10(XMIN)
    if IIII < 0:
        IIII = IIII - 1
    XMIN = 10.0 ** IIII
    X = np.zeros(8000)
    X[1] = XMIN
    w[1] = w0 + AAA * XMIN ** 2.0
    dw[1] = 2.0 * AAA * XMIN
    wwww = w[1]
    estr2 = w[1]
    qd = intgau(fcnqd, qd)
    rapp[1] = qd / rho0

    nfn = 2
    passo = X[nfn] / 10.0
    xmax = X[nfn] * 10.0

    while True:
        ABSERR = relerr * w[nfn]
        if ABSERR < 1e-12:
            ABSERR = 1e-12
        PRMT = np.zeros(5)
        PRMT[0] = X[nfn]
        PRMT[1] = xmax
        PRMT[2] = passo
        PRMT[3] = ABSERR
        y = np.zeros(2)
        dery = np.zeros(2)
        y[0] = w[nfn]
        y[1] = dw[nfn]
        dery[0] = 0.5
        dery[1] = 0.5
        IHLF = 0  # Placeholder for actual value
        DHPCG(PRMT, y, dery, 2, IHLF, FCT, OUTP, AUX, FCNQD)

        if IHLF < 11:
            print('bisezioni>10', nfn, w0, X[nfn])
            passo /= 10.0
            xmax = X[nfn] + (xmax - X[nfn]) / 10.0
            continue

        if nfn == 8000:
            break
        if w[nfn] < 0.0:
            continue

        passo *= 10.0
        xmax *= 10.0

        xtest = (dw[nfn - 1] * X[nfn - 1] - w[nfn - 1]) / dw[nfn - 1]
        if abs((xtest - X[nfn]) / xtest) < 1e-9:
            break
        xmax = X[nfn]
        passo = (X[nfn] - X[nfn - 1]) / 10.0
        nfn -= 1

    conc = (xtest + X[nfn]) / 2.0
    csi = np.zeros(nfn)
    for kk in range(nfn):
        csi[kk] = X[kk] / conc
        if kk == nfn - 1:
            csi[kk] = conc / conc

    xlogc = np.log10(conc)
    xmu1 = -(4.0 * pai / 9.0) * (X[nfn - 1] ** 2.0) * dw[nfn - 1]
    xmu2 = -(4.0 * pai / 9.0) * (X[nfn] ** 2.0) * dw[nfn]
    xmu = (xmu1 + xmu2) / 2.0
    C1 = -(X[nfn - 1] * dw[nfn - 1])
    C2 = -(X[nfn] * dw[nfn])
    C = (C1 + C2) / 2.0
    D = C + w0
    surfdens(fcnqd)
    v2mean(fcnqd, fcnqd2)
    calorespecifico(fcnqd1, fu1, fu2, fu3, fu4)
    xKb = 1.0
    Ct0_NtK = xKb * Ctot0
    Ct1_NtK = xKb * Ctot1
    Ct2_NtK = xKb * Ctot2
    energie(fu1, fu3, fcnqd1)

    if w0 == 8.20:
        xm_k = 0.54274
        for k in range(nfn):
            phi_g[k] = (-C - w[k]) / xm_k

    xmu_r = np.zeros(nfn)
    for l in range(nfn - 1):
        b1 = 0.0
        h1 = 0.0
        sum1 = 0.0
        if l > 1:
            for i in range(2, l + 1):
                b1 = rapp[i - 1] * X[i - 1] ** 2 + rapp[i] * X[i] ** 2
                h1 = X[i] - X[i - 1]
                sum1 += b1 * h1 * 0.5
        xmu_r[l] = sum1

        b2 = 0.0
        h2 = 0.0
        sum2 = 0.0
        for m in range(l, nfn - 1):
            b2 = rapp[m + 1] * X[m + 1] + rapp[m] * X[m]
            h2 = X[m + 1] - X[m]
            sum2 += b2 * h2 * 0.5
        secpezz[l] = sum2

        phi_g2[l] = -(9.0 / xm_k) * ((sum1 / X[l]) + sum2)

    for k in range(2, nfn - 1):
        eg2[k] = xnhat[k] * (-9.0) * (xmu_r[k] / X[k]) + secpezz[k]
    eg2[0] = 0.0
    eg2[nfn] = 0.0
    b = 0.0
    h = 0.0
    sum = 0.0
    for l in range(nfn - 1):
        b = eg2[l + 1] * csi[l + 1] ** 2 + eg2[l] * csi[l] ** 2
        h = csi[l + 1] - csi[l]
        sum += b * h * 0.5
    egt2 = sum * 4.0 * pai
    Shat = IW * 0.01
    
    with open('params.dat', 'a') as f:
        f.write(f"{w0} {conc} {xlogc} {xmu}\n")
    
    with open('profiles.dat', 'a') as f:
        f.write(f"x xi w rho/rho0 SD v2 logc={xlogc}\n")
        for jj in range(nfn):
            f.write(f"{X[jj]} {csi[jj]} {w[jj]} {rapp[jj]} {s_s0[jj]} {v2[jj]}\n")

    if w0 in [0.8, 1.35, 2.0, 3.5, 5.0, 6.0, 8.0, 8.2, 40.0]:
        with open('Cv.dat', 'a') as f:
            f.write(f"x xi w  dw Cv0 Cv Cvqt w0={w0}\n")
            for m in range(nfn):
                f.write(f"{X[m]} {csi[m]} {w[m]} {dw[m]} {Cv0[m] / xnhat[m]} {Cv[m] / xnhat[m]} {Cvqt[m] / xnhat[m]}\n")

    with open('CvNtK.dat', 'a') as f:
        f.write(f"{w0} {Ct0_NtK} {Ct1_NtK} {Ct2_NtK} {conc} {x0Cv} {xi0Cv}\n")

    with open('x0Cv.dat', 'a') as f:
        f.write(f"{w0} {x0Cv[0]} {x0Cv[1]} {x0Cv[2]} {xi0Cv[0]} {xi0Cv[1]} {xi0Cv[2]} {conc} {xlogc}\n")

    if w0 == 8.20:
        with open('phi.dat', 'a') as f:
            for k in range(nfn):
                f.write(f"{X[k]} {csi[k]} {w[k]} {dw[k]} {phi_g[k]}\n")

    espo = 4.0 / 3.0
    s2_v02 = (1.0 / (xMcap ** espo))
    Et_Mv02 = s2_v02 * (Etot1 / xMcap)
    xkt = s2_v02 * (xkt / xMcap)
    egt = s2_v02 * (egt / xMcap)
    ephit = s2_v02 * (ephit / xMcap)

    with open('Er.dat', 'a') as f:
        pass  # Write appropriate data

    with open('Etot.dat', 'a') as f:
        if w0 in [0.8, 1.35, 2.0, 3.5, 5.0, 6.0, 8.0, 8.2, 40.0]:
            for k in range(nfn):
                f.write(f"{X[k]} {csi[k]} {w[k]} {xk[k]} {eg[k]} {ephi[k]} {et_r[k]}\n")

        f.write(f"{w0} {xkt} {egt} {ephit} {Etot} {s2_v02} {Et_Mv02} {xMcap} {Vir2} {conc} {Etot1}\n")

    surf_dens_ekin(fu1)

    with open('Skin.dat', 'a') as f:
        f.write(f'NameID={NAMESID[iw]}\n')
        f.write(f'w0={w0}\n')
        f.write('x xi w Sk sk2\n')
        for jj in range(nfn):
            f.write(f"{X[jj]} {csi[jj]} {w[jj]} {Sk_Sk0[jj]} {sigekin[jj]}\n")

    T_imp = temp2 - temp1
    print('Tempo programma(s)=', T_imp, 'T(min)=', T_imp / 60.0)
    for j in range(2):
        print(chr(7))
        sleep(1)

if __name__ == "__main__":
    main()


