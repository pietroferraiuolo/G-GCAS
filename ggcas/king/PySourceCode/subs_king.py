import numpy as np

def surf_dens_ekin(fcn):
    global X, W, DW, RAPP, NFN, WWWW, ESTR1, ESTR2, csi, uk, sigekin, Sk_Sk0, R
    u1 = 0.0
    for l in range(1, NFN - 1):
        wwww = W[l]
        estr2 = W[l]
        u1 = intgau(fcn, u1)
        uk[l] = np.exp(wwww) * u1

    sum_ = 0.0
    b = 0.0
    h = 0.0
    xiR = 0.0
    for i in range(1, NFN + 1):
        xiR = csi[i]
        for k in range(i, NFN + 1):
            yy[k] = np.sqrt(csi[k]**2 - xiR**2)
            if k == i:
                sum_ = 0.0
            else:
                b = uk[k] + uk[k - 1]
                h = yy[k] - yy[k - 1]
                sum_ += b * h * 0.5

        if i == 1:
            sigekin0 = sum_
            if sum_ == 0:
                print('zero central Skin')

        sigekin[i] = sum_ * (1.0 / (1.0 + R**2))
        Sk_Sk0[i] = sigekin[i] / sigekin0

def caloric_curve():
    global X, W, DW, RAPP, NFN, W0, RHO0, xMhat, xnhat, x0Cv, xi0Cv, potch, xk, eg, ephi, xkt, egt, ephit, Etot, Vir2, et_r, Etot1
    fw0 = Etot / xMhat
    xincr = 7.0 / 1.0e4
    yEt_Mv02[1] = -2.0
    xincr2 = 1.5 / 1.0e4
    xs2_v02[1] = 0.5
    xs[1] = 0.0
    ys[1] = 0.0
    for i in range(2, 10001):
        xs2_v02[i] = xs2_v02[i - 1] + xincr2
        yEt_Mv02[i] = yEt_Mv02[i - 1] + xincr

    xsave = 0.0
    ysave = 0.0
    w0save = 0.0
    k = 0
    for i in range(1, 10001):
        xxx = 0.0
        xxx = xs2_v02[i]
        ex2 = 0.0
        yyy = 0.0
        xxxx = 0.0
        for j in range(1, 10001):
            yyy = yEt_Mv02[j]
            fw0 = Etot / xMhat
            xxxx = yyy / fw0
            if (abs((xxxx - xxx) / xxx) < abs((xs[i] - xxx) / xxx) and
                abs((xxxx - xxx) / xxx) <= 1.0e-7):
                k += 1
                xs[k] = xxxx
                ys[k] = yyy
                XRE[k] = abs((xs[k] - xxx) / xxx)
                ex2 = 1.0

    xre2 = np.zeros(k)
    xsv = np.zeros(k)
    ysv = np.zeros(k)
    for i in range(1, k + 1):
        xre2[i - 1] = xre[i]
        xsv[i - 1] = xs[i]
        ysv[i - 1] = ys[i]

    xresave = np.min(xre2)
    for i in range(1, k + 1):
        if xresave - xre2[i - 1] == 0.0:
            xsave = xsv[i - 1]
            ysave = ysv[i - 1]
            w0save = W0
            print(f'x={xsave}, y={ysave}')

def fcnqd2(t):
    global WWWW
    return np.exp(wwww - t) * (t**2.5)

def v2mean(fcn1, fcn2):
    global X, W, DW, RAPP, NFN, v2
    for l in range(1, NFN - 1):
        wwww = W[l]
        estr2 = W[l]
        qd = intgau(fcn1)
        v2_winf = qd
        qd2 = intgau(fcn2)
        v2[l] = 0.4 * (qd2 / v2_winf)
    v2[NFN] = v2[NFN - 1]

def surfdens(fcn):
    global X, W, DW, RAPP, NFN, t, phi, sigma, psi, S_S0, csi
    qd = 0.0
    for l in range(1, NFN - 1):
        wwww = W[l]
        estr2 = wwww
        qd = intgau(fcn, qd)
        psi[l] = qd

    sum_ = 0.0
    b = 0.0
    h = 0.0
    xiR = 0.0
    for i in range(1, NFN + 1):
        xiR = csi[i]
        for k in range(i, NFN + 1):
            yy[k] = np.sqrt(csi[k]**2 - xiR**2)
            if k == i:
                sum_ = 0.0
            else:
                b = psi[k] + psi[k - 1]
                h = yy[k] - yy[k - 1]
                sum_ += b * h * 0.5

        if i == 1:
            sigma0 = sum_

        sigma[i] = sum_
        S_S0[i] = sigma[i] / sigma0

def energie(f1, f3, fcn1):
    global ESTR1, ESTR2, W0, RHO0, INDGAU, PAI, X, W, DW, RAPP, NFN, csi, xk, eg, ephi, xkt, egt, ephit, Etot, Vir2, et_r, Etot1
    z = 0.0
    u1 = 0.0
    u3 = 0.0
    for m in range(1, NFN - 1):
        z = W[m]
        wwww = W[m]
        estr2 = wwww
        u1 = intgau(f1, u1)
        wwww = W[m]
        estr2 = wwww
        u3 = intgau(f3, u3)
        xk[m] = np.exp(z) * u1 * 4.0 * PAI * np.sqrt(2.0)
        eg[m] = (2.0 / 3.0) * (np.exp(z) * u1 + 0.4 * (z**(2.5))) * DW[m] * X[m] * 4.0 * PAI * np.sqrt(2.0)
        if abs(u3) > 1.0e-30:
            ephi[m] = -xk[m] - np.exp(z) * u3 * 4.0 * PAI * np.sqrt(2.0)
        else:
            ephi[m] = 0.0

        et_r[m] = xk[m] + eg[m] + ephi[m]

    xk[NFN] = 0.0
    eg[NFN] = 0.0
    ephi[NFN] = 0.0
    et_r[NFN] = 0.0

    wwww = W0
    estr2 = W0
    qd = 0.0
    qd = intgau(fcn1, qd)
    xnhat0 = qd
    rho0hat = (4.0 * np.sqrt(2.0) * PAI) * qd
    xMcap = ((9.0 / (4.0 * PAI))**(1.5)) * (1.0 / np.sqrt(rho0hat)) * xmu
    for i in range(1, NFN + 1):
        rhat[i] = np.sqrt(9.0 / (4.0 * PAI)) * (1.0 / np.sqrt(rho0hat)) * X[i]

    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    h = 0.0
    for k in range(1, NFN - 1):
        b1 = xk[k + 1] * rhat[k + 1]**2 + xk[k] * rhat[k]**2
        b2 = eg[k + 1] * rhat[k + 1]**2 + eg[k] * rhat[k]**2
        b3 = ephi[k + 1] * rhat[k + 1]**2 + ephi[k] * rhat[k]**2
        h = rhat[k + 1] - rhat[k]
        sum1 += b1 * h * 0.5
        sum2 += b2 * h * 0.5
        sum3 += b3 * h * 0.5

    xkt = 0.0
    egt = 0.0
    ephit = 0.0
    xkt = sum1 * 4.0 * PAI
    egt = sum2 * 4.0 * PAI
    ephit = sum3 * 4.0 * PAI
    Etot = sum1 + sum2 + sum3
    Etot1 = xkt + egt + ephit
    Vir = 2.0 * xkt + egt
    Vir2 = 2.0 * xkt / (-egt)

def calorespecifico(fcn1, f1, f2, f3, f4):
    global ESTR1, ESTR2, W0, RHO0, INDGAU, PAI, X, W, DW, RAPP, NFN, csi, Cv0, Cv, Cvqt, Ctot0, Ctot1, Ctot2, xMhat, xnhat0, conc
    for k in range(1, NFN + 1):
        Cv0[k] = 0.0
        Cv[k] = 0.0
        Cvqt[k] = 0.0
        xnhat[k] = 0.0

    for m in range(1, NFN - 1):
        z = 0.0
        zz = 0.0
        wwww = W[m]
        estr2 = wwww
        z = W[m]
        qd = 0.0
        u1 = 0.0
        u21 = 0.0
        u23 = 0.0
        u3 = 0.0
        u4 = 0.0
        zz = np.sqrt(z)
        qd = intgau(fcn1, qd)
        xnhat[m] = qd
        wwww = W[m]
        estr2 = wwww
        u1 = intgau(f1, u1)
        wwww = W[m]
        estr2 = wwww
        qq = 0.5
        u21 = intgau(f2, u21)

        qq = 1.5
        wwww = W[m]
        estr2 = wwww
        u23 = intgau(f2, u23)

        potch[m] = u21 * np.exp(z) / qd

        wwww = W[m]
        estr2 = wwww
        u3 = intgau(f3, u3)

        wwww = W[m]
        estr1 = 0.0
        estr2 = wwww
        u4 = intgau(f4, u4)
        Cv0[m] = ((u3 * u21 * np.exp(2.0 * z) / qd) - u4 * np.exp(z))
        Cv[m] = (Cv0[m] + X[m] * DW[m] * (qd + (2.0 / 3.0) * (z**(1.5)) * (2.5 - 1.4 * z) + 
            ((4.0 / 9.0) * (z**(4.0)) * (0.4 * z - 2.0) / qd) + 
            ((16.0 / 135.0) * (z**(6.5)) / (qd**2))))

        gw = np.exp(z) * erf(zz) * (-z + (np.sqrt(2.0) - 1.0) / 2.0) - (z + np.sqrt(2.0 * z)) / np.sqrt(PAI)
        hw = (zz / 4.0 + np.sqrt(PAI) * gw / 4.0 + (z**(1.5)) / 2.0) / (-zz / 2.0 + np.exp(z) * np.sqrt(PAI) * erf(zz) / 4.0 - (z**(1.5)) / 3.0)
        Cvqt[m] = (Cv0[m] + 0.5 * (1.5 - C + hw) * (-C - z) * xnhat[m])

    wwww = W0
    estr2 = W0
    qd = 0.0
    qd = intgau(fcn1, qd)
    xnhat0 = qd
    rho0hat = (4.0 * np.sqrt(2.0) * PAI) * qd
    xMcap = ((9.0 / (4.0 * PAI))**(1.5)) * (1.0 / np.sqrt(rho0hat)) * xmu

    qd = 0.0
    xnhat[NFN] = 0.0
    u1 = 0.0
    u21 = 0.0
    u23 = 0.0
    u3 = 0.0
    u4 = 0.0
    Cv0[NFN] = 0.0
    Cv[NFN] = 0.0
    Cvqt[NFN] = 0.0

    sum0 = 0.0
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    b0 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    h = 0.0
    for k in range(1, NFN - 1):
        b0 = (Cv0[k + 1] * csi[k + 1]**2 + Cv0[k] * csi[k]**2)
        b1 = (Cv[k + 1] * csi[k + 1]**2 + Cv[k] * csi[k]**2)
        b2 = (Cvqt[k + 1] * csi[k + 1]**2 + Cvqt[k] * csi[k]**2)
        b3 = (xnhat[k + 1] * csi[k + 1]**2 + xnhat[k] * csi[k]**2)
        h = csi[k + 1] - csi[k]
        sum0 += b0 * h * 0.5
        sum1 += b1 * h * 0.5
        sum2 += b2 * h * 0.5
        sum3 += b3 * h * 0.5

    Ctot0 = 0.0
    Ctot1 = 0.0
    Ctot2 = 0.0
    Ctot0 = sum0 / sum3
    Ctot1 = sum1 / sum3
    Ctot2 = sum2 / sum3

    xMhat = sum3
    for i in range(1, 11):
        x0Cv[i] = 0.0
        xi0Cv[i] = 0.0

    ex = 0.0
    ex2 = 0.0
    ex3 = 0.0
    for m in range(1, NFN - 1):
        if (Cv[m] >= 0.0 and Cv[m + 1] <= 0.0 and ex == 0.0):
            x0Cv[1] = (X[m] + X[m + 1]) * 0.5
            xi0Cv[1] = (csi[m] + csi[m + 1]) * 0.5
            ex = 1.0
        elif (Cv[m] <= 0.0 and Cv[m + 1] >= 0.0 and ex2 == 0.0):
            x0Cv[2] = (X[m] + X[m + 1]) * 0.5
            xi0Cv[2] = (csi[m] + csi[m + 1]) * 0.5
            ex2 = 1.0
        elif (Cv[m] >= 0.0 and Cv[m + 1] <= 0.0 and ex == 1.0 and ex3 == 0.0):
            x0Cv[3] = (X[m] + X[m + 1]) * 0.5
            xi0Cv[3] = (csi[m] + csi[m + 1]) * 0.5
            ex3 = 1.0

def funzerr(x):
    return np.exp(-(x**2))

def fcnqd(x):
    global WWWW
    return np.exp(wwww - x) * (x**1.5)

def fcnqd1(x):
    global WWWW
    return (np.exp(wwww - x) - 1.0) * (x**0.5)

def fu1(x):
    global WWWW
    return (np.exp(-x) - np.exp(-wwww)) * (x**1.5)

def fu2(x):
    global WWWW, qq
    return (x * np.exp(-x) - wwww * np.exp(-wwww)) * (x**qq)

def fu3(x):
    global WWWW
    return (np.exp(-x) - np.exp(-wwww)) * np.log(np.exp(-x) - np.exp(-wwww)) * (x**0.5)

def fu4(x):
    global WWWW
    return (x * np.exp(-x) - wwww * np.exp(-wwww)) * np.log(np.exp(-x) - np.exp(-wwww)) * (x**0.5)

def FCT(X, Y, DERY, FCN):
    global WWWW, ESTR1, ESTR2, W0, RHO0, RPP
    wwww = Y[0]
    estr2 = Y[0]

    if wwww > 0.0:
        rpp = 0.0
    else:
        rpp = 0.0

    qd = intgau(FCN)
    rpp = qd / rho0

    DERY[0] = Y[1]
    DERY[1] = -(2.0 / X) * Y[1] - 9.0 * rpp

def OUTP(XX, Y, DERY, IHLF, NDIM, PRMT):
    global X, W, DW, RAPP, NFN, RPP
    if Y[0] > 0.0:
        PRMT[4] = 1.0
    else:
        xrif = X[NFN] + PRMT[2]
        if XX < xrif:
            return

    NFN += 1
    if NFN == 8000:
        print('TROPPI PASSI')
        PRMT[4] = 1.00

    X[NFN] = XX
    W[NFN] = Y[0]
    DW[NFN] = Y[1]
    RAPP[NFN] = RPP

def INTGAU(fcn, area):
    global INDGAU
    if INDGAU == 1:
        GAUS20(fcn, area)
    elif INDGAU == 2:
        GAUS40(fcn, area)
    elif INDGAU == 3:
        GAUS80(fcn, area)
    elif INDGAU == 4:
        GAUS96(fcn, area)


