"""
Created on May 2024
    -Author: P.Ferraiuolo
"""

import numpy as np

# Define global variables
WWWW = None
ESTR1 = None
ESTR2 = None
W0 = None
RHO0 = 1.0  # Example initialization
PAI = np.pi
RPP = 0.0
NFN = 0
X = np.zeros(8000)
W = np.zeros(8000)
DW = np.zeros(8000)
RAPP = np.zeros(8000)
indgau = 1  # Example initialization, needs to be set appropriately

PARW0 = [
    1.00e-2, 3.00e-2, 5.00e-2, 7.00e-2, 9.00e-2,
    1.00e-1, 2.00e-1, 3.00e-1, 4.00e-1, 5.00e-1,
    6.00e-1, 7.00e-1, 8.00e-1, 9.00e-1, 1.00,
    1.10, 1.20, 1.30, 1.31, 1.32, 1.33, 1.34,
    1.35, 1.36, 1.37, 1.38, 1.39, 1.40, 1.50,
    1.60, 1.70, 1.80, 1.90, 2.00,
    2.10, 2.20, 2.30, 2.40, 2.50,
    2.60, 2.70, 2.80, 2.90, 3.00,
    3.10, 3.20, 3.30, 3.40, 3.50,
    3.60, 3.70, 3.80, 3.90, 4.00,
    4.10, 4.20, 4.30, 4.40, 4.50,
    4.60, 4.70, 4.80, 4.90, 5.00,
    5.10, 5.20, 5.30, 5.40, 5.50,
    5.60, 5.70, 5.80, 5.90, 6.00,
    6.10, 6.15, 6.20, 6.30, 6.40, 6.50,
    6.60, 6.70, 6.80, 6.90, 7.00,
    7.10, 7.20, 7.30, 7.40, 7.50, 7.58,
    7.60, 7.70, 7.80, 7.90, 8.00,
    8.10, 8.20, 8.30, 8.40, 8.50,
    8.60, 8.70, 8.80, 8.90, 9.00,
    9.10, 9.20, 9.30, 9.40, 9.50,
    9.60, 9.70, 9.80, 9.90, 1.00e+1,
    1.01e+1, 1.02e+1, 1.03e+1, 1.04e+1, 1.05e+1,
    1.06e+1, 1.07e+1, 1.08e+1, 1.09e+1, 1.10e+1,
    1.11e+1, 1.12e+1, 1.13e+1, 1.14e+1, 1.15e+1,
    1.16e+1, 1.17e+1, 1.18e+1, 1.19e+1, 1.20e+1,
    1.21e+1, 1.22e+1, 1.23e+1, 1.24e+1, 1.25e+1,
    1.26e+1, 1.27e+1, 1.28e+1, 1.29e+1, 1.30e+1,
    1.31e+1, 1.32e+1, 1.33e+1, 1.34e+1, 1.35e+1,
    1.36e+1, 1.37e+1, 1.38e+1, 1.39e+1, 1.40e+1,
    1.41e+1, 1.42e+1, 1.43e+1, 1.44e+1, 1.45e+1,
    1.46e+1, 1.47e+1, 1.48e+1, 1.49e+1, 1.50e+1,
    1.51e+1, 1.52e+1, 1.53e+1, 1.54e+1, 1.55e+1,
    1.56e+1, 1.57e+1, 1.58e+1, 1.59e+1, 1.60e+1,
    1.61e+1, 1.62e+1, 1.63e+1, 1.64e+1, 1.65e+1,
    1.66e+1, 1.67e+1, 1.68e+1, 1.69e+1, 1.70e+1,
    1.71e+1, 1.72e+1, 1.73e+1, 1.74e+1, 1.75e+1,
    1.76e+1, 1.77e+1, 1.78e+1, 1.79e+1, 1.80e+1,
    1.81e+1, 1.82e+1, 1.83e+1, 1.84e+1, 1.85e+1,
    1.86e+1, 1.87e+1, 1.88e+1, 1.89e+1, 1.90e+1,
    1.91e+1, 1.92e+1, 1.93e+1, 1.94e+1, 1.95e+1,
    1.96e+1, 1.97e+1, 1.98e+1, 1.99e+1, 2.00e+1,
    2.50e+1, 3.00e+1, 3.50e+1, 4.00e+1, 4.50e+1
]

NAMESID = [
    "NGC 104", "NGC 288", "NGC 362", "Whiting 1", "NGC 1261",
    "Pal 1", "AM1", "Eridanus", "Pal 2", "NGC 1851",
    "NGC 1904", "NGC 2298", "NGC 2419", "Ko 2", "NGC 2808",
    "E 3", "Pal 3", "NGC 3201", "Pal 4", "Ko 1",
    "NGC 4147", "NGC 4372", "Rup 106", "NGC 4590", "NGC 4833",
    "NGC 5024", "NGC 5053", "NGC 5139", "NGC 5272", "NGC 5286",
    "AM 4", "NGC 5466", "NGC 5634", "NGC 5694", "IC 4499",
    "NGC 5824", "Pal 5", "NGC 5897", "NGC 5904", "NGC 5927",
    "BH 176", "NGC 5986", "Lynga 7", "Pal 14", "NGC 6093",
    "NGC 6121", "NGC 6101", "NGC 6144", "NGC 6139", "Terzan 3",
    "NGC 6171", "1636-283", "NGC 6205", "NGC 6229", "NGC 6218",
    "FSR 1735", "NGC 6235", "NGC 6254", "Pal 15", "NGC 6266",
    "NGC 6273", "NGC 6287", "NGC 6304", "NGC 6316", "NGC 6341",
    "NGC 6333", "NGC 6356", "NGC 6352", "IC 1257", "NGC 6366",
    "Terzan 4", "NGC 6362", "Liller 1", "NGC 6380", "Ton 2",
    "NGC 6388", "NGC 6402", "NGC 6401", "Pal 6", "NGC 6426",
    "Djorg 1", "Terzan 5", "NGC 6440", "NGC 6441", "UKS 1",
    "NGC 6496", "Djorg 2", "NGC 6517", "Terzan 10", "NGC 6535",
    "NGC 6528", "NGC 6539", "NGC 6540", "NGC 6544", "NGC 6541",
    "2MS-GC01", "ESO-SC06", "NGC 6553", "2MS-GC02", "IC 1276",
    "Terzan 12", "NGC 6569", "BH 261", "GLIMPSE02", "NGC 6584",
    "NGC 6626", "NGC 6638", "NGC 6637", "NGC 6642", "NGC 6652",
    "NGC 6656", "Pal 8", "GLIMPSE01", "NGC 6712", "NGC 6715",
    "NGC 6717", "NGC 6723", "NGC 6749", "NGC 6760", "NGC 6779",
    "Terzan 7", "Pal 10", "Arp 2", "NGC 6809", "Terzan 8",
    "Pal 11", "NGC 6838", "NGC 6864", "NGC 6934", "NGC 6981",
    "NGC 7006", "NGC 7089", "Pal 12", "Pal 13", "NGC 7492"
]


def gaus_legendre(fcn, a, b, points):
    """
    Integrates the function fcn(x) between a and b using the Gauss-Legendre method.
    
    Parameters
    ----------
    fcn : callable
        The function to integrate.
    a : float
        The lower limit of integration.
    b : float
        The upper limit of integration.
    points : int
        The number of points to use for the integration (20, 40, 80, 96).
        
    Returns
    -------
    float
        The integral of fcn(x) from a to b.
    """
    
    if points == 20:
        # Gauss-Legendre 20-point nodes and weights
        x = np.array([0.076526521133497333755, 0.227785851141645078080,
                      0.373706088715419560673, 0.510867001950827098004,
                      0.636053680726515025453, 0.746331906460150792614,
                      0.839116971822218823395, 0.912234428251325905868,
                      0.963971927277913791268, 0.993128599185094924786])
        
        w = np.array([0.152753387130725850698, 0.149172986472603746788,
                      0.142096109318382051329, 0.131688638449176626898,
                      0.118194531961518417312, 0.101930119817240435037,
                      0.083276741576704748725, 0.062672048334109063570,
                      0.040601429800386941331, 0.017614007139152118312])
    elif points == 40:
        # Gauss-Legendre 40-point nodes and weights
        x = np.array([0.038772417506050821933, 0.116084070675255208483,
                      0.192697580701371099716, 0.268152185007253681141,
                      0.341994090825758473007, 0.413779204371605001525,
                      0.483075801686178712909, 0.549467125095128202076,
                      0.612553889667980237953, 0.671956684614179548379,
                      0.727318255189927103281, 0.778305651426519387695,
                      0.824612230833311663196, 0.865959503212259503821,
                      0.902098806968874296728, 0.932812808278676533361,
                      0.957916819213791655805, 0.977259949983774262663,
                      0.990726238699457006453, 0.998237709710559200350])
    
        w = np.array([0.077505947978424811264, 0.077039818164247965588,
                      0.076110361900626242372, 0.074723169057968264200,
                      0.072886582395804059061, 0.070611647391286779695,
                      0.067912045815233903826, 0.064804013456601038075,
                      0.061306242492928939167, 0.057439769099391551367,
                      0.053227846983936824355, 0.048695807635072232061,
                      0.043870908185673271992, 0.038782167974472017640,
                      0.033460195282547847393, 0.027937006980023401098,
                      0.022245849194166957262, 0.016421058381907888713,
                      0.010498284531152813615, 0.004521277098533191258])
    elif points == 80:
        # Gauss-Legendre 80-point nodes and weights
        x = np.array([0.019511383256793997654, 0.058504437152420668629,
                      0.097408398441584599063, 0.136164022809143886559,
                      0.174712291832646812559, 0.212994502857666132572,
                      0.250952358392272120493, 0.288528054884511853109,
                      0.325664370747701914619, 0.362304753499487315619,
                      0.398393405881969227024, 0.433875370831756093062,
                      0.468696615170544477036, 0.502804111888784987594,
                      0.536145920897131932020, 0.568671268122709784725,
                      0.600330622829751743155, 0.631075773046871966248,
                      0.660859898986119801736, 0.689637644342027600771,
                      0.717365185362099880254, 0.744000297583597272317,
                      0.769502420135041373866, 0.793832717504605449949,
                      0.816954138681463470371, 0.838831473580255275617,
                      0.859431406663111096977, 0.878722567678213828704,
                      0.896675579438770683194, 0.913263102571757654165,
                      0.928459877172445795953, 0.942242761309872674752,
                      0.954590766343634905493, 0.965485089043799251452,
                      0.974909140585727793386, 0.982848572738629070418,
                      0.989291302499755531027, 0.994227540965688277892,
                      0.997649864398237688900, 0.999553822651630629880])
        
        w = np.array([0.039017813656306654811, 0.038958395962769531199,
                      0.038839651059051968932, 0.038661759774076463327,
                      0.038424993006959423185, 0.038129711314477638344,
                      0.037776364362001397490, 0.037365490238730490027,
                      0.036897714638276008839, 0.036373749905835978044,
                      0.035794393953416054603, 0.035160529044747593496,
                      0.034473120451753928794, 0.033733214984611522817,
                      0.032941939397645401383, 0.032100498673487773148,
                      0.031210174188114701642, 0.030272321759557980661,
                      0.029288369583267847693, 0.028259816057276862397,
                      0.027188227500486380674, 0.026075235767565117903,
                      0.024922535764115491105, 0.023731882865930101293,
                      0.022505090246332461926, 0.021244026115782006389,
                      0.019950610878141998929, 0.018626814208299031429,
                      0.017274652056269306359, 0.015896183583725688045,
                      0.014493508040509076117, 0.013068761592401339294,
                      0.011624114120797826916, 0.010161766041103064521,
                      0.008683945269260858426, 0.007192904768117312753,
                      0.005690922451403198649, 0.004180313124694895237,
                      0.002663533589512681669, 0.001144950003186941534])
    elif points == 96:
        # Gauss-Legendre 96-point nodes and weights
        x = np.array([0.016276744849602969579, 0.048812985136049731112,
                      0.081297495464425558994, 0.113695850110665920911,
                      0.145973714654896941989, 0.178096882367618602759,
                      0.210031310460567203603, 0.241743156163840012328,
                      0.273198812591049141487, 0.304364944354496353024,
                      0.335208522892625422616, 0.365696861472313635031,
                      0.395797649828908603285, 0.425478988407300545365,
                      0.454709422167743008636, 0.483457973920596359768,
                      0.511694177154667673586, 0.539388108324357436227,
                      0.566510418561397168404, 0.593032364777572080684,
                      0.618925840125468570386, 0.644163403784967106798,
                      0.668718310043916153953, 0.692564536642171561344,
                      0.715676812348967626225, 0.738030643744400132851,
                      0.759602341176647498703, 0.780369043867433217604,
                      0.800308744139140817229, 0.819400310737931675539,
                      0.837623511228187121494, 0.854959033434601455463,
                      0.871388505909296502874, 0.886894517402420416057,
                      0.901460635315852341319, 0.915071423120898074206,
                      0.927712456722308690965, 0.939370339752755216932,
                      0.950032717784437635756, 0.959688291448742539300,
                      0.968326828463264212174, 0.975939174585136466453,
                      0.982517263563014677447, 0.988054126329623799481,
                      0.992543900323762624572, 0.995981842987209290650,
                      0.998364375863181677724, 0.999689503883230766828])
        
        w = np.array([0.032550614492363166242, 0.032516118713868835987,
                      0.032447163714064269364, 0.032343822568575928429,
                      0.032206204794030250669, 0.032034456231992663218,
                      0.031828758894411006535, 0.031589330770727168558,
                      0.031316425596861355813, 0.031010332586313837423,
                      0.030671376123669149014, 0.030299915420827593794,
                      0.029896344136328385984, 0.029461089958167905970,
                      0.028994614150555236543, 0.028497411065085385646,
                      0.027970007616848334440, 0.027412962726029242823,
                      0.026826866725591762198, 0.026212340735672413913,
                      0.025570036005349361499, 0.024900633222483610288,
                      0.024204841792364691282, 0.023483399085926219842,
                      0.022737069658329374001, 0.021966644438744349195,
                      0.021172939892191298988, 0.020356797154333324595,
                      0.019519081140145022410, 0.018660679627411467385,
                      0.017782502316045260838, 0.016885479864245172450,
                      0.015970562902562291381, 0.015038721026994938006,
                      0.014090941772314860916, 0.013128229566961572637,
                      0.012151604671088319635, 0.011162102099838498591,
                      0.010160770535008415758, 0.009148671230783386633,
                      0.008126876925698759217, 0.007096470791153865269,
                      0.006058545504235961683, 0.005014202742927517693,
                      0.003964554338444686674, 0.002910731817934946408,
                      0.001853960788946921732, 0.000796792065552012429])

    else:
        raise ValueError("Supported point values are 20, 40, 80, and 96.")
            
    # Initialize the area
    area = 0.0
    
    # Perform the integration
    for i in range(points):
        xi = (b - a) / 2.0 * x[i] + (b + a) / 2.0
        area += w[i] * fcn(xi)
        xi = -(b - a) / 2.0 * x[i] + (b + a) / 2.0
        area += w[i] * fcn(xi)
        
    area *= (b - a) / 2.0
    
    return area
# Define INTGAU function
def intgau(fcn, area):
    global indgau
    
    if indgau == 1:
        area[0] = gaus_legendre(fcn, ESTR1, ESTR2, 20)
    elif indgau == 2:
        area[0] = gaus_legendre(fcn, ESTR1, ESTR2, 40)
    elif indgau == 3:
        area[0] = gaus_legendre(fcn, ESTR1, ESTR2, 80)
    elif indgau == 4:
        area[0] = gaus_legendre(fcn, ESTR1, ESTR2, 96)
    else:
        raise ValueError(f"Invalid indgau value: {indgau}")

# Define OUTP function
def outp(xx, y, dery, ihlf, ndim, prmt):
    global NFN, X, W, DW, RAPP, RPP

    if y[0] > 0.0:
        xrif = X[NFN] + prmt[2]
        if xx < xrif:
            return
    else:
        prmt[4] = 1.0
        return

    NFN += 1
    if NFN == 8000:
        print('TROPPI PASSI')
        prmt[4] = 1.0

    X[NFN] = xx
    W[NFN] = y[0]
    DW[NFN] = y[1]
    RAPP[NFN] = RPP

# Define FCT function
def fct(x, y, dery, fcn):
    global WWWW, ESTR2, RHO0, RPP

    WWWW = y[0]
    ESTR2 = y[0]

    if WWWW > 0.0:
        qd = np.array([0.0])  # Placeholder for qd computation
        intgau(fcn, qd)
        RPP = qd[0] / RHO0
    else:
        RPP = 0.0

    dery[0] = y[1]
    dery[1] = -(2.0 / x) * y[1] - 9.0 * RPP

# Define FCNQD function
def fcnqd(x):
    global WWWW
    return np.exp(WWWW - x) * (x**1.5)

# Define FCNQD1 function
def fcnqd1(x):
    global WWWW
    return (np.exp(WWWW - x) - 1.0) * (x**0.5)

# Define FU1 function
def fu1(x):
    global WWWW
    return (np.exp(-x) - np.exp(-WWWW)) * (x**1.5)

# Define FU2 function
def fu2(x):
    global WWWW, qq
    return (x * np.exp(-x) - WWWW * np.exp(-WWWW)) * (x**qq)

# Define FU3 function
def fu3(x):
    global WWWW
    return (np.exp(-x) - np.exp(-WWWW)) * np.log(np.exp(-x) - np.exp(-WWWW)) * (x**0.5)

# Define FU4 function
def fu4(x):
    global WWWW
    return (x * np.exp(-x) - WWWW * np.exp(-WWWW)) * np.log(np.exp(-x) - np.exp(-WWWW)) * (x**0.5)

# Define DHPCG function
def dhpcg(prmt, y, dery, ndim, ihlf, fct, outp, aux, fcn):
    n = 1
    ihlf = 0
    x = prmt[0]
    h = prmt[2]
    prmt[4] = 0.0

    for i in range(ndim):
        aux[15, i] = 0.0
        aux[14, i] = dery[i]
        aux[0, i] = y[i]

    if h * (prmt[1] - x) < 0.0:
        ihlf = 13
    elif h * (prmt[1] - x) == 0.0:
        ihlf = 12
    else:
        ihlf = 0

    fct(x, y, dery, fcn)
    outp(x, y, dery, ihlf, ndim, prmt)

    if prmt[4] != 0.0:
        return

    if ihlf <= 0:
        for i in range(ndim):
            aux[7, i] = dery[i]
        isw = 1
    else:
        return

    while True:
        if isw == 1:
            x += h
            for i in range(ndim):
                aux[1, i] = y[i]

            ihlf += 1
            x -= h
            for i in range(ndim):
                aux[3, i] = aux[1, i]

            h *= 0.5
            n = 1
            isw = 2

        while True:
            if isw == 2:
                x += h
                fct(x, y, dery, fcn)
                n = 2
                for i in range(ndim):
                    aux[1, i] = y[i]
                    aux[8, i] = dery[i]
                isw = 3

            elif isw == 3:
                delt = 0.0
                for i in range(ndim):
                    delt += aux[14, i] * abs(y[i] - aux[3, i])
                delt *= 0.06666666666666667
                if delt <= prmt[3]:
                    isw = 4
                else:
                    if ihlf < 10:
                        isw = 2
                    else:
                        ihlf = 11
                        x += h
                        continue

            elif isw == 4:
                x += h
                fct(x, y, dery, fcn)
                for i in range(ndim):
                    aux[2, i] = y[i]
                    aux[9, i] = dery[i]
                n = 3
                isw = 5

            elif isw == 5:
                n = 1
                x += h
                fct(x, y, dery, fcn)
                x = prmt[0]
                for i in range(ndim):
                    aux[10, i] = dery[i]
                    y[i] = aux[0, i] + h * (0.375 * aux[7, i] + 0.7916666666666667 * aux[8, i]
                                            - 0.20833333333333333 * aux[9, i] + 0.041666666666666667 * dery[i])

                x += h
                n += 1
                fct(x, y, dery, fcn)
                outp(x, y, dery, ihlf, ndim, prmt)

                if prmt[4] != 0.0:
                    return

                if n < 4:
                    for i in range(ndim):
                        aux[n, i] = y[i]
                        aux[n + 6, i] = dery[i]
                else:
                    istop = False
                    while not istop:
                        for i in range(ndim):
                            delt = aux[8, i] + aux[9, i]
                            delt += delt + delt
                            y[i] = aux[0, i] + 0.3333333333333333 * h * (aux[7, i] + delt + aux[10, i])

                        n += 1
                        if n < 4:
                            for i in range(ndim):
                                aux[n, i] = y[i]
                                aux[n + 6, i] = dery[i]
                            istop = True

                        else:
                            for i in range(ndim):
                                delt = aux[8, i] + aux[9, i]
                                delt += delt + delt
                                y[i] = aux[0, i] + 0.375 * h * (aux[7, i] + delt + aux[10, i])

                            n += 1
                            if n < 4:
                                for i in range(ndim):
                                    aux[n, i] = y[i]
                                    aux[n + 6, i] = dery[i]
                                istop = True

                            else:
                                if istop:
                                    return

# Example main program
def main():
    global WWWW, ESTR1, ESTR2, W0, RHO0, indgau
    global X, W, DW, RAPP, NFN, RPP

    # Initialize variables and parameters
    PRMT = np.zeros(5)
    Y = np.zeros(2)
    DERY = np.zeros(2)
    AUX = np.zeros((16, 2))

    NFN = 0
    WWWW = 0.0

    # Set integration parameters
    PRMT[0] = 1.0  # Example starting value for x
    PRMT[1] = 10.  # Example end value for x
    PRMT[2] = 0.01  # Example step size
    PRMT[3] = 0.01  # Example tolerance
    PRMT[4] = 0.0  # Flag for termination

    indgau = 2  # Example setting for number of Gauss-Legendre points

    for i in range(6):
        W0 = 4.0 + i  # Example values for W0
        RHO0 = 1.0

        # Call DHPCG subroutine to integrate the King model
        dhpcg(PRMT, Y, DERY, 2, 0, fct, outp, AUX, fcnqd)

    # Example output
    print("X:", X[:NFN+1])
    print("W:", W[:NFN+1])
    print("DW:", DW[:NFN+1])
    print("RAPP:", RAPP[:NFN+1])

if __name__ == "__main__":
    main()
