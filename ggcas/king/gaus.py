"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import numpy as np

def gauss_legendre(fcn, a, b, points):
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