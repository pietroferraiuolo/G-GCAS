import os as _os
import shutil as _sh
import subprocess as _subp
from astropy.table import QTable as _QTable
from ggcas._utility import (
    get_file_list as _get_file_list,
    MCLUSTER_SOURCE_CODE as _MCSC,
    SIMULATION_FOLDER as _SF,
    _timestamp
)

_mcluster = _os.path.join(_MCSC, "mcluster")
_mcluster_sse = _os.path.join(_MCSC, "mcluster_sse")

def mcluster_run(SSE:bool=False, **arguments):
    """
    Run mcluster with the specified arguments.

    Parameters
    ----------
    SSE : bool
        If True, the SSE version of mcluster will be used. Default is False.
        This version of mcluster includes routines for stellar evolution.
    arguments : dict
        A dictionary containing the arguments to be passed to mcluster. The keys
        are the argument names and the values are the argument values.

        N <number> (number of stars)                                 
        M <value> (mass of cluster; specify either N or M)           
        P <0|1|2|3|-1> (density profile; 0= Plummer, 1= King (1966), 
                            2= Subr et al. (2007) mass-segregated,            
                            3= 2-dimensional EFF/Nuker template,              
                            -1= no density gradient)                          
        W <0.2-inf> (W0 parameter for King model)                       
        R <value> (half-mass radius [pc], ignored for P = 3;         
                    if -1, use Marks & Kroupa (2012) Mcl-Rh relation) 
        r <value> (scale radius of EFF/Nuker template [pc])          
        c <value> (cut-off radius of EFF/Nuker template [pc])        
        g <value> (power-law slope(s) of EFF/Nuker template; use     
                    once for EFF template; use three times for Nuker  
                    template (outer slope, inner slope, transition)   
        S <0.0-1.0> (degree of mass segregation; 0.0= no segregation)
        D <1.6-3.0> (fractal dimension; 3.0= no fractality)          
        T <value> (tcrit in N-body units,                            
                    in Myr if stellar evolution is on)                
        Q <value> (virial ratio)                                     
        C <0|1|3|5> (code; 0= Nbody6, 1= Nbody4, 3= table of stars, 5= Nbody6++)    
        A <value> (dtadj in N-body units)                            
        O <value> (deltat in N-body units)                           
        G <0|1> (GPU usage; 0= no GPU, 1= use GPU)                   
        o <name> (output name of cluster model)                      
        f <0|1|2|3|4> (IMF; 0= no IMF, 1= Kroupa (2001),             
                        2= user defined, 3= Kroupa (2001) with optimal sampling,
                        4= L3 IMF (Maschberger 2012))                           
        a <value> (IMF slope; for user defined IMF, may be used      
                    multiple times, from low mass to high mass;       
                    for L3 IMF use three times for alpha, beta and mu)
        m <value> (IMF mass limits, has to be used multiple times    
                    (at least twice), from low mass to high mass [Msun])
        B <number> (number of binary systems)                        
        b <value> (binary fraction, specify either B or b)           
        p <0|1|2|3> (binary pairing, 0= random, 1= ordered for M>5.0 Msun,
                        2= random but separate pairing for M>5.0 Msun)
                        3= random but use period distribution from Sana et al., (2012);
                        Oh, S., Kroupa, P., & Pflamm-Altenburg, J. (2015)
                        for M>5.0 Msun)
        s <number> (seed for randomization; 0= randomize by timer)   
        t <0|1|2|3> (tidal field; 0= no tidal field, 1= near-field,  
                    2= point-mass, 3= Milky-Way potential)           
        e <value> (epoch for stellar evolution [Myr])                
        Z <value> (metallicity [0.0001-0.03, 0.02 = solar])          
        X <value> (galactocentric radius vector, use 3x, [pc])       
        V <value> (cluster velocity vector, use 3x, [km/s])          
        x <value> (specify external (gas) Plummer potential, use 4x, 
                    1. gas mass [Msun], 2. Plummer radius [pc]         
                    3. decay time for gas expulsion [Myr], 4. delay    
                    time for start of gas expulsion [Myr])             
        u <0|1> (output units; 0= Nbody, 1= astrophysical)
    """
    args = []
    for key, value in arguments.items():
        args.append(f"-{key} {str(value)}")
    if SSE:
        if _mcluster_sse not in _os.listdir(_MCSC):
            if _mcluster in _os.listdir(_MCSC):
                _os.system(f"cd {_MCSC} && make clean")
            print("mcluster_see executable not found. Compiling code...")
            _os.system(f"cd {_MCSC} && make mcluster_sse")
        process = _subp.Popen([_mcluster_sse] + args, stdout=_subp.PIPE, stderr=_subp.PIPE, text=True)
    else:
        if _mcluster not in _os.listdir(_MCSC):
            if _mcluster_sse in _os.listdir(_MCSC):
                _os.system(f"cd {_MCSC} && make clean")
            print("mcluster executable not found. Compiling code...")
            _os.system(f"cd {_MCSC} && make mcluster")
        process = _subp.Popen([_mcluster] + args, stdout=_subp.PIPE, stderr=_subp.PIPE, text=True)
    # Print the program output in real-time
    for line in process.stdout:
        print(line, end="")
    # Wait for the process to finish
    process.wait()
    # Print the error message if the process failed
    if process.returncode != 0:
        for line in process.stderr:
            print(line, end="")
    else:
        data_path = _manage_output_files()
        return _QTable.read(data_path, format='ascii')

def _manage_output_files():
    """
    Move the output files to a new folder with a timestamp as the name.
    """
    tn = _timestamp()
    fold = _subp.run(['pwd'], capture_output=True, text=True)
    fold = fold.stdout.strip()
    _os.mkdir(_os.path.join(_SF, tn))
    output_data = _get_file_list(fold=fold, key='.txt')
    output_info = _get_file_list(fold=fold, key='.info')
    data_path = _os.path.join(_SF, tn, output_data.split('/')[-1])
    info_path = _os.path.join(_SF, tn, output_info.split('/')[-1])
    _sh.move(output_data, data_path)
    _sh.move(output_info, info_path)
    print(data_path+"\n"+info_path)
    return data_path