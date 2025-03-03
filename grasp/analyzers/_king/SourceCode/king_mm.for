c	Questo è un esempio di programma completo (King isotropo)
c	Leggilo bene e poi fallo girare e confronta i risultati in letteratura.
c	Infine trasforma il programma per la funzione di Boltzmann anisotropa
c	e confronta i risultati con l'articolo

C	CALCOLA COL FORMALISMO NEWTONIANO LE CONFIGURAZIONI DI EQUILIBRIO 
C	GRAVITAZIONALE DI UN GAS. LA FUNZIONE DI DISTRIBUZIONE DI KING E' 
C	TRATTATA COME NELL'ARTICOLO DEL 1966
C
C
c     !!!!! Questo programma è modificato per riprodurre i calcoli effettuati
c      da Da costa-Freeman 1976 (King multimassa con classi di massa)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 temp1,temp2
c      INTEGER, allocatable:: Nmasse
      COMMON/PARFCN/WWWW
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/VALINI/W0,RHO0
      COMMON/GAUSS/INDGAU
      COMMON/COSTANTI/PAI
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
     
      COMMON/VEL/v2_i(30,8000),v2(8000)
      DIMENSION Y(2),DERY(2),PRMT(5),AUX(16,2)
      DIMENSION PARW0(10000)
c      DIMENSION ymasse(10)
c      DIMENSION alfa(10)
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      EXTERNAL FCT,OUTP
      EXTERNAL FCNQD, FCNQD2
      EXTERNAL RHO,SURFDENS
      COMMON/SURF/t(8000),phi(8000),sigma(8000),psi(8000),S_S0(8000)
      COMMON/SURFI/sigma_i(30,8000),si_s0(30,8000),xni(30,8000),
     1     xn(8000)
      COMMON/TEST/C,D,P,phi0(30)
      DIMENSION SND_(80000)
      DIMENSION beta_in(10)
      DIMENSION ymu_in(10),Adj(30)
      DIMENSION yMagV(10),percmi(10),yNi(30),yMi(10)
      DIMENSION SurfBr(10,8000),SB_tot(8000),yLi(10),SB_tot2(8000)
c     COMMON/COST/yMtot
      DIMENSION Ahat(30),yK(30),Bi(30),SD_(8000)
      COMMON/SEGR/xMidr(30,8000),xNidr(30,8000)
      COMMON/SEGR2/xMiMtdr(30,8000),xNiNtdr(30,8000)
      COMMON/SEGR3/rhoidr(30,8000),xn_idr(30,8000)
      COMMON/SEGR4/rhoirhotdr(30,8000),xn_in_tdr(30,8000)
      DIMENSION phi_g(8000)
      COMMON/CV/Cv0(30,8000),Cv(30,8000),Cvqt(30,8000),Ctot0,Ctot1,Ctot2
     1     ,Cv0r(8000),Cvr(8000),Cv2r(8000)
      COMMON/MF/xMhat,xnhat_i(30,8000),xnhat(8000),xNiT(30)
      COMMON/PARFU2/qq
      COMMON/ENERGY/xk(8000),eg(8000),ephi(8000),xkt,egt,ephit,Etot,Vir2
      COMMON/ENERGY2/xki(30,8000),egi(30,8000),ephii(30,8000),
     1 xkti(30),egti(30),ephiti(30),Etoti(30),Viri(30),e_tot(8000),Etot2
      COMMON/ENERGY3/egi_2(30,8000),eg_2(8000),Vir3,egt_2,Etot3
      COMMON/IMF/xMmax,xMmin,alpha,xNtimf,xNi_imf(30)
      COMMON/CC/s2_v02,Et_Mv02,xmu,xMcap,xNhatT

      EXTERNAL fcnqd1,fu1,fu2,fu3,fu4,fdPdri
c      COMMON/NUMBER/xNiT(10)
      DIMENSION slopes(32), c_Harris(32), NGC(32)
      DIMENSION ymasse(10),alfa(10)
      DIMENSION xM_L_r(8000)
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/KINP/dPdr_i(30,8000),dPdr_i2(30,8000),dPdr1(8000),
     1     dPdr2(8000),Pr(8000),Pr2(8000)
      COMMON/CV2/x0Cv(10),xi0Cv(10),potch(8000)
      COMMON/ENERGY4/xkt_r(8000),egt_r(8000),ephit_r(8000),etot_r(8000)
c	per decidere in anticipo quali valori iniziali si vuole considerare
c	ci possono essere più parametri iniziali e quindi più DATA

c        DATA PARW0/1.00d-1,2.00d-1,3.00d-1,4.00d-1,5.00d-1,
c     1             6.00d-1,7.00d-1,8.00d-1,9.00d-1,1.00d+0,
c     1             1.10d+0,1.20d+0,1.30d+0,1.35d+0,1.40d+0,1.50d+0,
c     1             1.60d+0,1.70d+0,1.80d+0,1.90d+0,2.00d+0,
c     1             2.10d+0,2.20d+0,2.30d+0,2.40d+0,2.50d+0,
c     1             2.60d+0,2.70d+0,2.80d+0,2.90d+0,3.00d+0,
c     1             3.10d+0,3.20d+0,3.30d+0,3.40d+0,3.50d+0,
c     1             3.60d+0,3.70d+0,3.80d+0,3.90d+0,4.00d+0,
c     1             4.10d+0,4.20d+0,4.30d+0,4.40d+0,4.50d+0,
c     1             4.60d+0,4.70d+0,4.80d+0,4.90d+0,5.00d+0,
c     1             5.10d+0,5.20d+0,5.30d+0,5.40d+0,5.50d+0,
c     1             5.60d+0,5.70d+0,5.80d+0,5.90d+0,6.00d+0,
c     1             6.10d+0,6.15d+0,6.20d+0,6.30d+0,6.40d+0,6.50d+0,
c     1             6.60d+0,6.70d+0,6.80d+0,6.90d+0,7.00d+0,
c     1             7.10d+0,7.20d+0,7.30d+0,7.40d+0,7.50d+0,7.58d+0,
c     1             7.60d+0,7.70d+0,7.80d+0,7.90d+0,8.00d+0,
c     1             8.10d+0,8.20d+0,8.30d+0,8.40d+0,8.50d+0,
c     1             8.60d+0,8.70d+0,8.80d+0,8.90d+0,9.00d+0,
c     1             9.10d+0,9.20d+0,9.30d+0,9.40d+0,9.50d+0,
c     1             9.60d+0,9.70d+0,9.80d+0,9.90d+0,1.00d+1,
c     1             1.01d+1,1.02d+1,1.03d+1,1.04d+1,1.05d+1,
c     1             1.06d+1,1.07d+1,1.08d+1,1.09d+1,1.10d+1,
c     1             1.11d+1,1.12d+1,1.13d+1,1.14d+1,1.15d+1,
c     1             1.16d+1,1.17d+1,1.18d+1,1.19d+1,1.20d+1,
c     1       1.21d+1,1.22d+1,1.23d+1,1.24d+1,1.25d+1,
cc     1       1.253983561730D+1,
cc     1       1.254029d+1,
c     1       1.2540407837700014d+1, 
c     1             1.26d+1,1.27d+1,1.28d+1,1.29d+1,1.30d+1,
c     1             1.31d+1,1.32d+1,1.33d+1,1.34d+1,1.35d+1,
c     1             1.36d+1,1.37d+1,1.38d+1,1.39d+1,1.40d+1,
c     1             1.41d+1,1.42d+1,1.43d+1,1.44d+1,1.45d+1,
c     1             1.46d+1,1.47d+1,1.48d+1,1.49d+1,1.50d+1,
c     1             1.51d+1,1.52d+1,1.53d+1,1.54d+1,1.55d+1,
c     1             1.56d+1,1.57d+1,1.58d+1,1.59d+1,1.60d+1,
c     1             1.61d+1,1.62d+1,1.63d+1,1.64d+1,1.65d+1,
c     1             1.66d+1,1.67d+1,1.68d+1,1.69d+1,1.70d+1,
c     1             1.71d+1,1.72d+1,1.73d+1,1.74d+1,1.75d+1,
c     1             1.76d+1,1.77d+1,1.78d+1,1.79d+1,1.80d+1,
c     1             1.81d+1,1.82d+1,1.83d+1,1.84d+1,1.85d+1,
c     1             1.86d+1,1.87d+1,1.88d+1,1.89d+1,1.90d+1,
c     1             1.91d+1,1.92d+1,1.93d+1,1.94d+1,1.95d+1,
c     1             1.96d+1,1.97d+1,1.98d+1,1.99d+1,2.00d+1,4.0d+1/
        
c     Masse in unità solari (da mettere in unità di m1)
      
c     M3
        DATA ymasse/8.30d-1,7.70d-1,7.10d-1,6.30d-1,5.50d-1,
     1       4.70d-1,3.80d-1,3.00d-1,2.30d-1,1.20d-1/
        
c        DATA ymasse/1.d0,1.d0,1.d0,1.d0,1.d0, !!con masse tutte uguali ho king
c     1       1.d0,1.d0,1.d0,1.d0,1.d0/

c        DATA ymasse/1.9d0,1.7d0,1.5d0,1.3d0,1.1d0, !m \in [1.9,0.1]Msun
c     1       0.9d0,0.7d0,0.5d0,0.3d0,0.1d0/       

c        DATA ymasse/1.87d0,1.74d0,1.61d0,1.48d0,1.35d0, !m\in[1.87,0.7]Msun
c     1       1.22d0,1.09d0,0.96d0,0.83d0, 0.7d0/       

c        DATA ymasse/1.1d0,1.056d0,1.016d0,0.967d0,0.923d0, !m\in[1.1,0.7]Msun
c     1       0.88d0,0.834d0,0.789d0,0.745d0,0.70d0/       

c        DATA ymasse/0.5d0,0.456d0,0.411d0,0.366d0,0.322d0, !m\in[0.5,0.1]Msun
c     1       0.278d0,0.233d0,0.189d0,0.144d0,0.10d0/
c
c        DATA ymasse/2.0d0,1.83d0,1.66d0,1.50d0,1.33d0, !m\in[2.0,0.5]Msun
c     1       1.17d0,1.00d0,0.83d0,0.67d0,0.50d0/

c        DATA ymasse/1.0d1,8.5d0,6.0d0,4.5d0,3.0d0,2.0d0, !m\in[10,0.1]Msun
c     1       1.7d0,1.2d0,0.90d0,0.7d0/

c        DATA ymasse/0.2d0,0.19d0,0.18d0,0.17d0,0.16d0, !m\in[0.2,0.1]Msun
c     1       0.15d0,0.14d0,0.13d0,0.12d0,0.11d0/
c        DATA ymasse/0.3d0,0.28d0,0.26d0,0.24d0,0.22d0, !m\in[0.3,0.1]Msun
c     1       0.20d0,0.18d0,0.16d0,0.14d0,0.12d0/  
        
c        DATA ymasse/3.0d0,2.72d0,2.44d0,2.17d0,1.88d0, !m\in[3.0,0.5]Msun
c     1       1.61d0,1.33d0,1.06d0,0.78d0,0.5d0/

c        DATA ymasse/3.0d0,2.68d0,2.36d0,2.03d0,1.71d0, !m\in[3.0,0.1]Msun
c     1       1.39d0,1.07d0,0.74d0,0.42d0,0.1d0/       

c        DATA ymasse/2.0d0,1.86d0,1.71d0,1.57d0,1.42d0, !m\in[2.0,0.7]Msun
c     1       1.28d0,1.13d0,0.98d0,0.84d0,0.7d0/
        
c        DATA ymasse/1.1d0,1.03d0,0.967d0,0.90d0,0.83d0, !m\in[1.1,0.5]Msun
c     1       0.767d0,0.70d0,0.63d0,0.567d0,0.5d0/

c        DATA ymasse/1.1d0,1.01d0,0.92d0,0.83d0,0.74d0, !m\in[1.1,0.3]Msun
c     1       0.656d0,0.567d0,0.478d0,0.389d0,0.30d0/       

c        DATA ymasse/1.1d0,0.99d0,0.88d0,0.77d0,0.66d0, !m\in[1.1,0.1]Msun
c     1       0.55d0,0.44d0,0.33d0,0.22d0,0.1d0/     

c        DATA ymasse/1.1d0,1.0d0,0.9d0,0.8d0,0.7d0, !m\in[0.2,1.1]Msun
c     1       0.6d0,0.5d0,0.4d0,0.3d0,0.2d0/       

c        DATA ymasse/0.8d0,0.73d0,0.67d0,0.6d0,0.53d0, !m\in[0.2,0.8]M_sun 
c     1       0.47d0,0.40d0,0.33d0,0.27d0,0.2d0/    ! Ebrahimi2020
c
c      DATA ymasse/0.75d0,0.694d0,0.639d0,0.583d0,0.528d0, !m\in[0.25,0.75]
c     1     0.472d0,0.417d0,0.361d0,0.306d0,0.25d0/     

c      DATA ymasse/0.7d0,0.656d0,0.611d0,0.567d0,0.522d0,   !m\in[0.3,0.7]
c     1     0.478d0,0.433d0,0.389d0,0.344d0,0.30d0/
     
c      ymasse(1)=0.8d0
      Nmasse=10
c      ymstep=(0.8d0-0.2d0)/(Nmasse-1)
c      alfa(1)=1.000d0
    
c      DO k=2,Nmasse
c         ymasse(k)=ymasse(k-1)-ymstep!((0.8d0-0.2d0)/(Nmasse-1))
c         alfa(k)=1.000d0
c      ENDDO
c      Nmasse=10
      
c        DATA alfa/1.d0,1.d0,1.d0,1.d0,1.d0, 
c     1       1.d0,1.d0,1.d0,1.d0,1.d0/

c        DATA alfa/1.d+1,1.d-1,1.d-1,1.d-1,1.d-1, 
c     1       1.d-1,1.d-1,1.d-1,1.d-1,1.d-1/
        
c        DATA alfa/1.d0,0.d0,0.d0,0.d0,0.d0, 
c     1       0.d0,0.d0,0.d0,0.d0,0.d0/
        
c        Valori costanti alfa_i presenti nella funz di distribuzione
c
c        DATA alfa/1.d0,2.d0,2.5d0,5.d0,8.d0,1.4d1,2.1d1,4.0d1,
c     1       9.d1,1.d2/       
        
       DATA alfa/1.00d+0,2.86d+0,5.14d+0,1.528d+1,3.29d+1, !!!OKOKOKOK
     1       7.570d+1,1.3815d+2,3.069d+2,7.810d+2,1.050d+3/ !M3

c       DATA alfa/1.00d+0,2.86d+0,5.14d+0,1.528d+1,3.29d+1,  !range masse ridotto
c     1      0.d0,0.d0,0.d0,0.d0,0.d0/

c       DATA alfa/1.d0,1.25d0,1.60d0,2.30d0,3.47d0, !mass function alpha -2
c     1      5.57d0,10.57d0,21.52d0,47.82d0,337.45d0/      
        
c        DATA alfa/1.0000d0,2.85483376d0,5.12191381d0,1.52731683d+1,
c     1       3.29176522d+1,7.56233618d+1,1.38139141d+2,3.07569184d+2,
c     1       7.82051377d+2,1.05295421d+3/       

        
c        DATA alfa/1.00d+2,1.64d+2,1.0d+2,1.0d+2,7.56d+1, !!!OK2!
c     1       5.1d+1,3.35d+1,9.06d+0,8.37d+0,3.12d+0/ !mm14(riprodking)
        
c     1       3.788d+3,6.9075d+3,1.5345d+4,3.905d+4,5.250d+4/
        
c        DATA alfa/1.00d+5,1.00d+4,5.00d+3,1.00d+3,5.00d+2,
c     1       1.00d+2,5.00d+1,1.00d+1,5.00d+0,1.00d+0/
        
ccccc     Buone per far partire gli a_i
c         DATA alfa/1.00d-1,2.10d-1,3.00d-1,6.90d-1,1.30d+0, 
c     1       2.60d+0,4.10d+0,7.80d+0,2.22d+1,3.053d+1/
c        DATA alfa/1.d0,1.d0,1.d0,1.d0,1.d0,
c     1       1.d0,1.d0,1.d0,1.d0,1.d0/       
        
c        DATA alfa/1.00d+0,2.00d+0,2.5d0,5.0d0,7.9d0,1.40d+1,
c     1       2.00d+1,3.90d+1,8.63d+1,9.20d+1/       

       
c     Valori di beta_i e mu_i da riprodurre variando w0 e alfa(i)
       
c       DATA beta_in/6.900d+1,6.680d+1,6.400d+1,6.050d+1,5.650d+1,
c     1      5.180d+1,4.680d+1,4.140d+1,3.610d+1,2.670d+1/
       
c       DATA ymu_in/1.000d+0,1.660d+0,1.800d+0,2.890d+0,3.500d+0,
c     1      4.600d+0,4.410d+0,5.220d+0,6.930d+0,2.170d+0/

       DATA yMagV/1.27d+0,4.14d+0,5.10d+0,6.08d+0,7.36d+0,
     1      8.87d+0,1.04d+1,1.18d+1,13.3d+1,15.3d+1/      

c       DATA percmi/6.4d0,20.3d0,15.3d0,12.9d0,13.5d0,
c     1     10.2d0,8.4d0,5.3d0,4.8d0,2.9d0/

       DATA slopes/-0.45d0,-0.75d0,-0.79d0,-0.72d0,-0.74d0,
     1      -0.05d0,-1.22d0,-1.25d0,-1.21d0,-1.26d0,
     1      -1.02d0,-0.64d0,-1.14d0,-1.06d0,-0.81d0,
     1      -0.58d0,-0.16d0,-1.24d0,+0.02d0,-0.36d0,
     1      -0.57d0,-0.77d0,-0.58d0,-0.55d0,-0.78d0,
     1      -0.16d0,-0.43d0,-0.59d0,-0.83d0,-1.00d0,
     1      -0.72d0,-0.8d0/
       DATA c_Harris/2.07d0,0.99d0,1.76d0,1.16d0,1.86d0,
     1      1.38d0,1.29d0,1.41d0,1.72d0,0.74d0,
     1      1.89d0,1.41d0,1.04d0,0.86d0,1.73d0,
     1      1.23d0,1.68d0,0.80d0,1.55d0,1.34d0,
     1      1.38d0,1.68d0,1.09d0,1.86d0,1.47d0,
     1      1.11d0,2.50d0,1.38d0,0.93d0,2.29d0,
     1      1.59d0,2.50d0/
       DATA NGC/104,288,362,1261,1851,
     1      2298,3201,4590,5024,5053,
     1      5272,5286,5466,5897,5904,
     1      5986,6093,6101,6144,6218,
     1      6254,6341,6362,6541,6584,
     1      6723,6752,6779,6809,7078,
     1      7089,7099/


       
       OPEN(UNIT=1,FILE='nomefile1.dat',STATUS='replace') !,ACCESS='append')
       WRITE(1,310)'w0 c logc mu Mhat'
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='nomefile2.dat',STATUS='replace')!,ACCESS='append')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='multimass.dat',STATUS='replace')!,ACCESS='append')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='massfunc.dat',STATUS='replace')!,ACCESS='append')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='v2mean_i.dat',STATUS='replace') !,ACCESS='append')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='ai_conc_w0.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='sigma_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='rho_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='ai_Ahati_Ki.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Ahati_conc_w0.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Ki_conc_w0.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='mui_w0.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='segregation.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='segregation2.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='segregation3.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='segregation4.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='wi.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='phi.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='CvNtK.dat',STATUS='replace')
        WRITE(1,304)'w0 Cv0 Cv1 Cv2 c'
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Cv.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Er.dat',STATUS='replace')
        WRITE(1,304)'x xi w ukr egr ephir etr Kr Egr Ephir Etr'
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Etot.dat',STATUS='replace')
        WRITE(1,310)'w0 K Egr1 Egr2 Eeff Etot1 Etot3 Mhat',
     1  'Vir2 Vir3 conc s2_v02 Et_Mv02'
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Eri.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Etoti.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Et_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='K_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Egr_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Ephi_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Vir_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='betai.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='w0i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='Ni.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='SB_i.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='IMF.dat',STATUS='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='MF.dat',status='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='dPdr.dat',status='replace')
        CLOSE(UNIT=1,STATUS='keep')
        OPEN(UNIT=1,FILE='x0Cv.dat',STATUS='replace')
        WRITE(1,304)'w0 x1 x2 x3 xi1 xi2 xi3 conc c'
        CLOSE(UNIT=1,STATUS='keep')
                
 310    FORMAT(A37,A30)
c     PRINT*,'IW1,IW2,IWSTEP  1,200,1  W0 = PARW0(IW)'
c	READ*, IW1,IW2,IWSTEP
        IW1=1
        IW2=200!200
        IWSTEP=1
       
        temp1=0.0d0
        temp2=0.0d0
        call second(temp1)
c     w0 da 0.1 a 20
        PARW0(1)=0.1d0
        DO i=2,IW2
           PARW0(i)=PARW0(i-1)+0.1d0
        ENDDO
        xMu=ymasse(1)
        DO kk=1,Nmasse!10
           ym(kk) = ymasse(kk)/xMu!ymasse(1)
           Aalfa(kk) =(alfa(kk)/alfa(1))!1.00
           PRINT*,ym(kk), Aalfa(kk)
        END DO
c     ciclo per lo studio di un determinato cluster
        alpha=0.0d0
c      DO iii=6,32
c        alpha=slopes(iii)
c        c_in=c_Harris(iii)
c        ID=NGC(iii)
c        PRINT*,'Studio il cluster NGC',ID
c        PRINT*,'slope=',alpha,'c(Harris)=',c_in!,'w0=',w0,nfn
           
c     esponente della Initial Mass Function m^alpha
c        alpha=alpha+0.5d0           !-2.d0, -1.5d0, -1.0d0, -0.5d0, 0.0d0
        Kroupa=0.0d0       !0-> single powerlaw, !=0-> double powerlaw
      !ymasse(1)
c        PRINT*,'MF slope alpha=',alpha
        call InitialMassFunction(Kroupa,xMu)
        OPEN(UNIT=1,FILE='IMF.dat',STATUS='old',ACCESS='append')
        DO i=1,Nmasse!10
           WRITE(1,302)ym(i),xNi_imf(i),xNtimf
        ENDDO
        CLOSE(UNIT=1,STATUS='keep')
c     l'indice che precisa che subroutine per gli integrali si vuole usare
c     il valore 3 corrisponde a 80 punti
        
        indgau = 3

c     l'errore relativo si decide a priori e ci dà la precisione voluta	
        relerr = 1.d-9


        estr1  = 0.d0
        pai    = 3.141592654d0

        ex=0.0d0                !quando è 1 ho trovato gli ai e riparto
        ex1=0.0d0
 7766   CONTINUE
c        DO  IW=IW1,IW2,IWSTEP
c           w0 =parw0(iw)
          w0 = 16.0d0 
c           IF(ex==1.0d0)THEN
c              DO i=1,10
c                 Aalfa(i)=Bi(i)*dexp(-w0*ym(i)+D)
c                 Aalfa(i)=Ahat(i)*dexp(C*(ym(i)))!-1.d0))
c                 Aalfa(i)=Bi(i)*dexp(-w0*(ym(i)-1.d0))*
c     1            (1.d0-dexp(-w0))/(1.d0-dexp(-w0*ym(i)))
c     Aalfa(i)=Bi(i)*(dexp(w0)-1.d0)/(dexp(w0*ym(i))-1.d0)
c                 PRINT*, Aalfa(i), Bi(i)
c              ENDDO
c           ENDIF
           
c     5555      continue
           
c     PRINT*,w0
c     DO l=1, 10
c     wwww = w0*ym(l)
c     estr2=wwww
c     call intgau (fcnqd, qd)
c     psi0(l)= qd
c     ENDDO
          
 6666      continue
           w(1) = w0
           dw(1) = 0.d0
           wwww = w0
           estr2= w0
           
c     call intgau(fcnqd,qd)
        
c     Devo calcolare il rapporto Sum_i rho_i/Sum_i rho_i(0)
c     la funzione rho fa l'integrale di fcnqd e da il risultato in qd
c     la uso per calcolare int_0^(w0i) e sommare
        sum=0.0d+0
        rho0=0.0d+0
        
        call rho(fcnqd,sum)

        rho0=sum
        rapp(1) = sum/rho0

c     TEST per valutare l'impatto del termine exp(C)
        C=0.0d+0 
        D=0.0d0
C       INIZIALIZZAZIONE (per evitare la singolarità iniziale)
	AAA  = -1.5d0
	XMIN = (-relerr*w0/AAA)**0.5d0
	IIII = DLOG10(XMIN)
	if(IIII.lt.0) iiii = IIII-1
	XMIN = 10.d0**IIII

	X(2) = XMIN
	w(2) = w0+AAA*XMIN**2.d0
	dw(2)= 2.D0*AAA*XMIN

        wwww = w(2)
        estr2= w(2)
        
c        sum=0.0d+0
c     call intgau(fcnqd,qd)
        call rho(fcnqd,sum)
        rapp(2) = sum/rho0



C       FINE INIZIALIZZAZIONE
C
	nfn   = 2
        passo = x(nfn)/10.d0!200!10
	xmax  = x(nfn)*10.d0!2.d0!10

 1111	CONTINUE

	ABSERR  = relerr * w(nfn)
        if(abserr.lt.1.d-12) abserr=1.d-12

	PRMT(1) = X(NFN)	
	PRMT(2) = XMAX
	PRMT(3) = PASSO
	PRMT(4) = ABSERR

	y(1)     = w(nfn)
	y(2)     = dw(nfn)
	dery(1)  = 0.5d0
	dery(2)  = 0.5d0

	CALL DHPCG (PRMT,Y,DERY,2,IHLF,FCT,OUTP,AUX,FCNQD)

        if(ihlf.lt.11) goto 1122

c	print*,'bisezioni>10',nfn,x(nfn),w0!,w(nfn)
        passo = passo/10.d0!200!10
        xmax  = x(nfn)+(xmax-x(nfn))/10.d0!200!10
	goto 1111

 1122   continue

	if(nfn.eq.8000) goto 3333
	if(w(nfn).lt.0.d0)   GO TO 1234

        passo = passo*10.d0!2.d0!10
	xmax  = xmax*10.d0!2.d0!10

	goto 1111

 1234   continue
        xtest = (dw(nfn-1)*x(nfn-1)-w(nfn-1))/dw(nfn-1)
         
c	controllo precisione sulla routine di uscita per determinare il raggio

	if(dabs((xtest-x(nfn))/xtest).lt.1.d-9) goto 2222

        xmax  =  x(nfn)
        passo = (x(nfn)-x(nfn-1))/10.d0 !200!10
	nfn = nfn - 1
	goto 1111

 2222	continue   
        
c        print*,nfn

        conc  =  (xtest+x(nfn))/2.d0
c        PRINT*,xtest,x(nfn)

        DO kk=1, nfn
           csi(kk)=x(kk)/conc    
           IF(kk==nfn) csi(kk)=conc/conc
        END DO
        
	xlogc =   dlog10(conc) 
        xmu1  = -(4.d0*pai/9.d0)*(x(nfn-1)**2.d0)*dw(nfn-1)
        xmu2  = -(4.d0*pai/9.d0)*(x(nfn)**2.d0)*dw(nfn)
        xmu   =  (xmu1+xmu2)/2.d0

        C1 =-(x(nfn-1)*dw(nfn-1))
        C2= -(x(nfn)*dw(nfn))
        C = (C1+C2)/2.0d0         !-phi_R/sigma1^2
        D = C + w0                   !-phi_0/sigma1^2

c     Correzione delle a_i per ottenere la funzione di massa desiderata
cc     call massfunction(fcnqd1)
c        goto 8888 !inserire per bypass modifica a_i su IMF
c        IF((abs(xlogc-c_in).ge.0.1d0).AND.(ex1.eq.0.0d0))THEN
c           goto 3333
c        ELSEIF((ex1.eq.0.0d0).AND.(abs(xlogc-c_in).lt.0.1d0))THEN   !controllo solo su shape MF
cc           IF((abs((conc-c_in)).le.1.d-1).AND.ex1==0.0d0)THEN
           call massfunction(fcnqd1)


           corr = 1.0d-2
           DO j=1,Nmasse!10
cc              Adj(j)=corr
              Adj(j)=corr*(Aalfa(j))!/1.0d+1)
           ENDDO
           l=0
cc           PRINT*,'ADJUSTING: Ni/Ntot       Ni/Ntot(IMF)       a_i
cc     1diff  w0'
c           DO i=1,Nmasse!10
c               d=xNiT(i)-xNi_imf(i)
cc               PRINT*,xNiT(i),xNi_imf(i),Aalfa(i),d,w0
c            ENDDO
            
           DO i=2,Nmasse!10
              d=xNiT(i)-xNi_imf(i)
              IF(d.le.-1.0d-3)THEN
                 Aalfa(i) = Aalfa(i)+Adj(i)
              ELSEIF(d.ge.1.0d-3)THEN
                 Aalfa(i)=Aalfa(i)-Adj(i)
              ELSEIF(d.gt.-1.0d-3 .AND. d.le.-1.d-4)THEN
                 Aalfa(i) = Aalfa(i)+Adj(i)/1.d+1
              ELSEIF(d.lt.1.0d-3 .AND. d.ge.1.d-4)THEN
                 Aalfa(i) = Aalfa(i)-Adj(i)/1.d+1
              ELSEIF(d.gt.-1.0d-4 .AND. d.le.-1.d-5)THEN
                 Aalfa(i) = Aalfa(i)+Adj(i)/1.d+2
              ELSEIF(d.lt.1.0d-4 .AND. d.ge.1.d-5)THEN
                 Aalfa(i) = Aalfa(i)-Adj(i)/1.d+2
              ELSE
                 l=l+1
                 IF(l==(Nmasse-1))THEN
cc     PRINT*,k
                    
                    
cc     Se voglio fissare le a_i ai valori per cui riproduco la IMF a bassi w0
cc     devo porre ex1=1 (o diverso da 0)
cc                    PRINT*,'slope=',alpha,'c_in=',c_in,'c=',xlogc,
cc     1               'w0=',w0
c                    IF(abs(c_in-xlogc).le.1.0d-3)THEN
                       
cc     ex1=1.0d0
                       
c                       PRINT*,'slope=',alpha,'c=',xlogc,'deltac=',
c     1                      xlogc-c_in,'w0=',w0
                       goto 8888
c                    ELSEIF(xlogc-c_in.gt.0.05d0) THEN
c                       PRINT*,'Next Cluster'
c                       goto 33334                       
c                    ELSE
cc                       PRINT*,'prova w0=',w0
c                       goto 3333
c                    ENDIF
cc                    print*,'vai a 6666'
cc                    goto 6666
                 ENDIF
              ENDIF
           ENDDO
           goto 6666
cc        ELSEIF(abs(conc-c_in).gt.1.0d-2.AND.ex1==0.0d0)THEN
cc           goto 3333
cc           ELSEIF(ex1==1.0d0)THEN
cc           goto 3333
c        ENDIF

 8888   CONTINUE
        PRINT*,'Ni/Ntot     Ni/Ntot(IMF)     a_i, w0'
        DO i=1,Nmasse!10
           PRINT*,xNiT(i),xNi_imf(i),Aalfa(i), w0
        ENDDO
        
        
        
c     Calcolo massa relativa mu_i e concentrazione beta_i=r_t/r_c,i
c        IF(w0==1.2540407837700014d+1)THEN
                PRINT*,'         mu_i                    beta_i         
     1a_i'   
        call reltotmass(fcnqd)
           DO l=1,Nmasse!10
           beta(l)=conc*dsqrt(ym(l))!(ym(l)**(0.5d+0))
           PRINT*,ymu(l),beta(l),Aalfa(l)
        ENDDO
        GOTO 9999 !! mettere per bypass controllo beta_i e mu_i per M3 DF
        
        IF(ex==0.0d0)THEN
           corr = 1.0d-3
           DO j=1,Nmasse!10
              Adj(j)=corr*(Aalfa(j)/1.0d+1)
           ENDDO
c           call reltotmass(fcnqd) !mi calcola mu_i
           
           k=0
           DO j=1,Nmasse! 10
c     PRINT*,'Adjustment alfa_',j
              v=ymu(j)-ymu_in(j)
c     PRINT*,v
              IF(v.le.-1.0d-3)THEN
                 Aalfa(j) = Aalfa(j)+Adj(j)
              ELSEIF(v.ge.1.0d-3)THEN
                 Aalfa(j)=Aalfa(j)-Adj(j)
              ELSE
                 k=k+1
                 IF(k==Nmasse)THEN
c     PRINT*,k
                    goto 7777
                    
                 ENDIF
              ENDIF
           ENDDO
c     PRINT*,k
           goto 6666
        ENDIF
      
 7777   CONTINUE
        PRINT*,'         mu_i                    beta_i         
     1a_i'
 
        DO l=1,Nmasse!10
           beta(l)=conc*dsqrt(ym(l))!(ym(l)**(0.5d+0))
           PRINT*,ymu(l),beta(l),Aalfa(l)
        ENDDO
        
        vv1=beta(1)-beta_in(1)        
        IF(abs(vv1).le.0.5d+0.AND.abs(vv1).ge.1.0d-8.AND.ex==0.0d0)THEN
c     c              PRINT*,'Adjustment beta_',j
           IF(vv1.ge.0)THEN
              IF(abs(vv1).le.1.0d-6)THEN
                 w0=w0-1.0d-10
              ELSEIF(abs(vv1).le.1.0d-5.AND.abs(vv1).gt.1.0d-6)THEN
                 w0 = w0-1.00d-9
              ELSEIF(abs(vv1).le.1.0d-4.AND.abs(vv1).gt.1.0d-5)THEN
                 w0 = w0 -1.0d-7
              ELSEIF(abs(vv1).le.1.0d-3.AND.abs(vv1).gt.1.0d-4)THEN
                 w0=w0-5.0d-7
              ELSEIF(abs(vv1).le.1.0d-2.AND.abs(vv1).gt.1.0d-3)THEN
                 w0=w0-1.0d-6
              ELSEIF(abs(vv1).le.1.0d-1.AND.abs(vv1).gt.1.0d-2)THEN
                 w0=w0-1.0d-5
              ELSE
                 w0=w0-1.0d-4
              ENDIF
              goto 6666
           ELSEIF(vv1.le.0.0d0)THEN
              IF(abs(vv1).le.1.0d-6)THEN
                 w0=w0+1.0d-10
              ELSEIF(abs(vv1).le.1.0d-5.AND.abs(vv1).gt.1.0d-6)THEN
                 w0 = w0+1.00d-9
              ELSEIF(abs(vv1).le.1.0d-4.AND.abs(vv1).gt.1.0d-5)THEN
                 w0 = w0 +1.0d-7
              ELSEIF(abs(vv1).le.1.0d-3.AND.abs(vv1).gt.1.0d-4)THEN
                 w0=w0+5.0d-7
              ELSEIF(abs(vv1).le.1.0d-2.AND.abs(vv1).gt.1.0d-3)THEN
                 w0=w0+1.0d-6
              ELSEIF(abs(vv1).le.1.0d-1.AND.abs(vv1).gt.1.0d-2)THEN
                 w0=w0+1.0d-5
              ELSE
                 w0=w0+1.0d-4
              ENDIF
              goto 6666
           ENDIF
        ENDIF
        IF(abs(vv1).le.1.0d-8.AND.ex==0.0d0)THEN          
           w0fit=w0
           parw0(iw)=w0fit
           PRINT*,w0fit,parw0(iw)
           ex=1.0d0
            DO i=1,Nmasse!10
c               Bi(i)=Aalfa(i)*dexp(+w0*ym(i)-D) 
c               Ahat(i)=Aalfa(i)*dexp(C*(1.0d0-ym(i)))
c               Bi(i)=Aalfa(i)*dexp(w0*(ym(i)-1.d0))*
c     1          (1.d0-dexp(-w0*ym(i)))/(1.d0-dexp(-w0))
c  
c     Bi(i)=Aalfa(i)*(dexp(w0*ym(i))-1.d0)/(dexp(w0)-1.d0)
c               PRINT*,Bi(i)
            ENDDO
           goto 7766
        ELSEIF((abs(vv1).le.1.0d-8).AND.(ex==1.0d0))THEN
           w0fit=w0
           xxxx=1.0d0
           PRINT*,'prova',w0,w0fit
           GOTO 9999
        ENDIF
        
        
c 7766   CONTINUE
c     ENDIF
        
        IF(ex==0.0d0)THEN
           goto 3333
        ENDIF
 9999   CONTINUE

       

c     calcolo della funzione di massa e controllo su a_i
        call massfunction(fcnqd1)
        

        
c     subroutine per il calcolo della densità supericiale
       
c        call surfdens(fcnqd)

        
       
c     Calcolo velocità quadratica media
c        call v2mean(fcnqd,fcnqd2,fu1,fcnqd1)


     
        
c        yMtot = 3.30d+5
c        xNtot=0
        DO l=1,Nmasse!10
c           yNi(l)=yMtot*percmi(l)/(100.d0*ymasse(l)) !=(yMtot*ymu(l))/(ymasse(l)sum_l(ymu(l))
c           yMi(l)=yNi(l)*ymasse(l) !=yMtot*ymu(l)/sum_l ymu(l)
c           xNtot=xNtot+yNi(l)
           yLi(l)=(10.0d0)**(((4.83d0-yMagV(l))/2.5d0)) !Luminosità in solar units
        ENDDO
c        PRINT*,'Ntot=',xNtot !non ha senso, Mtot fissato e non dip da modello
c     sum=0.0d0
c        a1=1.364d0
c        b1=1.956d0
c        c1=1.079d0
c        a2=0.4439d0
c        b2=0.8863d0
c     c2=0.6608d0
        p1=0.06947d0!-0.00748d0
        p2=0.1632d0!30.54d0
        p3=-0.08146d0!8.32d0
        p4=-0.06745!68.07!68.16d0
        q1=-4.165d0!11.56d0
        q2=5.698d0!68.07d0

        DO i=1,nfn
c         yyy=a1*dexp(-((log10(x(i))-b1)/c1)**2.d0)+
c     1          a2*dexp(-((log10(x(i))-b2)/c2)**2.d0)
c     xM_L_r(i)=10.0d0**yyy
           xxxx=log10(x(i)*10.31d0)
           yyy=(p1*(xxxx**3.0d0)+p2*(xxxx**2.0d0)+ p3*xxxx+p4)/
     1          (xxxx**2.d0+q1*xxxx+q2)
c     xM_L_r(i)=yyy
           xM_L_r(i)=10.0d0**yyy
c     xL_M_r(i)=1.0d0/xM_L_r(i)

        ENDDO
cc        xM_L_r(nfn)=0.0d0
        DO n=1,nfn
          sum=0.0d0
           sum2=0.0d0
           sum3=0.0d0
           DO l=1,Nmasse!10
cc     SurfBr(l,n)=(yMi(l)/yMagV(l))*(yMagV(1)/yMi(1))*si_s0(l,n)
c              SurfBr(l,n)=(yLi(l)/yMi(l))*sigma_i(l,n)*ym(l) !sigma_i è numerica

              SurfBr(l,n)=(yLi(l)/ymu(l))*sigma_i(l,n)*ym(l)
              sum=sum+SurfBr(l,n)
              sum2=sum2+sigma_i(l,n)*ym(l) !SIGMA(R)=sum_i SIGMA_i(R)
              sum3=sum3+sigma_i(l,n) !N(R)=sum_i N_i(R) dens superf num
cc              PRINT*,SurfBr(l,n),sigma_i(l,n),sum
           ENDDO
          
           IF(n==1)THEN
              SurfBr_0=sum
              SD0=sum2
              SND0=sum3
           ENDIF
           SB_tot(n)=sum/SurfBr_0
           SD_(n)=sum2/SD0
           SB_tot2(n)=(1.d0/xM_L_r(n))*SD_(n)
           SND_(n)=sum3/SND0 !densità numerica superficiale
        ENDDO
c     PRINT*,xNtot      

        DO i=1,Nmasse!10
c           Ahat(i)=Aalfa(i)*dexp(C*(1.0d0-ym(i)))
           yK(i)=Ahat(i)/((ym(i))**3.0d0)
c           Bi(i)=Ahat(i)*dexp(D*(ym(i)-1.0d0))
c           Bi(i)=Aalfa(i)*dexp(w0*(ym(i)-1.d0))*
c     1          (1.d0-dexp(-w0*ym(i)))/(1.d0-dexp(-w0))
        ENDDO



        call segregation

        IF(w0==w0fit)THEN
           DO k=1,nfn
                phi_g(k)=(-C-w(k))!/ymasse(1) !=phi(r)/k\theta
c              PRINT*,phi_g(k),C,w(k)
           ENDDO
        ENDIF

c        call massfunction(fcnqd1)
c        call calorespecifico(fcnqd1,fu1,fu2,fu3,fu4,funzerr)
        
c        xKb=1.380649d-23 !J/K
c        xKb=1.380649d-16  !erg/K
c     xKb=8.61733262d-5  !eV/K
        xKb=1.0d0
        Ct0_NtK=xKb*Ctot0
        Ct1_NtK=xKb*Ctot1  !/xmu
        Ct2_NtK=xKb*Ctot2       !/xmu

        call energie(fcnqd1,fu1,fu3)
c        Shat=IW*0.01d0

        PRINT*,'Etot=',Etot

c        call CaloricCurve(fcnqd1)
c     scrive i parametri principali delle configurazioni in un unico file
c     w0, concentrazione, log10 concentrazione, mu (massa adimensionale)
        
c        call pgradient(fdPdri)
        
        OPEN(UNIT=1,FILE='nomefile1.dat',STATUS='old',ACCESS='append')
        WRITE(1,301) w0,conc,xlogc,xmu,xMhat
        CLOSE(UNIT=1, STATUS='keep') !,dispose='keep')
        
c        GOTO 999              !BYPASS STAMPE INUTILI
        
c     scrive i profili di ogni configurazione (uno di seguito all'altro)
c     IF(w0==3..OR.w0==5..OR.w0==8..OR.w0==12..OR.w0==12.4.OR.
c     1       w0==12.5.OR.w0==12.6.OR.w0==12.7)THEN
c     IF((abs(xlogc-1.28).le.0.01).OR.(abs(xlogc-1.48).le.0.01)
c     1  .OR.(abs(xlogc-1.68).le.0.01).OR.(abs(xlogc-1.88).le.0.01))THEN
       
c     IF(w0==6.15d0.OR.w0==7.58d0.OR.w0==8.50d0)THEN
c     IF(IW==1.OR.IW==25.OR.IW==50.OR.IW==75.OR.IW==100.OR.IW==125.OR.
c     1      IW==150.OR.IW==175.OR.IW==200.OR.w0==w0fit)THEN
      IF((abs(xlogc-0.5).le.1.d-2).OR.(abs(xlogc-1.0).le.1.d-2).OR.
     1       (abs(xlogc-1.5).le.1.d-2).OR.(abs(xlogc-2.0).le.1.d-2).OR.
     1       (abs(xlogc-2.5).le.1.d-2).OR.(abs(xlogc-3.00).le.1.d-2)
     1       .OR.((w0-40.0d0).lt.1.d-2).OR.((w0-4.0d0).lt.1.d-2))THEN
c     IF(IW==25.OR.IW==79.OR.IW==111.OR.IW==136.OR.IW==164)THEN
   
c        IF(w0.eq.w0fit)THEN
c     IF(w0-12.54.lt.1.d-2)THEN
c        IF(w0==1.0d0.OR.
          OPEN(UNIT=1,FILE='nomefile2.dat',STATUS='old',ACCESS='append')
      WRITE(1,305)'x xi w rpp SD SB SB2 SND v2 c',xlogc,'NGC',ID
          DO j=1,nfn
             WRITE(1,302)x(j),csi(j),w(j),rapp(j),SD_(j),SB_tot(j),
     1        SB_tot2(j),SND_(j),v2(j),xlogc
          END DO
          CLOSE(UNIT=1,STATUS='keep')
        ENDIF
      
c     scrive le masse e i pesi adimensionali
c        IF(IW==1.OR.IW==25.OR.IW==50.OR.IW==75.OR.IW==100.OR.IW==125.OR.
c     1      IW==150.OR.IW==175.OR.IW==200.OR.w0==w0fit) THEN
        OPEN(UNIT=1,FILE='multimass.dat',STATUS='old',ACCESS='append')
        WRITE(1,*)'m a_i c=',xlogc,'w0=',w0
           DO kk=1, Nmasse!10
              WRITE(1,302) ymasse(kk),Aalfa(kk),conc,xlogc,w0
           ENDDO
           CLOSE(UNIT=1, STATUS='keep')
c        ENDIF
c     scrive mu_i, m'_i,alfa_i, beta_i per ogni configurazione (w0)
        OPEN(UNIT=1,FILE='massfunc.dat',STATUS='old',ACCESS='append')
           WRITE(1,304) 'mu_i       m_i       a_i       beta    w0=',w0
           DO ll=1,Nmasse! 10
              WRITE(1,302)ymu(ll),ym(ll),Aalfa(ll),beta(ll)
           ENDDO
        CLOSE(UNIT=1,STATUS='keep')
       
c     scrive le velocità quadratiche medie di ciscuna classe di massa
c        IF(w0==3..OR.w0==5..OR.w0==8..OR.w0==12..OR.w0==12.4.OR.
c     1       w0==12.5.OR.w0==12.6.OR.w0==12.7)THEN
c        IF((conc-1.28.le.0.01).OR.(conc-1.48.le.0.01)
c     1       .OR.(conc-1.68.le.0.01).OR.(conc-1.88.le.0.01)) THEN
c        IF(w0==w0fit)THEN
          OPEN(UNIT=1,FILE='v2mean_i.dat',STATUS='old',ACCESS='append')
c     WRITE(1,304)'<v2>_i, w0=',w0
          DO j=1, nfn
          WRITE(1,306)v2_i(1,j),v2_i(2,j),v2_i(3,j),v2_i(4,j),v2_i(5,j),
     1             v2_i(6,j),v2_i(7,j),v2_i(8,j),v2_i(9,j),v2_i(10,j)       
            
       ENDDO
         
      CLOSE(UNIT=1,STATUS='keep')
c     ENDIF
      
c      IF(IW==1.OR.IW==25.OR.IW==50.OR.IW==75.OR.IW==100.OR.IW==125.OR.
c     1     IW==150.OR.IW==175.OR.IW==200.OR.w0==w0fit) THEN
      OPEN(UNIT=1,FILE='ai_conc_w0.dat',STATUS='old',ACCESS='append')
      WRITE(1,308)Aalfa(1),Aalfa(2),Aalfa(3),Aalfa(4),Aalfa(5),Aalfa(6),
     1     Aalfa(7),Aalfa(8),Aalfa(9),Aalfa(10),conc,w0
c      ENDIF
      CLOSE(UNIT=1,STATUS='keep')

c     IF(w0==w0fit)THEN
      IF((abs(w0-1.0d0).lt.1.d-5).OR.(abs(w0-3.0d0).lt.1.d-5).OR.
     1     (abs(w0-5.0d0).lt.1.d-5).OR.(abs(w0-7.0d0).lt.1.d-5).OR.
     1     (abs(w0-9.d0).lt.1.d-5).OR.(abs(w0-11.0d0).lt.1.d-5).OR.
     1     (abs(w0-12.54).lt.1.d-5).OR.(abs(w0-15.0d0).lt.1.d-5))THEN
         OPEN(UNIT=1,FILE='sigma_i.dat',STATUS='old',ACCESS='append')
         OPEN(UNIT=2,FILE='rho_i.dat',STATUS='old',ACCESS='append')
c     OPEN(UNIT=3,FILE='SB_i.dat',STATUS='old',ACCESS='append')
         WRITE(1,304)'s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 xi x w0=',w0
         WRITE(2,304)'r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 xi x w0=',w0
      DO j=1,nfn
         WRITE(1,306)si_s0(1,j),si_s0(2,j),si_s0(3,j),si_s0(4,j),
     1        si_s0(5,j),si_s0(6,j),si_s0(7,j),si_s0(8,j),si_s0(9,j),
     1        si_s0(10,j),csi(j),x(j)
         WRITE(2,306)rhoi(1,j),rhoi(2,j),rhoi(3,j),rhoi(4,j),rhoi(5,j),
     1        rhoi(6,j),rhoi(7,j),rhoi(8,j),rhoi(9,j),rhoi(10,j),csi(j),
     1        x(j)
c         WRITE(3,306)x(j),csi(j),SurfBr(1,j),SurfBr(2,j),SurfBr(3,j),
c     1    SurfBr(4,j),SurfBr(5,j),SurfBr(6,j),SurfBr(7,j),SurfBr(8,j),
c     1    SurfBr(9,j),SurfBr(10,j)
      ENDDO
      CLOSE(UNIT=1,STATUS='keep')
      CLOSE(UNIT=2,STATUS='keep')
c      CLOSE(UNIT=3,STATUS='keep')
      ENDIF
c      GOTO 999
c      IF
      OPEN(UNIT=1,FILE='ai_Ahati_Ki.dat',STATUS='old',ACCESS='append')
       WRITE(1,304)'mi_  ai_   Ahati_   Ki_    Bi_   C_  D_   w0=',w0
      DO k=1,Nmasse!10
         WRITE(1,309) ym(k),Aalfa(k),Ahat(k),yK(k),Bi(k),C,D
      ENDDO
      CLOSE(UNIT=1,STATUS='keep')

c      IF(IW==1.OR.IW==25.OR.IW==50.OR.IW==75.OR.IW==100.OR.IW==125.OR.
c     1     IW==150.OR.IW==175.OR.IW==200.OR.w0==w0fit) THEN
      OPEN(UNIT=1,FILE='Ahati_conc_w0.dat',STATUS='old',ACCESS='append')
      WRITE(1,308)Ahat(1),Ahat(2),Ahat(3),Ahat(4),Ahat(5),Ahat(6),
     1     Ahat(7),Ahat(8),Ahat(9),Ahat(10),conc,w0
c      ENDIF
      CLOSE(UNIT=1,STATUS='keep')
c      IF(IW==1.OR.IW==25.OR.IW==50.OR.IW==75.OR.IW==100.OR.IW==125.OR.
c     1     IW==150.OR.IW==175.OR.IW==200.OR.w0==w0fit) THEN
      OPEN(UNIT=1,FILE='Ki_conc_w0.dat',STATUS='old',ACCESS='append')
      OPEN(UNIT=2,FILE='mui_w0.dat',STATUS='old',ACCESS='append')
      WRITE(1,308)yK(1),yK(2),yK(3),yK(4),yK(5),yK(6),
     1     yK(7),yK(8),yK(9),yK(10),conc,w0
c     ENDIF
      WRITE(2,308)ymu(1),ymu(2),ymu(3),ymu(4),ymu(5),ymu(6),ymu(7),
     1     ymu(8),ymu(9),ymu(10),conc,w0
      CLOSE(UNIT=2,STATUS='keep')
      CLOSE(UNIT=1,STATUS='keep')

c     IF(w0.eq.w0fit)THEN
c      IF((w0.eq.0.1d0).OR.(w0.eq.1.d0).OR.(w0.eq.3.d0).OR.(w0.eq.5.d0)
c     1     .OR.(w0.eq.7.d0).OR.(w0.eq.9.d0).OR.(w0.eq.11.d0).OR.
c     2     (w0.eq.13.d0).OR.(w0.eq.15.d0))THEN
c     IF(abs(xlogc-1.63).lt.1.d-2)THEN
      IF((abs(w0-1.0d0).lt.1.d-5).OR.(abs(w0-4.0d0).lt.1.d-5).OR.
     1     (abs(w0-7.0d0).lt.1.d-5).OR.(abs(w0-10.0d0).lt.1.d-5).OR.
     1     (abs(w0-13.d0).lt.1.d-5).OR.(abs(w0-16.0d0).lt.1.d-5).OR.
     1     (abs(w0-0.1d0).lt.1.d-5).OR.(abs(w0-0.01d0).lt.1.d-5))THEN
c     1     (abs(w0-12.54).lt.1.d-5).OR.(abs(w0-15.0d0).lt.1.d-5))THEN
       OPEN(UNIT=1,FILE='segregation.dat',STATUS='old',ACCESS='append')
       OPEN(UNIT=2,FILE='segregation2.dat',STATUS='old',ACCESS='append')
       OPEN(UNIT=3,FILE='segregation3.dat',STATUS='old',ACCESS='append')
       OPEN(UNIT=4,FILE='segregation4.dat',STATUS='old',ACCESS='append')

       WRITE(3,304)'x xi w r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 c w0=',w0
       WRITE(4,304)'x xi w n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 c w0=',w0

       WRITE(2,304)'x xi w N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 c w0=',w0
      DO k=1, nfn
         WRITE(1,306)x(k),csi(k),w(k),xMiMtdr(1,k),xMiMtdr(2,k),
     1        xMiMtdr(3,k),xMiMtdr(4,k),xMiMtdr(5,k),xMiMtdr(6,k),
     2    xMiMtdr(7,k),xMiMtdr(8,k),xMiMtdr(9,k),xMiMtdr(10,k),w0
        
         WRITE(2,306)x(k),csi(k),w(k),xNiNtdr(1,k),xNiNtdr(2,k),
     1        xNiNtdr(3,k),xNiNtdr(4,k),xNiNtdr(5,k),xNiNtdr(6,k),
     2    xNiNtdr(7,k),xNiNtdr(8,k),xNiNtdr(9,k),xNiNtdr(10,k),xlogc,w0

         PRINT*,'test',rhoirhotdr(1,k)
         WRITE(3,306)x(k),csi(k),w(k),rhoirhotdr(1,k),rhoirhotdr(2,k),
     1        rhoirhotdr(3,k),rhoirhotdr(4,k),rhoirhotdr(5,k),
     2        rhoirhotdr(6,k),rhoirhotdr(7,k),rhoirhotdr(8,k),
     3        rhoirhotdr(9,k),rhoirhotdr(10,k),xlogc,w0
         
         WRITE(4,306)x(k),csi(k),w(k),xn_in_tdr(1,k),xn_in_tdr(2,k),
     1        xn_in_tdr(3,k),xn_in_tdr(4,k),xn_in_tdr(5,k),
     2        xn_in_tdr(6,k),xn_in_tdr(7,k),xn_in_tdr(8,k),
     3        xn_in_tdr(9,k),xn_in_tdr(10,k),xlogc,w0 
            
        
      ENDDO
      CLOSE(UNIT=1,STATUS='keep')
      CLOSE(UNIT=2,STATUS='keep')
      CLOSE(UNIT=3,STATUS='keep')
      CLOSE(UNIT=4,STATUS='keep')
      ENDIF

c      IF(w0==w0fit)THEN
         OPEN(UNIT=1,FILE='wi.dat',STATUS='old',ACCESS='append')
         DO k=1,nfn
            WRITE(1,307)x(k),w(k)*ym(1),w(k)*ym(2),w(k)*ym(3),w(k)*ym(4)
     1           ,w(k)*ym(5),w(k)*ym(6),w(k)*ym(7),w(k)*ym(8),
     1           w(k)*ym(9),w(k)*ym(10)
         ENDDO
         CLOSE(UNIT=1,STATUS='keep')
c      ENDIF

      IF(w0==w0fit)THEN
         OPEN(UNIT=1,FILE='phi.dat',status='old',access='append')
         DO k=1,nfn
            WRITE(1,302)x(k),csi(k),w(k),dw(k),phi_g(k)
         ENDDO
         CLOSE(UNIT=1,STATUS='keep')
      ENDIF
      
 999  CONTINUE
      
      IF((abs(w0-1.0d0).lt.1.d-5).OR.(abs(w0-3.0d0).lt.1.d-5).OR.
     1     (abs(w0-5.0d0).lt.1.d-5).OR.(abs(w0-7.0d0).lt.1.d-5).OR.
     1     (abs(w0-9.d0).lt.1.d-5).OR.(abs(w0-11.0d0).lt.1.d-5).OR.
     1     (abs(w0-12.54).lt.1.d-5).OR.(abs(w0-15.0d0).lt.1.d-5))THEN
         OPEN(UNIT=1,FILE='Cv.dat',STATUS='old',ACCESS='append')
         WRITE(1,304)'x xi w  dw Cv0 Cv Cvqt w0=',w0
         DO m=1,nfn
            WRITE(1,302)x(m),csi(m),w(m),dw(m),Cv0r(m)/xnhat(m),
     1           Cvr(m)/xnhat(m), Cv2r(m)/xnhat(m)
         ENDDO
         CLOSE(UNIT=1,STATUS='keep')

      ENDIF
      
      OPEN(UNIT=1,FILE='CvNtK.dat',status='old',access='append')
      WRITE(1,302) w0,Ct0_NtK,Ct1_NtK,Ct2_NtK,xlogc!conc
      CLOSE(UNIT=1,STATUS='keep')
      
     
c      goto 88888 !BYPASS STAMPE
      
      OPEN(UNIT=1,FILE='Er.dat',status='old',access='append')
      OPEN(UNIT=2,FILE='Etot.dat',status='old',access='append')
      OPEN(UNIT=3,FILE='Eri.dat',status='old',access='append')
      OPEN(UNIT=4,FILE='Etoti.dat',status='old',access='append')
     
c      IF((w0==0.8d0).OR.(w0==1.35d0).OR.(w0==2.0d0).OR.
c     1     (w0==3.5d0).OR.(w0==5.d0).OR.(w0==6.0d0).OR.
c     1     (w0==w0fit).OR.(w0==40.0d0))THEN
      IF((abs(w0-1.0d0).lt.1.d-5).OR.(abs(w0-3.0d0).lt.1.d-5).OR.
     1     (abs(w0-5.0d0).lt.1.d-5).OR.(abs(w0-7.0d0).lt.1.d-5).OR.
     1     (abs(w0-9.d0).lt.1.d-5).OR.(abs(w0-11.0d0).lt.1.d-5).OR.
     1     (abs(w0-6.0d0).lt.1.d-5).OR.(abs(w0-15.0d0).lt.1.d-5))THEN
         WRITE(1,304)'x xi w ukr egr ephir etr Kr Egr Ephir Etr w0',w0
         DO k=1,nfn
            WRITE(1,306)x(k),csi(k),w(k),xk(k)/e_tot(k),eg(k)/e_tot(k),
     1           ephi(k)/e_tot(k),e_tot(k),xkt_r(k),egt_r(k),
     1           ephit_r(k),etot_r(k)
         ENDDO
      ENDIF
      
      WRITE(2,306)w0,xkt,egt,egt_2,ephit,Etot,Etot3,xMcap,!hat(shat*Etot/(4.d0*pai*xMhat))
     1     Vir2,Vir3,xlogc,s2_v02,Et_Mv02
      DO i=1,Nmasse!10
         WRITE(4,306)w0,xkti(i),egti(i),ephiti(i),Etoti(i),Viri(i),conc,
     1        s2_v02
      ENDDO
      
      CLOSE(UNIT=1,STATUS='keep')
      CLOSE(UNIT=2,STATUS='keep')
      CLOSE(UNIT=3,STATUS='keep')
      CLOSE(UNIT=4,STATUS='keep')

c      GOTO 88888
      
      OPEN(UNIT=1,FILE='Et_i.dat',status='old',access='append')
      OPEN(UNIT=2,FILE='K_i.dat',status='old',access='append')
      OPEN(UNIT=3,FILE='Egr_i.dat',status='old',access='append')
      OPEN(UNIT=4,FILE='Ephi_i.dat',status='old',access='append')
      OPEN(UNIT=5,FILE='Vir_i.dat',status='old',access='append')
      
      WRITE(1,306)w0,Etoti(1),Etoti(2),Etoti(3),Etoti(4),Etoti(5),
     1     Etoti(6),Etoti(7),Etoti(8),Etoti(9),Etoti(10),s2_v02
      WRITE(2,306)w0,xkti(1),xkti(2),xkti(3),xkti(4),xkti(5),
     1     xkti(6),xkti(7),xkti(8),xkti(9),xkti(10),s2_v02
      WRITE(3,306)w0,egti(1),egti(2),egti(3),egti(4),egti(5),
     1     egti(6),egti(7),egti(8),egti(9),egti(10),s2_v02
      WRITE(4,306)w0,ephiti(1),ephiti(2),ephiti(3),ephiti(4),
     1     ephiti(5),ephiti(6),ephiti(7),ephiti(8),ephiti(9),ephiti(10),
     1     s2_v02
      WRITE(5,306)w0,Viri(1),Viri(2),Viri(3),Viri(4),Viri(5),Viri(6),
     1 Viri(7),Viri(8),Viri(9),Viri(10),s2_v02
     
      
      
      CLOSE(UNIT=1,STATUS='keep')
      CLOSE(UNIT=2,STATUS='keep')
      CLOSE(UNIT=3,STATUS='keep')
      CLOSE(UNIT=4,STATUS='keep')
      CLOSE(UNIT=5,STATUS='keep')
            

      OPEN(UNIT=1,FILE='betai.dat',status='old',access='append')
      WRITE(1,307)beta(1),beta(2),beta(3),beta(4),beta(5),beta(6),
     1     beta(7),beta(8),beta(9),beta(10)
      CLOSE(UNIT=1,status='keep')
       OPEN(UNIT=1,FILE='w0i.dat',status='old',access='append')
      WRITE(1,307)w0*ym(1),w0*ym(2),w0*ym(3),w0*ym(4),w0*ym(5),w0*ym(6),
     1      w0*ym(7),w0*ym(8),w0*ym(9),w0*ym(10)
      CLOSE(UNIT=1,status='keep')


      OPEN(UNIT=1,FILE='Ni.dat',status='old',access='append')
      WRITE(1,307)w0,xNiT(1),xNiT(2),xNiT(3),xNiT(4),xNiT(5),xNiT(6),
     1 xNiT(7),xNiT(8),xNiT(9),xNiT(10) 
      CLOSE(UNIT=1,status='keep')

c      IF(w0==0.1d0.OR.w0==1.0d0.or.w0==4.0d0.or.w0==7.0d0.or.
c     1     w0==10.0d0.or.w0==13.0d0.or.w0==16.0d0.or.w0==19.0d0)THEN     
      OPEN(UNIT=1,FILE='MF.dat',status='old',access='append')
      DO i=1,Nmasse!10
         WRITE(1,307)w0, xNiT(i),ym(i),ymasse(i)
      ENDDO
      CLOSE(UNIT=1,status='keep')
c      ENDIF

      IF(abs(w0-12.54d0).lt.1.d-5)THEN
         OPEN(UNIT=1,FILE='dPdr.dat',status='old',access='append')
         DO m=1,nfn-1
            WRITE(1,307)x(m),csi(m),dPdr1(m),dPdr2(m),Pr(m),Pr2(m),w0
         ENDDO
      ENDIF

88888 CONTINUE
      
      OPEN(UNIT=1,FILE='x0Cv.dat',status='old',access='append')
      WRITE(1,302) w0,x0Cv(1),x0Cv(2),x0Cv(3),xi0Cv(1),xi0Cv(2),xi0Cv(3)
     1     ,conc,xlogc
      CLOSE(UNIT=1,STATUS='keep')

      
 301  FORMAT(1pd9.2,1pd10.3,1x,1p2d15.8,1p2d15.8)
 302  FORMAT(1p10d16.8)
 303  FORMAT(1p4d16.8)          !,1pd9.8,1pd12.8)
 304  FORMAT(A50,1pd16.8)
 305  FORMAT(A30,1pd9.2,A10,I5)
 306  FORMAT(1p15d16.8)
 307  FORMAT(1p11d16.8)
 308  FORMAT(1p12d16.8)
 309  FORMAT(1p7d15.8)
      
       
 3333 CONTINUE
c      ENDDO      ! enddo ciclo w0
      
33334 CONTINUE
c      ENDDO      !enddo ciclo su 32 GCs
      
      CALL SECOND(TEMP2)
      T_imp=0.0d0
      T_imp=(temp2-temp1)
      PRINT*,'Tempo programma(s)=',T_imp,'T(min)=',T_imp/60.0
      DO j=1,2
         print *, char(7)
         call sleep(1)
      ENDDO
      
      
      STOP
      END

      SUBROUTINE pgradient(fdPdri)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/VALINI/W0,RHO0
      COMMON/GAUSS/INDGAU
      COMMON/COSTANTI/PAI
      COMMON/PARFCN/WWWW
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
c      COMMON/SURF/t(8000),phi(8000),sigma(8000),psi(8000),S_S0(8000)
      COMMON/SURFI/sigma_i(30,8000),si_s0(30,8000),xni(30,8000),
     1     xn(8000)
      COMMON/MF/xMhat,xnhat_i(30,8000),xnhat(8000),xNiT(30)
      COMMON/PARFU2/qq
      COMMON/KINP/dPdr_i(30,8000),dPdr_i2(30,8000),dPdr1(8000),
     1     dPdr2(8000),Pr(8000),Pr2(8000)
c      PRINT*,'prova'
      DO m=1,nfn-1
         sum=0.0d0
         sum2=0.0d0
         DO i=1,Nmasse
            wwww=w(m)*ym(i)
            estr2=wwww
            qd=0.0d0
            call intgau(fdPdri,qd5)
            dPdr_i(i,m)=Aalfa(i)*(1.d0/dsqrt(ym(i)))*qd5*(2.d0/3.d0)
            dPdr_i2(i,m)=Aalfa(i)*(1.d0/dsqrt(ym(i)))*xnhat_i(i,m)
            sum=sum+dPdr_i(i,m)
            sum2=sum2+dPdr_i2(i,m)
         ENDDO
         dPdr1(m)=dw(m)*sum
         dPdr2(m)=dw(m)*sum
         differr=abs(dPdr1(m)-dPdr2(m))/dPdr1(m)      
c     PRINT*,dPdr1(m),dPdr2(m),differr
         sum3=0.0d0
         sum4=0.0d0
         h=0.0d0
         b=0.0d0
         b2=0.0d0
         DO l=m,nfn-1
            b=dPdr1(m)+dPdr1(m+1)
            b2=dPdr2(m)+dPdr2(m+1)
            h=x(m+1)-x(m)
            sum3=sum3+b*h*0.5d0
            sum4=sum4+b2*h*0.5d0

         ENDDO
         IF(m.eq.1)THEN
            Pr0=-sum3
            Pr20=-sum4
         ENDIF
         Pr(m)=-sum3/Pr0
         Pr2(m)=-sum4/Pr20
c         PRINT*,Pr(m),Pr2(m),(Pr(m)-Pr2(m))/Pr(m)
      ENDDO
      
     
         
         

      
      RETURN
      END
      
      
      
c     subroutine per il calcolo della curva calorica
      SUBROUTINE CaloricCurve(fcnqd1)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/VALINI/W0,RHO0
      COMMON/GAUSS/INDGAU
      COMMON/COSTANTI/PAI
      COMMON/PARFCN/WWWW
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/ENERGY/xk(8000),eg(8000),ephi(8000),xkt,egt,ephit,Etot,Vir2
      COMMON/ENERGY2/xki(30,8000),egi(30,8000),ephii(30,8000),
     1 xkti(30),egti(30),ephiti(30),Etoti(30),Viri(30),e_tot(8000),Etot2
      COMMON/ENERGY3/egi_2(30,8000),eg_2(8000),Vir3,egt_2,Etot3
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      COMMON/MF/xMhat,xnhat_i(30,8000),xnhat(8000),xNiT(30)
      COMMON/CC/s2_v02,Et_Mv02,xmu,xMcap,xNhatT
      DIMENSION rhat(8000)
       sum=0.0d0
      DO i=1,Nmasse!10
         z=0.0d0
         wwww=w0*ym(i)
         estr2=wwww
         z=w0*ym(i)
         
         qd=0.0d0
         call intgau(fcnqd1,qd) !qd=nhat_i(w0_i)
         sum=sum+Aalfa(i)*(1.d0/dsqrt(ym(i)))*qd !xnhat_i(i,m)
c     sto sommando le rho0hat_i
      ENDDO
      
      rho0hat=sum
      xMcap=((9.d0/(4.d0*PAI))**(1.5d0))*(1.d0/dsqrt(rho0hat))*xmu

      DO i=1,nfn
         rhat(i)=dsqrt(9.0d0/(4.d0*PAI))*(1.d0/dsqrt(rho0hat))*x(i)
      ENDDO
       sum1=0.0d0
      sum2=0.0d0
      sum3=0.0d0
      sum4=0.0d0
      b1=0.0d0
      b2=0.0d0
      b3=0.0d0
      b4=0.0d0
      h=0.0d0
      sum11=0.0d0
      sum2_2=0.d0
      b2_2=0.0d0
      b11=0.0d0
      DO k=1,nfn-1
         b1=xk(k+1)*rhat(k+1)*rhat(k+1)+xk(k)*rhat(k)*rhat(k)
         b2=eg(k+1)*rhat(k+1)*rhat(k+1)+eg(k)*rhat(k)*rhat(k)
         b3=ephi(k+1)*rhat(k+1)*rhat(k+1)+ephi(k)*rhat(k)*rhat(k)
         b4=xnhat(k+1)*rhat(k+1)*rhat(k+1)+xnhat(k)*rhat(k)*rhat(k)
         h=rhat(k+1)-rhat(k)
         sum1=sum1+b1*h*0.5d0
         sum2=sum2+b2*h*0.5d0
         sum3=sum3+b3*h*0.5d0
         sum4=sum4+b4*h*0.5d0
        
c         b2_2=eg_2(k+1)*csi(k+1)*csi(k+1)+eg_2(k)*csi(k)*csi(k)
c         b11=e_tot(k+1)*csi(k+1)*csi(k+1)+e_tot(k)*csi(k)*csi(k)
c         sum2_2=sum2_2+b2_2*h*0.5d0       
c         sum11=sum11+b11*h*0.5d0
c     PRINT*,sum2
      ENDDO
      xkt=0.0d0
      egt=0.0d0
c      egt_2=0.d0
      ephit=0.0d0
      xkt=sum1*4.d0*pai
      egt=sum2*4.d0*pai
c      egt_2=sum2_2*4.d0*pai
      ephit=sum3*4.d0*pai
      Etot=xkt+egt+ephit
      
      xNhatT=(4.0d0*dsqrt(2.0d0)*pai*sum4)*4.0d0*pai
c      PRINT*,xNhatT
      sum5=0.0d0
      DO i=1, Nmasse
         sum5=sum5+ym(i)*xNiT(i)
      ENDDO
      
      xmhatmean=sum5

      
      espo=4.d0/3.d0
c      s2_v02=(1.d0/(xMcap**espo))
      s2_v02=(1.d0/(xNhatT**espo))
c     Et_Mv02=s2_v02*((4.d0*pai*(conc**3.d0)*Etot)/(rho0hat*xmu))
      Et_Mv02=s2_v02*(Etot/(xMhatmean*xNhatT))
      xkt=s2_v02*(xkt/(xMhatmean*xNhatT))
      egt=s2_v02*(egt/(xMhatmean*xNhatT))
      ephit=s2_v02*(ephit/(xMhatmean*xNhatT))
c     Et_Mv02=s2_v02*(Etot/xMcap)

        DO i=1,Nmasse!10
         b4=0.0d0
         b5=0.0d0
         b6=0.0d0
         h2=0.0d0
         sum4=0.0d0
         sum5=0.0d0
         sum6=0.0d0
         DO k=1,nfn-1
            b4=xki(i,k+1)*rhat(k+1)*rhat(k+1)+xki(i,k)*rhat(k)*rhat(k)
            b5=egi(i,k+1)*rhat(k+1)*rhat(k+1)+egi(i,k)*rhat(k)*rhat(k)
            b6=ephii(i,k+1)*rhat(k+1)*rhat(k+1)+ephii(i,k)*rhat(k)*
     1           rhat(k)
            h2=rhat(k+1)-rhat(k)
            sum4=sum4+b4*h2*0.5d0
            sum5=sum5+b5*h2*0.5d0
            sum6=sum6+b6*h2*0.5d0
         ENDDO
         xkti(i)=sum4*4.d0*pai
         egti(i)=sum5*4.d0*pai
         ephiti(i)=sum6*4.d0*pai
         
c     Etoti(i)=(sum4+sum5+sum6)
         Etoti(i)=xkti(i)+egti(i)+ephiti(i)
         Viri(i)=2.0d0*xkti(i)/(-egti(i))
      ENDDO

      DO i=1,Nmasse
         xkti(i)=s2_v02*(xkti(i)/xNhatT)
         egti(i)=s2_v02*(egti(i)/xNhatT)
         ephiti(i)=s2_v02*(ephiti(i)/xNhatT)
         Etoti(i)=s2_v02*(Etoti(i)/xNhatT)
      ENDDO
      

      
      RETURN
      END
      
c     subroutine per il calcolo della InitialMassFunction
      SUBROUTINE InitialMassFunction(Kroupa,xMu)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/VALINI/W0,RHO0
      COMMON/GAUSS/INDGAU
      COMMON/COSTANTI/PAI
      COMMON/PARFCN/WWWW
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/IMF/xMmax,xMmin,alpha,xNtimf,xNi_imf(30)
      DIMENSION xm(30)
      xMt=0.5d0/xMu             !=0.5/m_1
c     Per avere una IMF continua definita per casi, devo riscalarmi le masse
c     con la massa di transizione
    
      DO i=1,Nmasse 
         xm(i)=ym(i)/xMt        !(m_i/m_1)*(m_1/0.5) = m_i/0.5
c         PRINT*,'masse adimensionali IMF',xm(i)
      ENDDO
      xMmax=xm(1)
      xMmin=xm(Nmasse)
c     Dunque la IMF cambia pendenza per xm(i)=1
      
      
      IF(Kroupa.eq.0.d0)THEN
     
      z=0.d0
      z=alpha+1.d0
      qd=0.d0
      xNtimf=0.d0
      IF(z.eq.0.0d0)THEN
         qd=log(xMmax)-log(xMmin)
c         PRINT*,'qd=',qd
         ELSE
         qd=((xMmax**z)-(xMmin**z))/z
      ENDIF
      
      xNtimf=qd
c      PRINT*,'Ntot IMF=',xNtimf
      somma=0.d0
      DO i=1,Nmasse!10
         somma=somma+(xm(i)**alpha)
      ENDDO
      Deltam=0.d0
      Deltam=xNtimf/somma !assumo un Deltam costante tc sum_i Ni/Ntot=1
c      sum=0.d0
      DO i=1,Nmasse!10
         xNi_imf(i)=(xm(i)**alpha)*Deltam/xNtimf
c         sum=sum+xNi_imf(i)
      ENDDO
c      PRINT*,sum
      ELSE
        
         z1=0.d0
         z2=0.0d0
         z1=alpha+2.0d0
         z2=alpha+1.0d0
         qd1=0.0d0
         qd2=0.0d0
         IF(z1.eq.0.0d0)THEN
            qd1=log(1.0d0)-log(xMmin)
         ELSE
            qd1=((1.0d0**(z1))-(xMmin**(z1)))/(z1)
         ENDIF
         IF(z2.eq.0.0d0)THEN
            qd2=log(xMmax)-log(1.0d0)
         ELSE
            qd2=((xMmax**z2)-(1.0d0**z2))/z2
         ENDIF
         xNtimf=qd1+qd2
         
         somma=0.d0
         DO i=1,Nmasse!10
            IF(xm(i).le.1.0d0)THEN
               somma=somma+(xm(i)**(alpha+1.0d0))
            ELSEIF(xm(i).ge.1.0d0)THEN
               somma=somma+(xm(i)**alpha)
            ENDIF
            
         ENDDO
         Deltam=0.d0
         Deltam=xNtimf/somma    !assumo un Deltam costante tc sum_i Ni/Ntot=1
c     sum=0.d0
         DO i=1,Nmasse!10
            IF(xm(i).le.1.0d0)THEN
               xNi_imf(i)=(xm(i)**(alpha+1.0d0))*Deltam/xNtimf
            ELSEIF(xm(i).ge.1.0d0)THEN
               xNi_imf(i)=(xm(i)**alpha)*Deltam/xNtimf
            ENDIF
            
c            sum=sum+xNi_imf(i)
         ENDDO
c         PRINT*,sum
      ENDIF
      
      
      RETURN
      END
      
c     subroutine per il calcolo delle energie adimensionali
      SUBROUTINE energie(fcnqd1,fu1,fu3)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/VALINI/W0,RHO0
      COMMON/GAUSS/INDGAU
      COMMON/COSTANTI/PAI
      COMMON/PARFCN/WWWW
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/PARFU2/qq
      COMMON/ENERGY/xk(8000),eg(8000),ephi(8000),xkt,egt,ephit,Etot,Vir2
      COMMON/ENERGY2/xki(30,8000),egi(30,8000),ephii(30,8000),
     1 xkti(30),egti(30),ephiti(30),Etoti(30),Viri(30),e_tot(8000),Etot2
      COMMON/ENERGY3/egi_2(30,8000),eg_2(8000),Vir3,egt_2,Etot3
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      COMMON/TEST/C,D,P,phi0(30)
      COMMON/ENERGY4/xkt_r(8000),egt_r(8000),ephit_r(8000),etot_r(8000)
c     EXTERNAL fu1,fu3
    
      z=0.0d0
      
      DO m=1,nfn-1
         somma1=0.0d0
         somma2=0.0d0
         somma3=0.0d0
         somma4=0.0d0
         u1=0.d0
         u3=0.0d0
         qd=0.d0
         DO i=1,Nmasse!10
            z=w(m)*ym(i)
            wwww=w(m)*ym(i)
            estr2=wwww
            call intgau(fu1,u1)
            wwww=w(m)*ym(i)
            estr2=wwww
            call intgau(fu3,u3)
            wwww=w(m)*ym(i)
            estr2=wwww
            call intgau(fcnqd1,qd)
           
            xki(i,m)=dexp(z)*u1*Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))
     1           *4.d0*pai*dsqrt(2.0d0)
            egi(i,m)=(2.d0/3.d0)*(dexp(z)*u1+0.4d0*(z**(2.5d0)))
     1           *ym(i)*dw(m)*x(m)*Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))
     1           *4.d0*pai*dsqrt(2.0d0)
            egi_2(i,m)=Aalfa(i)*(1.d0/dsqrt(ym(i)))*qd*4.d0*pai
     1           *dsqrt(2.0d0)  !rho(r)
c             PRINT*,qd,egi(i,m)
            ephii(i,m)=-xki(i,m)-dexp(z)*u3*Aalfa(i)*
     1       (1.d0/(ym(i)**(1.5d0)))*4.d0*pai*dsqrt(2.0d0) !-dexp(z)*(u1+u3)
c            PRINT*,u3
            somma1=somma1+xki(i,m)!*Aalfa(i)*(ym(i)**(-1.5d0))
            somma2=somma2+egi(i,m)!*Aalfa(i)*(ym(i)**(-1.5d0))
            somma3=somma3+ephii(i,m) !*Aalfa(i)*(ym(i)**(-1.5d0))
            somma4=somma4+egi_2(i,m)
         ENDDO
         xk(m)=somma1
         eg(m)=somma2
         eg_2(m)=somma4*(-C-w(m))*0.5d0
         ephi(m)=somma3
         e_tot(m)=somma1+somma2+somma3
      ENDDO
      DO i=1,Nmasse!10
         xki(i,nfn)=0.0d0
         egi(i,nfn)=0.0d0
         egi_2(i,nfn)=0.d0
         ephii(i,nfn)=0.0d0
      ENDDO
      xk(nfn)=0.0d0
      eg(nfn)=0.0d0
      eg_2(nfn)=0.0d0
      ephi(nfn)=0.0d0
      e_tot(nfn)=0.0d0

     
      sum1=0.0d0
      sum2=0.0d0
      sum3=0.0d0
      b1=0.0d0
      b2=0.0d0
      b3=0.0d0
      h=0.0d0
      sum11=0.0d0
      sum2_2=0.d0
      b2_2=0.0d0
      b11=0.0d0
      DO k=1,nfn-1
         b1=xk(k+1)*csi(k+1)*csi(k+1)+xk(k)*csi(k)*csi(k)
         b2=eg(k+1)*csi(k+1)*csi(k+1)+eg(k)*csi(k)*csi(k)
         b2_2=eg_2(k+1)*csi(k+1)*csi(k+1)+eg_2(k)*csi(k)*csi(k)
         b3=ephi(k+1)*csi(k+1)*csi(k+1)+ephi(k)*csi(k)*csi(k)
         b11=e_tot(k+1)*csi(k+1)*csi(k+1)+e_tot(k)*csi(k)*csi(k)
         h=csi(k+1)-csi(k)
         sum1=sum1+b1*h*0.5d0
         sum2=sum2+b2*h*0.5d0
         sum2_2=sum2_2+b2_2*h*0.5d0
         sum3=sum3+b3*h*0.5d0
         sum11=sum11+b11*h*0.5d0
c     PRINT*,sum2
         xkt_r(k)=sum1
         egt_r(k)=sum2
         ephit_r(k)=sum3
         etot_r(k)=sum11
      ENDDO
      xkt=0.0d0
      egt=0.0d0
      egt_2=0.d0
      ephit=0.0d0
      xkt=sum1*4.d0*pai
      egt=sum2*4.d0*pai
      egt_2=sum2_2*4.d0*pai
      ephit=sum3*4.d0*pai
      Etot=0.0d0
      Etot=(sum1+sum2+sum3)
      Etot2=sum11*4.d0*pai
      Etot3=sum1+sum2_2+sum3
     
      Vir=2.d0*xkt+egt              !dovrebbe essere nullo
c      PRINT*,'2Ktot+Egr=',Vir,2.d0*xkt,egt
      Vir2=2.d0*xkt/(-egt)          !dovrebbe essere 1
c      PRINT*,'2Ktot/-Egr=',Vir2
c     PRINT*,Etot,sum1,sum2,sum3
      Vir3=2.d0*xkt/(-egt_2)
      DO i=1,Nmasse!10
         b4=0.0d0
         b5=0.0d0
         b6=0.0d0
         h2=0.0d0
         sum4=0.0d0
         sum5=0.0d0
         sum6=0.0d0
         DO k=1,nfn-1
            b4=xki(i,k+1)*csi(k+1)*csi(k+1)+xki(i,k)*csi(k)*csi(k)
            b5=egi(i,k+1)*csi(k+1)*csi(k+1)+egi(i,k)*csi(k)*csi(k)
            b6=ephii(i,k+1)*csi(k+1)*csi(k+1)+ephii(i,k)*csi(k)*csi(k)
            h2=csi(k+1)-csi(k)
            sum4=sum4+b4*h2*0.5d0
            sum5=sum5+b5*h2*0.5d0
            sum6=sum6+b6*h2*0.5d0
         ENDDO
         xkti(i)=sum4*4.d0*pai
         egti(i)=sum5*4.d0*pai
         ephiti(i)=sum6*4.d0*pai
         
         Etoti(i)=(sum4+sum5+sum6)
         Viri(i)=2.0d0*xkti(i)/(-egti(i))
      ENDDO
      
      
      return
      end

      SUBROUTINE massfunction(fcnqd1)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/VALINI/W0,RHO0
      COMMON/GAUSS/INDGAU
      COMMON/COSTANTI/PAI
      COMMON/PARFCN/WWWW
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/MF/xMhat,xnhat_i(30,8000),xnhat(8000),xNiT(30)
c      COMMON/PARFU2/qq
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      EXTERNAL FCT,OUTP
      EXTERNAL FCNQD, FCNQD2
      EXTERNAL RHO,SURFDENS
c      COMMON/SURF/t(8000),phi(8000),sigma(8000),psi(8000),S_S0(8000)
c      COMMON/SURFI/sigma_i(10,8000),si_s0(10,8000),xni(10,8000),
c     1   xn(8000)
c      COMMON/TEST/C,D,P,phi0(10)
c      COMMON/NUMBER/xNiT(10)

      
      DO m=1,nfn-1
         sum4=0.d0
         DO i=1,Nmasse!10
            z=0.0d0
            wwww=w(m)*ym(i)
            estr2=wwww
            z=w(m)*ym(i)
            
            qd=0.0d0
            call intgau(fcnqd1,qd)
            xnhat_i(i,m)=qd
            sum4=sum4+Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))*xnhat_i(i,m)
                     
         ENDDO
         xnhat(m)=sum4
         
      ENDDO
      
      DO i=1,Nmasse!10
         xnhat_i(i,nfn)=0.0d0
      ENDDO
      
      xnhat(nfn)=0.0d0
      somma3=0.d0
      b3=0.d0
      h=0.d0
      DO k=1,nfn-1
         b3=(xnhat(k+1)*x(k+1)*x(k+1)+xnhat(k)*x(k)*x(k))
         h=x(k+1)-x(k)
         somma3=somma3+b3*h*0.5d0
      ENDDO

      xMhat=somma3

       DO i=1,Nmasse!10
         b4=0.d0
         h4=0.d0
         sum4=0.d0
         DO m=1,nfn-1
            b4=xnhat_i(i,m+1)*x(m+1)*x(m+1)+xnhat_i(i,m)*x(m)*x(m)
            h4=x(m+1)-x(m)
            sum4=sum4+b4*h4*0.5d0
         ENDDO
         xNiT(i)=(Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))*sum4)/xMhat
      ENDDO
      
      
      return
      end
      
c     Subroutine per il calcolo del calore specifico adimensionale
      SUBROUTINE calorespecifico(fcnqd1,fu1,fu2,fu3,fu4,funzerr)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/VALINI/W0,RHO0
      COMMON/GAUSS/INDGAU
      COMMON/COSTANTI/PAI
      COMMON/PARFCN/WWWW
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/CV/Cv0(30,8000),Cv(30,8000),Cvqt(30,8000),Ctot0,Ctot1,Ctot2
     1 ,Cv0r(8000),Cvr(8000),Cv2r(8000)
      COMMON/MF/xMhat,xnhat_i(30,8000),xnhat(8000),xNiT(30)
      COMMON/PARFU2/qq
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      EXTERNAL FCT,OUTP
      EXTERNAL FCNQD, FCNQD2
      EXTERNAL RHO,SURFDENS
      COMMON/SURF/t(8000),phi(8000),sigma(8000),psi(8000),S_S0(8000)
      COMMON/SURFI/sigma_i(30,8000),si_s0(30,8000),xni(30,8000),
     1   xn(8000)
      COMMON/TEST/C,D,P,phi0(30)
      COMMON/CV2/x0Cv(10),xi0Cv(10),potch(8000)

c      COMMON/NUMBER/xNiT(10)
      
     
      DO k=1,nfn
         DO i=1,Nmasse!10
            Cv0(i,k)=0.0d0
            Cv(i,k)=0.0d0
            Cvqt(i,k)=0.0d0
c            xnhat(i,k)=xni(i,k)
         ENDDO
      ENDDO
      
      
      DO m=1,nfn-1
         sum1=0.0d0
         sum2=0.0d0
         sum3=0.0d0
c         sum4=0.0d0
         DO i=1,Nmasse!10
            z=0.0d0
            zz=0.0d0
            wwww=w(m)*ym(i)
            estr2=wwww
            z=w(m)*ym(i)
            
            qd=0.0d0
            u1=0.0d0
            u21=0.0d0
            u23=0.0d0
            u3=0.0d0
            u4=0.0d0
c         zz=dsqrt(z)
c         call intgau(funzerr,zzz)
c         qd=(dsqrt(PAI)/2.0d0)*dexp(z)*(derf(zz))-
c     1        (z**(0.5d0))*(1.0d0+(2.0d0*z/3.0d0))
            call intgau(fcnqd1,qd)
      
c         u1=dexp(-z)*(1.5d0*qd-(2.0d0/5.0d0)*(z**(2.5d0)))
c         u21=dexp(-z)*(1.5d0*qd-(2.0d0/3.0d0)*(z**(2.5d0)))
c     u23=dexp(-z)*(3.75d0*qd-(z**(2.5d0))*(1.0d0+0.4d0*z))

c         u21=u1-(4.d0/15.d0)*dexp(-z)*(z**(2.5d0))
c         u23=2.5d0*u1-0.4d0*dexp(-z)*(z**(3.5d0))
c     PRINT*,qd,qd1,qd2,qd-qd1

c            xnhat_i(i,m)=qd
            wwww=w(m)*ym(i)
            estr2=wwww
            call intgau(fu1,u1)
c       
c     qq=1/2 oppure 3/2 a seconda di cosa compare in fu2
            wwww=w(m)*ym(i)
            estr2=wwww
            qq=0.5d0
            call intgau(fu2,u21)
            
            qq=1.5d0
            wwww=w(m)*ym(i)
            estr2=wwww
            call intgau(fu2,u23)
c     PRINT*,u2_3,u2_33,u2_3-u2_33

            
            wwww=w(m)*ym(i)
            estr2=wwww
            call intgau(fu3,u3)
            
            
            wwww=w(m)*ym(i)
            estr2=wwww
            call intgau(fu4,u4)
            
         
c      Cv0(m)=dexp(wwww)*(u3*(1.5d0-((2.0d0/3.0d0)*(wwww**(2.5d0))/qd))
c     1        -u4)
            Cv0(i,m)=(((u3*u21*dexp(2.0d0*z))/qd)-u4*dexp(z))
c         PRINT*,u3

      
            Cv(i,m)=(Cv0(i,m)+x(m)*dw(m)*ym(i)*(qd+(2.0d0/3.0d0)
     1      *(z**(1.5d0))*(2.5d0-1.4d0*z)+((4.0d0/9.0d0)*(z**(4.0d0))*
     1      (0.4d0*z-2.0d0)/qd)+ ((16.d0/135.d0)*(z**(6.5d0))/(qd*qd))))

c      Cvqt(m)=2*PAI*(Cv0(m)-3.0d0*qd-2.0d0*(wwww**(1.5d0))+
c     1     +4.0d0*(wwww**(3.5d0))/5.0d0
c     1     +8.0d0*(wwww**(4.0d0))/(9.0d0*qd))
            Cvqt(i,m)=(Cv0(i,m)-2.d0*dexp(z)*u23+
     1           (2.d0*dexp(2.d0*z)*u1*u21/qd))
c     PRINT*,Cv0(m),Cv(m),Cvqt(m),xnhat(m),wwww
         

c      IF(w(1).eq.0.7d0)THEN
c         PRINT*,'problema strano'
c         PRINT*,w(m),qd,u1,u2_1,u2_3,u3,u4
c     ENDIF
         sum1=sum1+Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))*Cv0(i,m)
         sum2=sum2+Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))*Cv(i,m)
         sum3=sum3+Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))*Cvqt(i,m)
c         sum4=sum4+Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))*xnhat_i(i,m)
         
      ENDDO
      Cv0r(m)=sum1
      Cvr(m)=sum2
      Cv2r(m)=sum3
c      xnhat(m)=sum4
      ENDDO
       
      qd=0.0d0
      u1=0.0d0
      u21=0.0d0
      u23=0.0d0
      u3=0.0d0
      u4=0.0d0
      DO i=1,Nmasse!10
         Cv0(i,nfn)=0.0d0
         Cv(i,nfn)=0.0d0
         Cvqt(i,nfn)=0.0d0
c         xnhat_i(i,nfn)=0.0d0
      ENDDO
      Cv0r(nfn)=0.0d0
      Cvr(nfn)=0.0d0
      Cv2r(nfn)=0.0d0
c      xnhat(nfn)=0.0d0

      somma0=0.0d0
      somma1=0.0d0
      somma2=0.0d0
c      somma3=0.0d0
      b0=0.0d0
      b1=0.0d0
      b2=0.0d0
c      b3=0.0d0
      h=0.0d0      
      DO k=1,nfn-1
         b0=(Cv0r(k+1)*x(k+1)*x(k+1)+Cv0r(k)*x(k)*x(k))
         b1=(Cvr(k+1)*x(k+1)*x(k+1)+Cvr(k)*x(k)*x(k))
         b2=(Cv2r(k+1)*x(k+1)*x(k+1)+Cv2r(k)*x(k)*x(k))
c         b3=(xnhat(k+1)*x(k+1)*x(k+1)+xnhat(k)*x(k)*x(k))
         h=x(k+1)-x(k)
         somma0=somma0+b0*h*0.5d0
         somma1=somma1+b1*h*0.5d0
         somma2=somma2+b2*h*0.5d0
c         somma3=somma3+b3*h*0.5d0
      ENDDO
      Ctot0=0.0d0
      Ctot1=0.0d0
      Ctot2=0.0d0
      Ctot0=somma0/xMhat !somma3
      Ctot1=somma1/xMhat !somma3
      Ctot2=somma2/xMhat !somma3
     
c      xMhat=somma3

c      DO i=1,10
c         b4=0.d0
c         h4=0.d0
c         sum4=0.d0
c         DO m=1,nfn-1
c            b4=xnhat_i(i,m+1)*x(m+1)*x(m+1)+xnhat_i(i,m)*x(m)*x(m)
c            h4=x(m+1)-x(m)
c            sum4=sum4+b4*h4*0.5d0
c         ENDDO
c         xNiT(i)=(Aalfa(i)*(1.d0/(ym(i)**(1.5d0)))*sum4)/xMhat
c      ENDDO

      ex11=0.0d0
      ex22=0.0d0
      ex33=0.0d0
      DO m=1,nfn-1
         IF((Cvr(m).ge.0.0d0).AND.(Cvr(m+1).le.0.0d0)
     1        .AND.(ex11.eq.0.0d0))THEN
            x0Cv(1)=(x(m)+x(m+1))*0.5d0
            xi0Cv(1)=(csi(m)+csi(m+1))*0.5d0
            ex11=1.0d0
            PRINT*,'primo zero Cv: r/rk=',x0Cv(1),'r/rt=',xi0Cv(1),
     1           'w0=', w0,m
         ELSEIF((Cvr(m).le.0.0d0).AND.(Cvr(m+1).ge.0.0d0)
     1           .AND.(ex22.eq.0.0d0))THEN
            x0Cv(2)=(x(m)+x(m+1))*0.5d0
            xi0Cv(2)=(csi(m)+csi(m+1))*0.5d0
            ex22=1.0d0
            PRINT*,'secondo zero Cv: r/rk=',x0Cv(2),'r/rt=',xi0Cv(2),
     1           'w0=', w0,m
         ELSEIF((Cvr(m).ge.0.0d0).AND.(Cvr(m+1).le.0.0d0)
     1           .AND.(ex.eq.1.0d0).AND.(ex33.eq.0.0d0))THEN
               x0Cv(3)=(x(m)+x(m+1))*0.5d0
               xi0Cv(3)=(csi(m)+csi(m+1))*0.5d0
               ex33=1.0d0
               PRINT*,'terzo zero Cv: r/rk=',x0Cv(3),'r/rt=',xi0Cv(3),
     1              'w0=', w0,m
c            GOTO 9999
         ENDIF
      ENDDO
 9999 CONTINUE
           
      
      RETURN
      END

       function fcnqd1(x)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARFCN/WWWW
      fcnqd1=(dexp(wwww-x)-1.0d0)*(x**(0.5d0))
      return
      end
      
      function fu1(x)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARFCN/WWWW
      fu1 = (dexp(-x)-dexp(-wwww))*(x**1.5d0) 
      return
      end

      function fu2(x)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARFCN/WWWW
      COMMON/PARFU2/qq
      fu2 = (x*dexp(-x)-wwww*dexp(-wwww))*(x**(qq)) 
      return
      end
      
      function fu3(x)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PARFCN/WWWW
      fu3=(dexp(-x)-dexp(-wwww))*dlog(dexp(-x)-dexp(-wwww))*(x**(0.5d0))
      return
      end

      function fu4(x)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PARFCN/WWWW
      fu4=(x*dexp(-x)-wwww*dexp(-wwww))*dlog(dexp(-x)-dexp(-wwww))*
     1     (x**(0.5d0))
      return
      end

      function fdPdri(x)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PARFCN/WWWW
      fdPdri=dexp(wwww-x)*(x**(1.5d0))
      return
      end
      

      
c     Subroutine per il calcolo della densità superficiale Sigma/Sigma0
      
      SUBROUTINE surfdens(fcnqd)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/SURF/t(8000),phi(8000),sigma(8000),psi(8000),S_S0(8000)
      COMMON/SURFI/sigma_i(30,8000),si_s0(30,8000),xni(30,8000),xn(8000)
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/PARFCN/WWWW
      COMMON/VALINI/W0,RHO0
      DIMENSION yy(8000)
      qd=0.0d+0   
      DO l=1,nfn
c         wwww=w(l)
c        estr2=wwww
c         call rho(fcnqd,qd)
c         PRINT*,qd
c         psi(l)=qd
c         IF(l==nfn)psi(l)=0.0d0
c     psi(l)=rapp(l)*rho0
         somma=0.0d0
         DO n=1,Nmasse!10              !calcolo la densità numerica per ogni classe
            xni(n,l)=rhoi(n,l)/ym(n)
            somma=somma+xni(n,l)
         ENDDO
         IF(l==1)xn0=somma
         xn(l)=somma/xn0
         psi(l)=xn(l)           !proietto la densità numerica
      END DO
c      somma=0.0d0
c      base=0.0d0
c      altezza=0.0d0
c      DO m=1,nfn-1
c         base=psi(m)+psi(m+1)
c         altezza=csi(m+1)-csi(m)
c         somma=somma+base*altezza*0.5d0
c      ENDDO
c      sigma0_2=2.0d0*somma
      
c      sum=0.0d+0
c      b=0.0d0
c      h=0.0d0
c      R=0.0d0
c      t0=0.0d0
c      t1=0.0d0
cc     t2=0.0d0
c      phi0=0.0d0
c      phi1=0.0d0
c      phi2=0.0d0
c      sigma0=0.0d0
c      DO k=1, nfn            !nfn-1
c         R=csi(k)
c         sum=0.0d0
c         sum_2=0.0d0
c         DO l=k,nfn
c            t(l)=(csi(l)*csi(l)-R*R)/(1+csi(l)*csi(l))

c            phi(l)=((1+csi(l)*csi(l))**(1.5d0))*psi(l)
cc            phi(l)=psi(l)/(((1-t(l))**(1.5d0)))
cc           phi(l)=(((1+R*R)/(1-t(l)))**(1.5d0))*psi(l)
c            IF(l==k)THEN
c               sum=0.0d0
c               phi0=phi(l)     !ovvero rho(k)(1+R^2)^(3/2)
c               t0=t(l)          !che è zero
cc            ELSEIF(l==k+1)THEN
cc c              t2=t(l)
cc               phi2=phi(l)
cc     phi2=(phi1+phi0)/2.0d0
c           ELSEIF(l==k+1)THEN
c               
c               t1=t(l)
c               t2=(t0+t1)/2.0d0 ! t1/2->phi(t2)=rho(k)*((1+R^2)/(1+t1/2))^3/2
c               phi2=phi0/((1.0d0+t1/2.0d0)**(1.5d0))
c               phi1=phi(l)
c               sum=sum+(4.0d0/15.0d0)*sqrt(t1)*((5.0d0-(t1/t2))*phi0
c     1         + phi1*(3.0d0*t1-5.0d0*t2)/(2.0d0*(t1-t2))
c     1         + phi2*(t1*t1)/(t2*(t1-t2)))
c            ELSE
c             b=phi(l)*(1.0d0/dsqrt(t(l)))+phi(l-1)*(1.0d0/dsqrt(t(l-1)))
c               h=t(l)-t(l-1)
c               phi(l)=(psi(l)*csi(l)/sqrt(csi(l)*csi(l)-R*R))
c               b_2=(psi(l)*csi(l)/sqrt(csi(l)*csi(l)-R*R))+
c     1              (psi(l-1)*csi(l-1)/sqrt(csi(l-1)*csi(l-1)-R*R))
c               h_2=csi(l)-csi(l-1)
c               sum=sum+b*h*0.5d+0               
c               sum_2=sum_2+b_2*h_2*0.5
c            ENDIF
c         END DO
c         IF(k==1)THEN
c            sigma0=sum
c            sigma0=sum_2*2.0d0
c            PRINT*,'deltasigma0',sigma0-sigma0_2
c            sigma0=sum*2.0d0
c         ENDIF
c         sigma(k)=sum*(1.0d0/(1.0d0+R*R))
c         sigma(k)=sum_2*2.0d0
c         s_s0(k)=sigma(k)/sigma0
c      END DO

      
c      somma=0.0d0
c      base=0.0d0
c      altezza=0.0d0
c      DO m=1,nfn-1
c         base=rhoi(1,m)+rhoi(1,m+1)
c         altezza=csi(m+1)-csi(m)
c         somma=somma+base*altezza*0.5d0
c      ENDDO
c      sigma0_1=2.0d0*somma
      
      sigma0_1=0.0d0
      DO m=1,Nmasse!10
         
         sum=0.0d+0
         b=0.0d0
         h=0.0d0
         xiR=0.0d0
c         t0=0.0d0
c         t1=0.0d0
c         t2=0.0d0
c         phi0=0.0d0
c         phi1=0.0d0
c         phi2=0.0d0
         DO i=1,nfn
            xiR=csi(i)
            DO k=i,nfn
               yy(k)=dsqrt((csi(k)**2.d0)-(xiR**2.d0))
               IF(k==i)THEN
                  sum=0.0d0
               ELSE
                  b=xni(m,k)+xni(m,k-1)
                  h=yy(k)-yy(k-1)
                  sum=sum+b*h*0.5d0
               ENDIF
            ENDDO
            IF(i==1.and.m==1)THEN
               sigma0_1=sum
            ENDIF
            sigma_i(m,i)=sum 
            si_s0(m,i)=sigma_i(m,i)/sigma0_1
         ENDDO
      ENDDO
      
         
c         DO k=1, nfn            !nfn-1
c            R=csi(k)
c            DO l=k,nfn
c               t(l)=(csi(l)*csi(l)-R*R)/(1+csi(l)*csi(l))
c               phi(l)=((1+csi(l)*csi(l))**(1.5d0))*xni(m,l)
cc     phi(l)=psi(l)/(((1-t(l))**(1.5d0)))
cc     phi(l)=(((1+R*R)/(1-t(l)))**(1.5d0))*psi(l)
c               IF(l==k)THEN
c                  sum=0.0d0
c                  phi0=phi(l)   !ovvero rho(k)(1+R^2)^(3/2)
c                  t0=t(l)       !che è zero             
c               ELSEIF(l==k+1)THEN              
c                  t1=t(l)
c                  t2=(t0+t1)/2.0d0 ! t1/2->phi(t2)=rho(k)*((1+R^2)/(1+t1/2))^3/c2
c                  phi2=phi0/((1.0d0+t1/2.0d0)**(1.5d0))
c                  phi1=phi(l)
c                  sum=sum+(4.0d0/15.0d0)*dsqrt(t1)*((5.0d0-(t1/t2))*phi0
c     1                 + phi1*(3.0d0*t1-5.0d0*t2)/(2.0d0*(t1-t2))
c     1                 + phi2*(t1*t1)/(t2*(t1-t2)))
c               ELSE
c                  b=phi(l)*(1.0d0/dsqrt(t(l)))+
c     1             phi(l-1)*(1.0d0/dsqrt(t(l-1)))
c                  h=t(l)-t(l-1)
c                  sum=sum+b*h*0.5d+0
c               ENDIF
c            ENDDO
c            IF(k==1.AND.m==1)THEN
c               sigma0_1=sum     !normalizzo a sigma_1(0)
c            ENDIF
c            sigma_i(m,k)=sum*(1.0d0/(1.0d0+R*R)) 
c            si_s0(m,k)=sigma_i(m,k)/sigma0_1
c         ENDDO
c      ENDDO
      
    
      return
      end
      
c     Subroutine per il calcolo della velocità quadratica media (unità 3sigma1)
      
      SUBROUTINE v2mean(fcnqd,fcnqd2,fu1,fcnqd1)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/VEL/v2_i(30,8000),v2(8000)
      COMMON/PARFCN/WWWW
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse

      estr1=0.0d0
      DO l=1, nfn-1
         sum =0.0d0
         sum2=0.0d0
         qd=0.0d0
         qd2=0.0d0
         qd1=0.0d0
         w1=w(l)
         DO k=1,Nmasse! 10
            wwww=w1*ym(k)
            estr2=w1*ym(k)
            qd=0.0d0
            call intgau(fcnqd,qd)
            v2_i_winf=qd
            call intgau(fcnqd2,qd2)
            v2_i(k,l)=(0.4d0)*(qd2/v2_i_winf)
            if(l==nfn)v2_i(k,l)=v2_i(k,l-1)
c            qd1=0.0d0
c            wwww=w1*ym(k)
c            estr2=w1*ym(k)
c            call intgau(fu1,qd1)
            sum =sum + Aalfa(k)*(ym(k)**(-2.5d0))*qd2
            sum2=sum2+Aalfa(k)*(ym(k)**(-1.5d0))*qd!*qd1!*dexp(wwww)
            
         ENDDO
         v2(l) = 0.4d0*sum/sum2!(2.0d0*sum)/(3.0d0*sum2)
c         PRINT*,v2(l)
      ENDDO
      v2(nfn)=v2(nfn-1)
      wwww=w1
      estr2=wwww
      return
      end
   
c	la funzione da integrare (ce ne possono essere più di una)
c
c
        function fcnqd (x)
	IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/PARFCN/WWWW
        fcnqd = (dexp(wwww-x))*(x**1.5d0)
c        fcnqd= dexp(-(x*x))
	return
        end
c
      function fcnqd2(t)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARFCN/WWWW
        fcnqd2 = (dexp(wwww-t))*(t**2.5d0) 
	return
        end
c	la subroutine con la equazione differenziale
c
	SUBROUTINE FCT (X,Y,DERY,fcnqd)
	IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/PARFCN/WWWW
        COMMON/ESTREMI/ESTR1,ESTR2
        COMMON/VALINI/W0,RHO0
        COMMON/COSTANTI/PAI
	COMMON/FCTOUTP/RPP
	DIMENSION Y(2),DERY(2)
        COMMON/TEST/C,D,P,psi0(30)
        wwww = y(1)
        estr2= y(1)

        if(wwww.gt.0.d0) go to 1
        rpp = 0.d0
        go to 2
 1      continue
c       call intgau (fcnqd,qd)
        call rho(fcnqd,sum)
        rpp = sum/rho0
 2      continue

	dery(1) = y(2)
	dery(2) = -(2.d0/x)*y(2)-9.d0*rpp

	return
	end

      SUBROUTINE segregation
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      COMMON/SEGR/xMidr(30,8000),xNidr(30,8000)
      COMMON/SEGR2/xMiMtdr(30,8000),xNiNtdr(30,8000)
      DIMENSION xMitot(8000),xNitot(8000)
      COMMON/SEGR3/rhoidr(30,8000),xn_idr(30,8000)
      COMMON/SEGR4/rhoirhotdr(30,8000),xn_in_tdr(30,8000)
      DIMENSION rhotdr(8000),xntdr(8000)

c      DO i=1,10
c         base1=0.0d0
c         base2=0.0d0
c         h=0.0d0
c         somma1=0.0d0
c         somma2=0.0d0
c         DO j=1,nfn-1
c            base1=rhoi(i,j+1)*x(j+1)*x(j+1)+rhoi(i,j)*x(j)*x(j)
c            base2=base1/ym(i)
c            h=x(j+1)-x(j)
c            xMidr(i,j)=base1*h*0.5d0
c            xNidr(i,j)=base2*h*0.5d0 
c         ENDDO
c         somma1=somma1+xMidr(i,j)
c            somma2=somm2+xNidr(i,j)
c      ENDDO
c      DO j=1, nfn-1
c         xMiMtdr(i,j)=xMidr(i,j)/somma1
c         xNiNtdr(i,j)=xNidr(i,j)/somma2
c      ENDDO
      
      DO j=1,nfn
         somma1=0.0d0
         somma2=0.0d0
         base1=0.0d0
         base2=0.0d0
         h=0.0d0
         DO i=1,Nmasse!10
           base1=rhoi(i,j)*x(j)*x(j)+rhoi(i,j+1)*x(j+1)*x(j+1)
            base2=base1/ym(i)
            h=x(j+1)-x(j)
            xMidr(i,j)=base1*h*0.5d0
            xNidr(i,j)=base2*h*0.5d0
            somma1=somma1+xMidr(i,j)
            somma2=somma2+xNidr(i,j)
         ENDDO
         xMitot(j)=somma1
         xNitot(j)=somma2
      ENDDO
     
      

      DO j=1,nfn
         DO i=1,Nmasse!10
            xMiMtdr(i,j)=xMidr(i,j)/xMitot(j)
            xNiNtdr(i,j)=xNidr(i,j)/xNitot(j)
         ENDDO
      ENDDO

c     Calcolo rho_i(r)/rho_tot(r) (e stesso con densità numerica) per vedere segregazione
      DO j=1,nfn-1
         somma3=0.0d0
         somma4=0.0d0
         DO i=1,Nmasse
            somma3=somma3+rhoi(i,j)
            somma4=somma4+rhoi(i,j)/ym(i)
         ENDDO
         rhotdr(j)=somma3
         xntdr(j)=somma4
      ENDDO

      DO j=1,nfn-1
         DO i=1,Nmasse
            rhoirhotdr(i,j)=rhoi(i,j)/rhotdr(j)
            xn_in_tdr(i,j)=rhoi(i,j)/(ym(i)*xntdr(j))
         ENDDO
      ENDDO
      
            


         
      RETURN
      END
      
c     Subroutine per il calcolo della massa relativa mu_i e beta_i=r_t/r_c,i
      
      SUBROUTINE reltotmass (fcnqd)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/PARFCN/WWWW
      COMMON/VALINI/W0,RHO0
      COMMON/COSTANTI/PAI
      COMMON/MULTIMASS/Aalfa(30), ym(30),Nmasse
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/RESULTS2/rhoi(30,8000),ymu(30),csi(8000),beta(30),conc
      
c     Calcolo densità per ogni classe di massa e step

      qd=0.d+0
      z1=0.0d+0
      z2=0.0d+0
      z3=0.0d+0
      z4=0.0d+0
      DO m=1, nfn
         w_1 = w(m)
         DO l=1,Nmasse! 10
            wwww=w_1*ym(l)
            estr2=dsqrt(w_1*ym(l))
c            estr2=wwww
c            call intgau(fcnqd,qd)
c      PRINT*,qd
c            qd=DERF(estr2)*(dsqrt(pai)/2.0d0)
c     ! Calcolo densità con formula 7 (senza termine exp(C))
c              
            z1=Aalfa(l)*(1.0d0/dsqrt(ym(l)))
            z2=((2.d0*PAI)**(1.5d0))*dexp(wwww)*(DERF(estr2))
            z3=-((2.d0)**(2.5d0))*PAI*estr2
            z4=-(((2.d0)**(3.5d0))*PAI*((wwww)**(1.5d0)))/3.0d0
            rhoi(l,m)= z1*(z2+z3+z4)
              
c     ! Calcolo densità con integrale su w_i
              
c              rhoi(l,m) =qd*Aalfa(l)*(1.0d0/dsqrt(ym(l)))
c          rhoi(l,m)=qd*Aalfa(l)*(1.d0/sqrt(ym(l)))*dexp(-w0*ym(l)) !variante
              IF(m==nfn) rhoi(l,m)=0.0d0     
           END DO
        END DO
        
c     Ciclo per il calcolo delle variabili mu_i e beta_i
        
        DO l=1,Nmasse! 10
         somma = 0.0d+0
         DO n=1, nfn-1
            b=(rhoi(l,n+1)*csi(n+1)*csi(n+1))+(rhoi(l,n)*csi(n)*csi(n))
            h=csi(n+1) - csi(n)
            somma = somma + b*h*0.5d+0
         END DO
         IF(l==1) Den=somma
         ymu(l)=somma/Den
      END DO
      wwww=w_1
      estr2=wwww
      return
      end

c     Subroutine per il calcolo della densità
      
      SUBROUTINE rho (fcnqd,sum)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ESTREMI/ESTR1,ESTR2
      COMMON/MULTIMASS/Aalfa(30),ym(30),Nmasse
      COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
      COMMON/PARFCN/WWWW
      COMMON/COSTANTI/PAI
      COMMON/VALINI/W0,RHO0
      COMMON/TEST/D,C,P,psi0(30)
      
      estr1 = 0.0d+0
      sum = 0.0d+0
      w_1 = wwww
      qd=0.0d+0
      z1=0.0d+0
      z2=0.0d+0
      z3=0.0d+0
      z4=0.0d+0
      DO l=1,Nmasse! 10
         wwww = w_1*ym(l)
         estr2 =dsqrt(w_1*ym(l))
c        estr2=w_1*ym(l)
c         call intgau (fcnqd, qd)       
c         sum=sum+Aalfa(l)*qd*(1.0d0/dsqrt(ym(l)))*dexp(-(w0*ym(l)))
         
         z1=Aalfa(l)*(1.0d0/dsqrt(ym(l)))
         z2=((2.0d0*PAI)**(1.5d0))*dexp(wwww)*(DERF(estr2))
         z3=-((2.0d0)**(2.5d0))*PAI*dsqrt(wwww)
         z4=-((2.0d0)**(3.5d0))*PAI*((wwww)**(1.5d0))/3.0d0
         sum = sum + z1*(z2+z3+z4)
         
c     IF(w_1==w0)THEN 
c     psi0(l)= qd
c     psi0(l)=z1*(z2+z3+z4)
c     ENDIF
c     sum = sum + ((z1*(z2+z3+z4))/psi0(l))*Aalfa(l)
c     sum=sum+(qd/psi0(l))*Aalfa(l) !con sum_i Alfa(i) =1
c     PRINT*,psi_0, sum
         
      END DO
      wwww=w_1
      estr2=wwww
      return
      end
      
c     la subroutine dei valori di output per l'integrazione
      
	SUBROUTINE OUTP (XX,Y,DERY,IHLF,NDIM,PRMT)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION Y(2),DERY(2),PRMT(5)
        COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
	COMMON/FCTOUTP/RPP

	if (y(1).gt.0.d0) goto 2

 1      PRMT(5)=1.d0
        goto 3

 2      xrif = x(nfn) + prmt(3)
        if(xx.lt.xrif) return
 
 3      nfn=nfn+1
	if (nfn.eq.8000) PRINT*,'TROPPI PASSI'
	if (nfn.eq.8000) PRMT(5)=1.00

        
	    x(nfn)  = xx
	    w(nfn)  = y(1)
           dw(nfn)  = y(2)
	 rapp(nfn)  = rpp

	return
	end

c	la subroutine per il calcolo degli integrali
c	
        subroutine INTGAU (fcn,area)
C
C       INTEGRA CON IL METODO DI GAUSS-LEGENDRE LA FUNZIONE FCN TRA A E B
C       I PARAMETRI DELLA FCN SONO PASSATI CON UN COMMON TRA PROGRAMMA 
C       CHIAMANTE LA SUBROUTINE INTGAU E LA FUNCTION FCN(X). 
C       I LIMITI DI INTEGRAZIONE CON UN COMMON TRA PROGRAMMA CHIAMANTE E 
C       LE SUBROUTINES GAUS10.FOR, GAUS20.FOR, GAUS40.FOR, GAUS80.FOR.
C       
C
	IMPLICIT REAL*8 (A-H,O-Z)
        common/gauss/indgau
C
        goto (1,2,3,4) ,indgau
 1	call GAUS20 (fcn,area)
        return
 2	call GAUS40 (fcn,area)
        return
 3	call GAUS80 (fcn,area)
        return
 4	call GAUS96 (fcn,area)
        return
	end

c
c	la subroutine dhpcg che è il cuore del programma

      SUBROUTINE DHPCG(PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX,FCNQD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PRMT(5),Y(2),DERY(2),AUX(16,2)
      N=1
      IHLF=0
      X=PRMT(1)
      H=PRMT(3)
      PRMT(5)=0.D0
      DO 1 I=1,NDIM
         AUX(16,I)=0.D0
         AUX(15,I)=DERY(I)
         AUX(1,I)=Y(I)
 1    CONTINUE
      IF(H*(PRMT(2)-X).lt.0.0d0)THEN
         GOTO 3
      ELSE IF(H*(PRMT(2)-X).eq.0.0d0)THEN
         GOTO 2
      ELSE
         GOTO 4
      ENDIF
 2    IHLF=12
      GOTO 4
 3    IHLF=13
 4    CONTINUE
      CALL FCT(X,Y,DERY,fcnqd)
      CALL OUTP (X,Y,DERY,IHLF,NDIM,PRMT)
      IF((PRMT(5).lt.0.0d0).OR.(PRMT(5).gt.0.0d0))THEN
        GOTO 6
      ELSE
         GOTO 5
      ENDIF
      
 5    IF(IHLF.le.0.0d0)THEN
        GOTO 7
      ELSE
         GOTO 6
      ENDIF
 6    RETURN
 7    DO 8 I=1,NDIM
         AUX(8,I)=DERY(I)
 8    CONTINUE
      ISW=1
      GOTO 100
 9    X=X+H
      DO 10 I=1,NDIM
         AUX(2,I)=Y(I)
 10   CONTINUE
 11   IHLF=IHLF+1
      X=X-H
      DO 12 I=1,NDIM
         AUX(4,I)=AUX(2,I)
 12   CONTINUE
      H=.5D0*H
      N=1
      ISW=2
      GOTO 100
 13   X=X+H
      CALL FCT(X,Y,DERY,fcnqd)
      N=2
      DO 14 I=1,NDIM
         AUX(2,I)=Y(I)
         AUX(9,I)=DERY(I)
 14   CONTINUE
      ISW=3
      GOTO 100
 15   DELT=0.D0
      DO 16 I=1,NDIM
         DELT=DELT+AUX(15,I)*DABS(Y(I)-AUX(4,I))
 16   CONTINUE
      DELT=.066666666666666667D0*DELT
      IF((DELT-PRMT(4)).le.0.0d0)THEN
         GOTO 19
      ELSE
         GOTO 17
      ENDIF
      
 17   IF((IHLF-10).lt.0.0d0)THEN
         GOTO 11
      ELSE
         GOTO 18
      ENDIF
      
 18   IHLF=11
      X=X+H
      GOTO 4
 19   X=X+H
      CALL FCT(X,Y,DERY,fcnqd)
      DO 20 I=1,NDIM
         AUX(3,I)=Y(I)
         AUX(10,I)=DERY(I)
 20   CONTINUE
      N=3
      ISW=4
      GOTO 100
 21   N=1
      X=X+H
      CALL FCT(X,Y,DERY,fcnqd)
      X=PRMT(1)
      DO 22 I=1,NDIM
	AUX(11,I)=DERY(I)
	Y(I)=AUX(1,I)+H*(.375D0*AUX(8,I)+.7916666666666667D0*AUX(9,I)
     1   -.20833333333333333D0*AUX(10,I)+.041666666666666667D0*DERY(I))
 22   CONTINUE
 23   X=X+H
      N=N+1
      CALL FCT(X,Y,DERY,fcnqd)
      CALL OUTP (X,Y,DERY,IHLF,NDIM,PRMT)
      IF((PRMT(5).lt.0.0d0).OR.(PRMT(5).gt.0.0d0))THEN
         GOTO 6
      ELSE
         GOTO 24
      ENDIF
      
 24   IF((N-4).lt.0.0d0)THEN
         GOTO 25
      ELSE
         GOTO 200
      ENDIF
      
 25   DO 26 I=1,NDIM
         AUX(N,I)=Y(I)
         AUX(N+7,I)=DERY(I)
 26   CONTINUE
      IF((N-3).lt.0.0d0)THEN
         GOTO 27
      ELSEIF( (N-3).eq.0.0d0)THEN
         GOTO 29
      ELSE
         GOTO 200
      ENDIF
 27   DO 28 I=1,NDIM
	DELT=AUX(9,I)+AUX(9,I)
	DELT=DELT+DELT
        Y(I)=AUX(1,I)+.33333333333333333D0*H*(AUX(8,I)+DELT+AUX(10,I))
 28   CONTINUE
      GOTO 23
 29   DO 30 I=1,NDIM
         DELT=AUX(9,I)+AUX(10,I)
         DELT=DELT+DELT+DELT
         Y(I)=AUX(1,I)+.375D0*H*(AUX(8,I)+DELT+AUX(11,I))
 30   CONTINUE
      GOTO 23
 100  DO 101 I=1,NDIM
         Z=H*AUX(N+7,I)
         AUX(5,I)=Z
         Y(I)=AUX(N,I)+.4D0*Z
 101  CONTINUE
      Z=X+.4D0*H
      CALL FCT(Z,Y,DERY,fcnqd)
      DO 102 I=1,NDIM
         Z=H*DERY(I)
         AUX(6,I)=Z
      Y(I)=AUX(N,I)+.29697760924775360D0*AUX(5,I)+.15875964497103583D0*Z
 102  CONTINUE
	Z=X+.45573725421878943D0*H
	CALL FCT(Z,Y,DERY,fcnqd)
	DO 103 I=1,NDIM
	Z=H*DERY(I)
	AUX(7,I)=Z
	Y(I)=AUX(N,I)+.21810038822592047D0*AUX(5,I)-3.0509651486929308D0*
     1       AUX(6,I)+3.8328647604670103D0*Z
 103  CONTINUE
      Z=X+H
      CALL FCT(Z,Y,DERY,fcnqd)
      DO 104 I=1,NDIM
       Y(I)=AUX(N,I)+.17476028226269037D0*AUX(5,I)-.55148066287873294D0*
	1AUX(6,I)+1.2055355993965235D0*AUX(7,I)+.17118478121951903D0*
	2H*DERY(I)
 104  CONTINUE
      GOTO(9,13,15,21),ISW
 200  ISTEP=3
 201  IF(((N-8).lt.0.0d0).OR.((N-8).gt.0.0d0))THEN
         GOTO 204
      ELSE
         GOTO 202
      ENDIF
      
 202  DO 203 N=2,7
         DO  I=1,NDIM
            AUX(N-1,I)=AUX(N,I)
            AUX(N+6,I)=AUX(N+7,I)
         ENDDO
 203  CONTINUE
         N=7
 204     N=N+1
         DO 205 I=1,NDIM
            AUX(N-1,I)=Y(I)
            AUX(N+6,I)=DERY(I)
 205     CONTINUE
         X=X+H
 206     ISTEP=ISTEP+1
         DO 207 I=1,NDIM
	DELT=AUX(N-4,I)+1.3333333333333333D0*H*(AUX(N+6,I)+AUX(N+6,I)-
	1AUX(N+5,I)+AUX(N+4,I)+AUX(N+4,I))
	Y(I)=DELT-.9256198347107438D0*AUX(16,I)
        AUX(16,I)=DELT
 207    CONTINUE
	CALL FCT(X,Y,DERY,fcnqd)
	DO 208 I=1,NDIM
	DELT=.125D0*(9.D0*AUX(N-1,I)-AUX(N-3,I)+3.D0*H*(DERY(I)+AUX(N+6,I)
	1+AUX(N+6,I)-AUX(N+5,I)))
	AUX(16,I)=AUX(16,I)-DELT
        Y(I)=DELT+.07438016528925620D0*AUX(16,I)
 208    CONTINUE
	DELT=0.D0
	DO 209 I=1,NDIM
           DELT=DELT+AUX(15,I)*DABS(AUX(16,I))
 209    CONTINUE
        IF((DELT-PRMT(4)).lt.0.0d0)THEN
           GOTO 210
        ELSE
           GOTO 222
        ENDIF
        
210	CONTINUE
        CALL FCT(X,Y,DERY,fcnqd)
        CALL OUTP (X,Y,DERY,IHLF,NDIM,PRMT)
	IF((PRMT(5).lt.0.0d0).OR.(PRMT(5).gt.0.0d0))THEN
           GOTO 212
        ELSE
           GOTO 211
        ENDIF
 211    IF((IHLF-11).lt.0.0d0)THEN
           GOTO 213
        ELSE
           GOTO 212
        ENDIF
 212    RETURN
 213    IF(H*(X-PRMT(2)).lt.0.0d0)THEN
           GOTO 214
        ELSE
           GOTO 212
        ENDIF
 214    IF(DABS(X-PRMT(2))-.1D0*DABS(H).lt.0.0d0)THEN
           GOTO 212
        ELSE
           GOTO 215
        ENDIF
 215    IF(DELT-.02D0*PRMT(4).le.0.0d0)THEN
           GOTO 216
        ELSE
           GOTO 201
        ENDIF
 216    IF(IHLF.le.0.0d0)THEN
           GOTO 201
        ELSE
           GOTO 217
        ENDIF
 217    IF(N-7.lt.0.0d0)THEN
           GOTO 201
        ELSE
           GOTO 218
        ENDIF
 218    IF(ISTEP-4.lt.0.0d0)THEN
           GOTO 201
        ELSE
           GOTO 219
        ENDIF
 219    IMOD=ISTEP/2
	IF(ISTEP-IMOD-IMOD.ne.0.0d0)THEN
           GOTO 201
        ELSE
           GOTO 220
        ENDIF
220	H=H+H
	IHLF=IHLF-1
	ISTEP=0
	DO 221 I=1,NDIM
	AUX(N-1,I)=AUX(N-2,I)
	AUX(N-2,I)=AUX(N-4,I)
	AUX(N-3,I)=AUX(N-6,I)
	AUX(N+6,I)=AUX(N+5,I)
	AUX(N+5,I)=AUX(N+3,I)
	AUX(N+4,I)=AUX(N+1,I)
	DELT=AUX(N+6,I)+AUX(N+5,I)
	DELT=DELT+DELT+DELT
	AUX(16,I)=8.962962962962963D0*(Y(I)-AUX(N-3,I))
	1-3.3611111111111111D0*H*(DERY(I)+DELT+AUX(N+4,I))
 221    CONTINUE
        GOTO 201
222	IHLF=IHLF+1
	IF(IHLF-10.le.0.0d0)THEN
           GOTO 223
        ELSE
           GOTO 210
        ENDIF
 223    H=.5D0*H
	ISTEP=0
	DO 224 I=1,NDIM
	Y(I)=.390625D-2*(8.D1*AUX(N-1,I)+135.D0*AUX(N-2,I)+4.D1*AUX(N-3,I)
	1+AUX(N-4,I))-.1171875D0*(AUX(N+6,I)-6.D0*AUX(N+5,I)-AUX(N+4,I))*H
	AUX(N-4,I)=.390625D-2*(12.D0*AUX(N-1,I)+135.D0*AUX(N-2,I)+
	1108.D0*AUX(N-3,I)+AUX(N-4,I))-.0234375D0*(AUX(N+6,I)+
	218.D0*AUX(N+5,I)-9.D0*AUX(N+4,I))*H
	AUX(N-3,I)=AUX(N-2,I)
        AUX(N+4,I)=AUX(N+5,I)
 224    CONTINUE
	X=X-H
	DELT=X-(H+H)
	CALL FCT(DELT,Y,DERY,fcnqd)
	DO 225 I=1,NDIM
	AUX(N-2,I)=Y(I)
	AUX(N+5,I)=DERY(I)
        Y(I)=AUX(N-4,I)
 225    CONTINUE
	DELT=DELT-(H+H)
	CALL FCT(DELT,Y,DERY,fcnqd)
	DO 226 I=1,NDIM
	DELT=AUX(N+5,I)+AUX(N+4,I)
	DELT=DELT+DELT+DELT
	AUX(16,I)=8.962962962962963D0*(AUX(N-1,I)-Y(I))
	1-3.3611111111111111D0*H*(AUX(N+6,I)+DELT+DERY(I))
        AUX(N+3,I)=DERY(I)
 226    CONTINUE
	GOTO 206
	END                                                              
