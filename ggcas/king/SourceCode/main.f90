      PROGRAM MAIN
        IMPLICIT REAL*8(A-H,O-Z)
        REAL*4 temp1,temp2
        COMMON/PARFCN/WWWW
        COMMON/ESTREMI/ESTR1,ESTR2
        COMMON/VALINI/W0,RHO0
        COMMON/GAUSS/INDGAU
        COMMON/COSTANTI/PAI
        COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
        COMMON/VEL/v2(8000)
        DIMENSION Y(2),DERY(2),PRMT(5),AUX(16,2)
        EXTERNAL FCT,OUTP, INTGAU!, SURFDENS
        external :: fcnqd, fcnqd1, fcnqd2
        external :: fu1, fu2, fu3, fu4, funzerr, funzione
       ! external :: FU1,FU2,FU3,FU4,FCNQD1,funzerr, funzione
        COMMON/SURF/t(8000),phi(8000),sigma(8000),psi(8000),S_S0(8000)
        COMMON/SURF2/csi(8000)
        COMMON/CV/Cv0(8000),Cv(8000),Cvqt(8000),Ctot0,Ctot1,Ctot2
        COMMON/CV2/xMhat,xnhat(8000),x0Cv(10),xi0Cv(10),potch(8000)
        COMMON/PARFU2/qq
        DIMENSION xN(8000)
        DIMENSION phi_g(8000)
        DIMENSION phi_g2(8000),eg2(8000),xmu_r(8000),secpezz(8000)
        COMMON/ENERGY/xk(8000),eg(8000),ephi(8000),xkt,egt,ephit,Etot,Vir2,et_r(8000),Etot1
        COMMON/TEST/C,D
        COMMON/CC/yEt_Mv02(10000),xs2_v02(10000),fw0,xsave,ysave,w0save
        COMMON/Mrho0hat/rho0hat,xMcap,xmu,xnhat0,conc
        COMMON/SURFEKIN/uk(8000),sigekin(8000),Sk_Sk0(8000)
        DIMENSION parw0(6)!135)
   DATA PARW0/6.18d0,6.22d0,7.41d0,6.48d0,&
   2.09d0,6.29d0/

   DATA NAMESID/"NGC 4372","NGC 5139","NGC 6121","NGC 6656",&
   "Pal 5","GLIMPSE02"/


   !PRINT*,"Analyzed the following GCs: ",NAMESID
   OPEN(UNIT=1,FILE='params.dat',STATUS='replace') !,ACCESS='append')
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='profiles.dat',STATUS='replace')!,ACCESS='append')
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='CvNtK.dat',STATUS='replace')
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='Cv.dat',STATUS='replace')
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='phi.dat',STATUS='replace')
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='Er.dat',STATUS='replace')
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='Etot.dat',STATUS='replace')
   WRITE(1,303)'w0 K Egr Eeff Etot s2_v02 Et_Mv02 Mcap&
   &Vir conc Etot1'
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='CalCurve.dat',STATUS='replace')
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='x0Cv.dat',STATUS='replace')
   WRITE(1,303)'w0 x1 x2 x3 xi1 xi2 xi3 conc c'
   CLOSE(UNIT=1,STATUS='keep')
   OPEN(UNIT=1,FILE='Skin.dat',STATUS='replace')
   ! WRITE(1,303)'x xi w Sk'
   CLOSE(UNIT=1,STATUS='keep')



   temp1=0.0d0
   temp2=0.0d0
   call second(temp1)


   PRINT*,'IW1,IW2,IWSTEP  1,200,1  W0 = PARW0(IW)'
   IW1=1
   IW2=6!135!160    !221
   IWSTEP=1
!     Costruisco un ciclo che va da w0=0.1 a w0=45
!        PARW0(1)=0.1d0
!        DO i=2,IW2
!           PARW0(i)=PARW0(i-1)+0.1d0
!        ENDDO
!	l'indice che precisa che subroutine per gli integrali si vuole usare
!	il valore 3 corrisponde a 80 punti
   indgau = 4
!	l'errore relativo si decide a priori e ci dï¿½ la precisione voluta
   relerr = 1.d-9
   estr1  = 0.d0
   pai    = 3.141592654d0
      w0 = 6.19d0 !! parw0(iw)
      w(1) = w0
      dw(1) = 0.d0
      wwww = w0
      estr2= w0
      call prova(funzione, "PROVA", 2.0D0)
      call intgau(fcnqd,qd)
      print*,"PROVA"
      rho0 = qd
         rapp(1) = qd/rho0
   !       INIZIALIZZAZIONE (per evitare la singolaritï¿½ iniziale)
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
      call intgau(fcnqd,qd)
      rapp(2) = qd/rho0
!       FINE INIZIALIZZAZIONE
      nfn   = 2
      passo = x(nfn)/10.0d0!10!2
      xmax  = x(nfn)*10.0d0!10!2
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
      CALL DHPCG(PRMT,Y,DERY,2,IHLF,FCT,OUTP,AUX,FCNQD)
      if(ihlf.lt.11) goto 1122
      print*,'bisezioni>10',nfn,w0!,x(nfn)
      passo = passo/10.d0!10!200
      xmax  = x(nfn)+(xmax-x(nfn))/10.d0!10!200
      goto 1111
 1122   continue
      if(nfn.eq.8000) goto 3333
      if(w(nfn).lt.0.d0)   GO TO 1234
      passo = passo*10.d0!10!2
      xmax  = xmax *10.d0!10!2
      goto 1111
 1234   continue
      xtest = (dw(nfn-1)*x(nfn-1)-w(nfn-1))/dw(nfn-1)
   !	controllo precisione sulla routine di uscita per determinare il raggio
      if(dabs((xtest-x(nfn))/xtest).lt.1.d-9) goto 2222
      xmax  =  x(nfn)
      passo = (x(nfn)-x(nfn-1))/10.d0!10!200
      nfn = nfn - 1
      goto 1111
 2222	continue
      conc  =  (xtest+x(nfn))/2.d0
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
      call surfdens(fcnqd)
      call v2mean(fcnqd,fcnqd2)
      call calorespecifico(fcnqd1,fu1,fu2,fu3,fu4)!,funzerr)
      xKb=1.0d0
      Ct0_NtK=xKb*Ctot0
      Ct1_NtK=xKb*Ctot1!/xmu
      Ct2_NtK=xKb*Ctot2       !/xmu
      call energie(fu1,fu3,fcnqd1)
      IF(w0.eq.8.20d0)THEN
         xm_k=0.54274d0
         DO k=1,nfn
            phi_g(k)=(-C-w(k))!/xm_k !=phi(r)/k\theta
         ENDDO
      ENDIF
      xmu_r(1)=0.d0
      DO l=1,nfn-1
         b1=0.d0
         h1=0.d0
         sum1=0.d0
         IF(l.gt.1)THEN
            DO i=2,l
               b1=rapp(i-1)*x(i-1)*x(i-1)+rapp(i)*x(i)*x(i)
               h1=x(i)-x(i-1)
               sum1=sum1+b1*h1*0.5d0
            ENDDO
         ENDIF
         xmu_r(l)=sum1
         b2=0.d0
         h2=0.d0
         sum2=0.d0
         DO m=l,nfn-1
            b2=rapp(m+1)*x(m+1)+rapp(m)*x(m)
            h2=x(m+1)-x(m)
            sum2=sum2+b2*h2*0.5d0
         ENDDO
         secpezz(l)=sum2

      phi_g2(l)=-(9.d0/xm_k)*((sum1/x(l))+sum2)
      ENDDO
      DO k=2,nfn-1
         eg2(k)=xnhat(k)*(-9.d0)*((xmu_r(k)/x(k)))!+secpezz(k))!*(-C-w(k))
      ENDDO
      eg2(1)=0.d0!0.5d0*xnhat(1)*(-9.d0)*(secpezz(1))
      eg2(nfn)=0.0d0
      b=0.0d0
      h=0.d0
      sum=0.0d0
      DO l=1,nfn-1
         b=eg2(l+1)*csi(l+1)*csi(l+1)+eg2(l)*csi(l)*csi(l)
         h=csi(l+1)-csi(l)
         sum=sum+b*h*0.5d0
      ENDDO
      egt2=sum*4.d0*pai
      Shat=IW*0.01d0
      OPEN(UNIT=1,FILE='params.dat',STATUS='old',ACCESS='append')
      WRITE(1,301) w0,conc,xlogc,xmu
      CLOSE(UNIT=1,STATUS='keep')
      OPEN(UNIT=1,FILE='profiles.dat',STATUS='old',ACCESS='append')
      WRITE(1,303)'x xi w rho/rho0 SD v2 logc=',xlogc
         DO jj=1,nfn
            WRITE(1,306) x(jj),csi(jj),w(jj),rapp(jj),s_s0(jj),v2(jj)
         ENDDO
      CLOSE(UNIT=1,STATUS='keep') !,dispose='keep')
      IF((w0==0.8d0).OR.(w0==1.35d0).OR.(w0==2.0d0).OR.&
      (w0==3.5d0).OR.(w0==5.d0).OR.(w0==6.0d0).OR.&
      (w0==8.0d0).OR.(w0==8.2d0).OR.(w0==40.0d0))THEN
         OPEN(UNIT=1,FILE='Cv.dat',STATUS='old',ACCESS='append')
         WRITE(1,303)'x xi w  dw Cv0 Cv Cvqt w0=',w0
         DO m=1,nfn
            WRITE(1,302)x(m),csi(m),w(m),dw(m),Cv0(m)/xnhat(m),&
            Cv(m)/xnhat(m), Cvqt(m)/xnhat(m)
         ENDDO
         CLOSE(UNIT=1,STATUS='keep')
      ENDIF
      OPEN(UNIT=1,FILE='CvNtK.dat',status='old',access='append')
      WRITE(1,302) w0,Ct0_NtK,Ct1_NtK,Ct2_NtK,conc,x0Cv,xi0Cv
      CLOSE(UNIT=1,STATUS='keep')
      OPEN(UNIT=1,FILE='x0Cv.dat',status='old',access='append')
      WRITE(1,302) w0,x0Cv(1),x0Cv(2),x0Cv(3),xi0Cv(1),xi0Cv(2),xi0Cv(3)&
      ,conc,xlogc
      CLOSE(UNIT=1,STATUS='keep')
      IF(w0==8.20d0)THEN
         OPEN(UNIT=1,FILE='phi.dat',status='old',access='append')
         DO k=1,nfn
            WRITE(1,302)x(k),csi(k),w(k),dw(k),phi_g(k)
         ENDDO
         CLOSE(UNIT=1,STATUS='keep')
      ENDIF
      espo=4.d0/3.d0
      s2_v02=(1.d0/(xMcap**espo))
      Et_Mv02=s2_v02*(Etot1/xMcap)
      xkt=s2_v02*(xkt/xMcap)
      egt=s2_v02*(egt/xMcap)
      ephit=s2_v02*(ephit/xMcap)
      OPEN(UNIT=1,FILE='Er.dat',status='old',access='append')
      OPEN(UNIT=2,FILE='Etot.dat',status='old',access='append')
      IF((w0==0.8d0).OR.(w0==1.35d0).OR.(w0==2.0d0).OR.&
      (w0==3.5d0).OR.(w0==5.d0).OR.(w0==6.0d0).OR.&
      (w0==8.0d0).OR.(w0==8.2d0).OR.(w0==40.0d0))THEN
         DO k=1,nfn
            WRITE(1,302)x(k),csi(k),w(k),xk(k),eg(k),ephi(k),et_r(k)
         ENDDO
      ENDIF
      WRITE(2,302)w0,xkt,egt,ephit,Etot,s2_v02,Et_Mv02,xMcap,& !hat,!(shat*Etot/(4.d0*pai*xMhat))
      Vir2,conc, Etot1
      CLOSE(UNIT=1,STATUS='keep')
      CLOSE(UNIT=2,STATUS='keep')
      call surf_dens_ekin(fu1)
      OPEN(UNIT=1,FILE='Skin.dat',status='old',access='append')
      WRITE(1,305)'NameID=',NAMESID(iw)
      WRITE(1,303)'w0=',w0
      WRITE(1,303)'x xi w Sk sk2'
      DO jj=1,nfn
         WRITE(1,302) x(jj),csi(jj),w(jj),Sk_Sk0(jj),sigekin(jj)
      ENDDO
      CLOSE(UNIT=1,STATUS='keep') !,dispose='keep')
301   FORMAT(1pd9.2,1pd10.3,1x,1p2d15.8)
302   FORMAT(1p11d16.8)
303   FORMAT(A70,1p3d15.8)
304   FORMAT(1p5d16.8)
305   FORMAT(A30,A30)
306   FORMAT(1p6e15.8)
 3333   CONTINUE
      CALL SECOND(TEMP2)
      T_imp=0.0d0
      T_imp=(temp2-temp1)
      PRINT*,'Tempo programma(s)=',T_imp,'T(min)=',T_imp/60.0
      DO j=1,2
         print *, char(7)
         call sleep(1)
      ENDDO
      STOP
END PROGRAM
