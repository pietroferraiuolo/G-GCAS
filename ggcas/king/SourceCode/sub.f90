   SUBROUTINE surf_dens_ekin(fcn)
   IMPLICIT REAL*8(A-H,O-Z)
   COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
   COMMON/PARFCN/WWWW
   COMMON/ESTREMI/ESTR1,ESTR2
   COMMON/SURF2/csi(8000)
   COMMON/SURFEKIN/uk(8000),sigekin(8000),Sk_Sk0(8000)
   DIMENSION yy(8000)
   external fcn
   u1=0.0d0
   DO l=1, nfn-1
      wwww=w(l)
      estr2=w(l)
      call intgau(fcn,u1)
      uk(l)=dexp(wwww)*u1
   ENDDO
   sum=0.0d0
   b=0.0d0
   h=0.0d0
   xiR=0.0d0
   DO i=1,nfn
      xiR=csi(i)
      DO k=i,nfn
         yy(k)=dsqrt((csi(k)**2.d0)-(xiR**2.d0))
         IF(k==i)THEN
            sum=0.0d0
         ELSE
            b=uk(k)+uk(k-1)
            h=yy(k)-yy(k-1)
            sum=sum+b*h*0.5d0
!               PRINT*,sum,b,h
         ENDIF
      ENDDO
      IF(i==1)THEN
         sigekin0=sum
         IF(sum==0)THEN
            PRINT*,'zero central Skin'
         ENDIF

      ENDIF
      sigekin(i)=sum!!*(1.0d0/(1.0d0+R*R))

      Sk_Sk0(i)=sigekin(i)/sigekin0
   END DO
   return
   end
   SUBROUTINE caloric_curve()
   IMPLICIT REAL*8(A-H,O-Z)
   COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
   COMMON/PARFCN/WWWW
   COMMON/ESTREMI/ESTR1,ESTR2
   COMMON/VALINI/W0,RHO0
   COMMON/CV2/xMhat,xnhat(8000),x0Cv(10),xi0Cv(10),potch(8000)
   COMMON/ENERGY/xk(8000),eg(8000),ephi(8000),xkt,egt,ephit,Etot,Vir2&
   ,et_r(8000),Etot1
   COMMON/CC/yEt_Mv02(10000),xs2_v02(10000),fw0,xsave,ysave,w0save
   DIMENSION xs(10000),ys(10000),xre(10000)
   real, DIMENSION(:), allocatable:: xre2
   real, DIMENSION(:), allocatable:: xsv
   real, DIMENSION(:), allocatable:: ysv
   fw0=Etot/xMhat
   xincr=7.0d0/1.d+4
   yEt_Mv02(1)=-2.0d0
   xincr2= 1.5d0/1.0d+4
   xs2_v02(1)=0.5d0
   xs(1)=0.0d0
   ys(1)=0.0d0
   DO i=2,10000
      xs2_v02(i)=xs2_v02(i-1)+xincr2        !\in[0.5,2]
      yEt_Mv02(i)=yEt_Mv02(i-1)+xincr !\in[-2,5]
   ENDDO
   xsave=0.0d0
   ysave=0.d0
   w0save=0.d0
   k=0
   DO i=1,10000
      xxx=0.0d0
      xxx=xs2_v02(i)
      ex2=0.0d0
      yyy=0.0d0
      xxxx=0.0d0
      DO j=1,10000
         yyy=yEt_Mv02(j)
         fw0=Etot/xMhat
         xxxx=yyy/fw0
         IF((abs((xxxx-xxx)/xxx)).lt.(abs((xs(i)-xxx)/xxx)).AND.&
         (abs((xxxx-xxx)/xxx).le.1.0d-7))THEN
               k=k+1
            xs(k)=xxxx
            ys(k)=yyy
            XRE(k)=abs((xs(k)-xxx)/xxx)
            ex2=1.0d0
         ENDIF
      ENDDO
   ENDDO
   allocate (xre2(1:k))
   allocate (xsv(1:k))
   allocate (ysv(1:k))
   DO i=1,k
      xre2(i)=xre(i)
      xsv(i)=xs(i)
      ysv(i)=ys(i)
   ENDDO
   xresave=MINVAL(xre2)          !xs(k)
   DO i=1,k
      IF(xresave-xre2(i).eq.0.0d0)THEN
         xsave=xsv(i)
         ysave=ysv(i)
         w0save=w0
         PRINT*,'x=',xsave,'y=',ysave
      ENDIF
   ENDDO
   RETURN
   END

   function fcnqd2(t)
   IMPLICIT REAL*8 (A-H,O-Z)
   COMMON/PARFCN/WWWW
   fcnqd2 = (dexp(wwww-t))*(t**2.5d0)
   return
   end

   SUBROUTINE v2mean(fcn1,fcn2)
   IMPLICIT REAL*8(A-H,O-Z)
   COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
   COMMON/VEL/v2(8000)
   COMMON/PARFCN/WWWW
   COMMON/ESTREMI/ESTR1,ESTR2
   EXTERNAL fcn1,fcn2
         DO l=1, nfn-1
         wwww=w(l)
         estr2=w(l)
         call intgau(fcn1,qd)
         v2_winf=qd
         call intgau(fcn2,qd2)
         v2(l)=(0.4d0)*(qd2/v2_winf)
      ENDDO
      v2(nfn)=v2(nfn-1)
      return
      end

SUBROUTINE surfdens(fcn)
   IMPLICIT REAL*8(A-H,O-Z)
   COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
   COMMON/SURF/t(8000),phi(8000),sigma(8000),psi(8000),S_S0(8000)
   COMMON/SURF2/csi(8000)
   COMMON/ESTREMI/estr1,estr2
   COMMON/PARFCN/WWWW
   DIMENSION yy(8000)
   EXTERNAL fcn
   qd=0.0d+0
   DO l=1,nfn-1
      wwww=w(l)
      estr2=wwww
      call intgau(fcn,qd)
      psi(l)=qd
   END DO
   sum=0.0d0
   b=0.0d0
   h=0.0d0
   xiR=0.0d0
   DO i=1,nfn
      xiR=csi(i)
      DO k=i,nfn
         yy(k)=dsqrt((csi(k)**2.d0)-(xiR**2.d0))
         IF(k==i)THEN
            sum=0.0d0
         ELSE
            b=psi(k)+psi(k-1)
            h=yy(k)-yy(k-1)
            sum=sum+b*h*0.5d0
         ENDIF
      ENDDO
      IF(i==1)THEN
         sigma0=sum
      ENDIF
      sigma(i)=sum
      s_s0(i)=sigma(i)/sigma0
   END DO
   return
   end

SUBROUTINE energie(f1,f3,fcn1)
   IMPLICIT REAL*8(A-H,O-Z)
   COMMON/ESTREMI/ESTR1,ESTR2
   COMMON/VALINI/W0,RHO0
   COMMON/GAUSS/INDGAU
   COMMON/COSTANTI/PAI
   COMMON/PARFCN/WWWW
   COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
   COMMON/SURF2/csi(8000)
   COMMON/PARFU2/qq
   COMMON/ENERGY/xk(8000),eg(8000),ephi(8000),xkt,egt,ephit,Etot,Vir2,&
   et_r(8000),Etot1
   COMMON/Mrho0hat/rho0hat,xMcap,xmu,xnhat0,conc
   DIMENSION rhat(8000)
   EXTERNAL f1,f3, fcn1
   z=0.0d0
   u1=0.0d0
   u3=0.0d0
   DO m=1,nfn-1
      z=w(m)
      wwww=w(m)
      estr2=wwww
      call intgau(f1,u1)
      wwww=w(m)
      estr2=wwww
      call intgau(f3,u3)
      xk(m)=dexp(z)*u1*4.d0*pai*dsqrt(2.d0)
      eg(m)=(2.d0/3.d0)*(dexp(z)*u1+0.4d0*(z**(2.5d0)))*dw(m)*x(m)*&
      4.d0*pai*dsqrt(2.d0)
      IF(abs(u3).gt.1.d-30)THEN
         ephi(m)=-xk(m)-dexp(z)*u3*4.d0*pai*dsqrt(2.d0)
      ELSE
         ephi(m)=0.0d0
         ENDIF
      et_r(m)=xk(m)+eg(m)+ephi(m)
   ENDDO
   xk(nfn)=0.0d0
   eg(nfn)=0.0d0
   ephi(nfn)=0.0d0
   et_r(nfn)=0.0d0

   wwww=w0
   estr2=w0
   qd=0.0d0
   call intgau(fcn1,qd)
   xnhat0=qd
   rho0hat=(4.d0*dsqrt(2.d0)*PAI)*qd
   xMcap=((9.d0/(4.d0*PAI))**(1.5d0))*(1.d0/dsqrt(rho0hat))*xmu
   DO i=1,nfn
      rhat(i)=dsqrt(9.0d0/(4.d0*PAI))*(1.d0/dsqrt(rho0hat))*x(i)
   ENDDO
   sum1=0.0d0
   sum2=0.0d0
   sum3=0.0d0
   b1=0.0d0
   b2=0.0d0
   b3=0.0d0
   h=0.0d0
   DO k=1,nfn-1
      b1=xk(k+1)*rhat(k+1)*rhat(k+1)+xk(k)*rhat(k)*rhat(k)
      b2=eg(k+1)*rhat(k+1)*rhat(k+1)+eg(k)*rhat(k)*rhat(k)
      b3=ephi(k+1)*rhat(k+1)*rhat(k+1)+ephi(k)*rhat(k)*rhat(k)
      h=rhat(k+1)-rhat(k)
      sum1=sum1+b1*h*0.5d0
      sum2=sum2+b2*h*0.5d0
      sum3=sum3+b3*h*0.5d0
   ENDDO
   xkt=0.0d0
   egt=0.0d0
   ephit=0.0d0
   xkt=sum1*4.d0*pai
   egt=sum2*4.d0*pai
   ephit=sum3*4.d0*pai
   Etot=0.0d0
   Etot=(sum1+sum2+sum3)
   Etot1=xkt+egt+ephit
   Vir=2.d0*xkt+egt
   Vir2=2.d0*xkt/(-egt)
   return
   end


SUBROUTINE calorespecifico(fcn1,f1,f2,f3,f4)
   IMPLICIT REAL*8(A-H,O-Z)
   COMMON/ESTREMI/ESTR1,ESTR2
   COMMON/VALINI/W0,RHO0
   COMMON/GAUSS/INDGAU
   COMMON/COSTANTI/PAI
   COMMON/PARFCN/WWWW
   COMMON/RESULTS/X(8000),W(8000),DW(8000),RAPP(8000),NFN
   COMMON/SURF2/csi(8000)
   COMMON/CV/Cv0(8000),Cv(8000),Cvqt(8000),Ctot0,Ctot1,Ctot2
   COMMON/CV2/xMhat,xnhat(8000),x0Cv(10),xi0Cv(10),potch(8000)
   COMMON/PARFU2/qq
   COMMON/TEST/C,D
   COMMON/Mrho0hat/rho0hat,xMcap,xmu,xnhat0,conc
   EXTERNAL fcn1,f1,f2,f3,f4,funzerr
   DO k=1,nfn
      Cv0(k)=0.0d0
      Cv(k)=0.0d0
      Cvqt(k)=0.0d0
      xnhat(k)=0.0d0
   ENDDO
   DO m=1,nfn-1
      z=0.0d0
      zz=0.0d0
      wwww=w(m)
      estr2=wwww
      z=w(m)
      qd=0.0d0
      u1=0.0d0
      u21=0.0d0
      u23=0.0d0
      u3=0.0d0
      u4=0.0d0
      zz=dsqrt(z)
      call intgau(fcn1,qd)
      xnhat(m)=qd
      wwww=w(m)
      estr2=wwww
      call intgau(f1,u1)
      wwww=w(m)
      estr2=wwww
      qq=0.5d0
      call intgau(f2,u21)

      qq=1.5d0
      wwww=w(m)
      estr2=wwww
      call intgau(f2,u23)

      potch(m)=u21*dexp(z)/qd

      wwww=w(m)
      estr2=wwww
      call intgau(f3,u3)

      wwww=w(m)
      estr1=0.0d0
      estr2=wwww
      call intgau(f4,u4)
      Cv0(m)=(((u3*u21*dexp(2.0d0*z))/qd)-u4*dexp(z))
      Cv(m)=(Cv0(m)+x(m)*dw(m)*(qd+(2.0d0/3.0d0)*&
      (z**(1.5d0))*(2.5d0-1.4d0*z)+((4.0d0/9.0d0)*(z**(4.0d0))*&
      (0.4d0*z-2.0d0)/qd)+ ((16.d0/135.d0)*(z**(6.5d0))/(qd*qd))))
      gw=dexp(z)*erf(zz)*(-z+(dsqrt(2.d0)-1.d0)/2.d0)-&
      (z+dsqrt(2.d0*z))/dsqrt(pai)
      hw=(zz/4.d0+dsqrt(pai)*gw/4.d0+(z**(1.5d0))/2.d0)/&
      ( -zz/2.d0+dexp(z)*dsqrt(pai)*erf(zz)/4.d0-(z**(1.5d0))/3.d0)
      Cvqt(m)=(Cv0(m)+0.5d0*(1.5d0-C+hw)*(-C-z)*xnhat(m))
   ENDDO
   wwww=w0
   estr2=w0
   qd=0.0d0
   call intgau(fcn1,qd)
   xnhat0=qd
   rho0hat=(4.d0*dsqrt(2.d0)*PAI)*qd
   xMcap=((9.d0/(4.d0*PAI))**(1.5d0))*(1.d0/dsqrt(rho0hat))*xmu

   qd=0.0d0
   xnhat(nfn)=0.0d0
   u1=0.0d0
   u21=0.0d0
   u23=0.0d0
   u3=0.0d0
   u4=0.0d0
   Cv0(nfn)=0.0d0
   Cv(nfn)=0.0d0
   Cvqt(nfn)=0.0d0

   sum0=0.0d0
   sum1=0.0d0
   sum2=0.0d0
   sum3=0.0d0
   b0=0.0d0
   b1=0.0d0
   b2=0.0d0
   b3=0.0d0
   h=0.0d0
   DO k=1,nfn-1
      b0=(Cv0(k+1)*csi(k+1)*csi(k+1)+Cv0(k)*csi(k)*csi(k))
      b1=(Cv(k+1)*csi(k+1)*csi(k+1)+Cv(k)*csi(k)*csi(k))
      b2=(Cvqt(k+1)*csi(k+1)*csi(k+1)+Cvqt(k)*csi(k)*csi(k))
      b3=(xnhat(k+1)*csi(k+1)*csi(k+1)+xnhat(k)*csi(k)*csi(k))
      h=csi(k+1)-csi(k)
      sum0=sum0+b0*h*0.5d0
      sum1=sum1+b1*h*0.5d0
      sum2=sum2+b2*h*0.5d0
      sum3=sum3+b3*h*0.5d0
   ENDDO
   Ctot0=0.0d0
   Ctot1=0.0d0
   Ctot2=0.0d0
   Ctot0=sum0/sum3
   Ctot1=sum1/sum3
   Ctot2=sum2/sum3

   xMhat=sum3
   DO i=1,10
      x0Cv(i)=0.0d0
      xi0Cv(i)=0.0d0
   ENDDO

   ex=0.0d0
   ex2=0.0d0
   ex3=0.0d0
   DO m=1,nfn-1
      IF((Cv(m).ge.0.0d0).AND.(Cv(m+1).le.0.0d0)&
         .AND.(ex.eq.0.0d0))THEN
         x0Cv(1)=(x(m)+x(m+1))*0.5d0
         xi0Cv(1)=(csi(m)+csi(m+1))*0.5d0
         ex=1.0d0
      ELSEIF((Cv(m).le.0.0d0).AND.(Cv(m+1).ge.0.0d0)&
         .AND.(ex2.eq.0.0d0))THEN
         x0Cv(2)=(x(m)+x(m+1))*0.5d0
         xi0Cv(2)=(csi(m)+csi(m+1))*0.5d0
         ex2=1.0d0
      ELSEIF((Cv(m).ge.0.0d0).AND.(Cv(m+1).le.0.0d0)&
         .AND.(ex.eq.1.0d0).AND.(ex3.eq.0.0d0))THEN
            x0Cv(3)=(x(m)+x(m+1))*0.5d0
            xi0Cv(3)=(csi(m)+csi(m+1))*0.5d0
            ex3=1.0d0
      ENDIF
   ENDDO
 9999 CONTINUE

   RETURN
   END

   function funzerr(x)
   IMPLICIT REAL*8 (A-H,O-Z)
   funzerr= dexp(-(x**2.0d0))
   return
   end

   function fcnqd(x)
   IMPLICIT REAL*8 (A-H,O-Z)
   COMMON/PARFCN/WWWW
   fcnqd = (dexp(wwww-x))*(x**1.5d0)
   return
   end

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
   fu4=(x*dexp(-x)-wwww*dexp(-wwww))*dlog(dexp(-x)-dexp(-wwww))*&
   (x**(0.5d0))
   return
   end

   SUBROUTINE FCT (X,Y,DERY, FCN)
   IMPLICIT REAL*8 (A-H,O-Z)
   COMMON/PARFCN/WWWW
   COMMON/ESTREMI/ESTR1,ESTR2
   COMMON/VALINI/W0,RHO0
   COMMON/COSTANTI/PAI
   COMMON/FCTOUTP/RPP
   DIMENSION Y(2),DERY(2)
   EXTERNAL :: fcn
   wwww = y(1)
   estr2= y(1)

   if(wwww.gt.0.d0) go to 1
   rpp = 0.d0
   go to 2
1      continue
   call intgau(fcn,qd)
   rpp = qd/rho0
2      continue

   dery(1) = y(2)
   dery(2) = -(2.d0/x)*y(2)-9.d0*rpp

   return
   end

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

   subroutine INTGAU(fcn,area)
	IMPLICIT REAL*8 (A-H,O-Z)
   real*8, external :: fcn
   common/gauss/indgau
   goto (1,2,3,4) ,indgau
 1	call GAUS20 (fcn,area)
   return
 2	call GAUS40 (fcn,area)
   return
 3	call GAUS80 (fcn,area)
   return
4 	call GAUS96 (fcn,area)
   return

	end
