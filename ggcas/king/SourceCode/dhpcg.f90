SUBROUTINE DHPCG(PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX,FCN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PRMT(5),Y(2),DERY(2),AUX(16,2)
      EXTERNAL fcn
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
      CALL FCT(X,Y,DERY,fcn)
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
      CALL FCT(X,Y,DERY,fcn)
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
      CALL FCT(X,Y,DERY,fcn)
      DO 20 I=1,NDIM
         AUX(3,I)=Y(I)
         AUX(10,I)=DERY(I)
 20   CONTINUE
      N=3
      ISW=4
      GOTO 100
 21   N=1
      X=X+H
      CALL FCT(X,Y,DERY,fcn)
      X=PRMT(1)
      DO 22 I=1,NDIM
	AUX(11,I)=DERY(I)
	Y(I)=AUX(1,I)+H*(.375D0*AUX(8,I)+.7916666666666667D0*AUX(9,I)&
      -.20833333333333333D0*AUX(10,I)+.041666666666666667D0*DERY(I))
 22   CONTINUE
 23   X=X+H
      N=N+1
      CALL FCT(X,Y,DERY,fcn)
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
      CALL FCT(Z,Y,DERY,fcn)
      DO 102 I=1,NDIM
         Z=H*DERY(I)
         AUX(6,I)=Z
      Y(I)=AUX(N,I)+.29697760924775360D0*AUX(5,I)+.15875964497103583D0*Z
 102  CONTINUE
	Z=X+.45573725421878943D0*H
	CALL FCT(Z,Y,DERY,fcn)
	DO 103 I=1,NDIM
	Z=H*DERY(I)
	AUX(7,I)=Z
	Y(I)=AUX(N,I)+.21810038822592047D0*AUX(5,I)-3.0509651486929308D0*&
            AUX(6,I)+3.8328647604670103D0*Z
 103  CONTINUE
      Z=X+H
      CALL FCT(Z,Y,DERY,fcn)
      DO 104 I=1,NDIM
       Y(I)=AUX(N,I)+.17476028226269037D0*AUX(5,I)-.55148066287873294D0*&
	   AUX(6,I)+1.2055355993965235D0*AUX(7,I)+.17118478121951903D0*&
	   H*DERY(I)
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
	DELT=AUX(N-4,I)+1.3333333333333333D0*H*(AUX(N+6,I)+AUX(N+6,I)-&
	AUX(N+5,I)+AUX(N+4,I)+AUX(N+4,I))
	Y(I)=DELT-.9256198347107438D0*AUX(16,I)
        AUX(16,I)=DELT
 207    CONTINUE
	CALL FCT(X,Y,DERY,fcn)
	DO 208 I=1,NDIM
	DELT=.125D0*(9.D0*AUX(N-1,I)-AUX(N-3,I)+3.D0*H*(DERY(I)+AUX(N+6,I)&
	+AUX(N+6,I)-AUX(N+5,I)))
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
        CALL FCT(X,Y,DERY,fcn)
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
	AUX(16,I)=8.962962962962963D0*(Y(I)-AUX(N-3,I))&
	-3.3611111111111111D0*H*(DERY(I)+DELT+AUX(N+4,I))
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
	Y(I)=.390625D-2*(8.D1*AUX(N-1,I)+135.D0*AUX(N-2,I)+4.D1*AUX(N-3,I)&
	+AUX(N-4,I))-.1171875D0*(AUX(N+6,I)-6.D0*AUX(N+5,I)-AUX(N+4,I))*H
	AUX(N-4,I)=.390625D-2*(12.D0*AUX(N-1,I)+135.D0*AUX(N-2,I)+&
	108.D0*AUX(N-3,I)+AUX(N-4,I))-.0234375D0*(AUX(N+6,I)+&
	18.D0*AUX(N+5,I)-9.D0*AUX(N+4,I))*H
	AUX(N-3,I)=AUX(N-2,I)
        AUX(N+4,I)=AUX(N+5,I)
 224    CONTINUE
	X=X-H
	DELT=X-(H+H)
	CALL FCT(DELT,Y,DERY,fcn)
	DO 225 I=1,NDIM
	AUX(N-2,I)=Y(I)
	AUX(N+5,I)=DERY(I)
        Y(I)=AUX(N-4,I)
 225    CONTINUE
	DELT=DELT-(H+H)
	CALL FCT(DELT,Y,DERY,fcn)
	DO 226 I=1,NDIM
	DELT=AUX(N+5,I)+AUX(N+4,I)
	DELT=DELT+DELT+DELT
	AUX(16,I)=8.962962962962963D0*(AUX(N-1,I)-Y(I))&
	-3.3611111111111111D0*H*(AUX(N+6,I)+DELT+DERY(I))
        AUX(N+3,I)=DERY(I)
 226    CONTINUE
	GOTO 206
	END
