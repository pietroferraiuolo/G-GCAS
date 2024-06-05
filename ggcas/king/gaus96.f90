	SUBROUTINE GAUS96 (FCN,AREA)
!
!       INTEGRA TRA A E B LA FUNZIONE FCN(X) CON IL METODO
!       DI GAUSS-LEGENDRE. SONO USATI 96 PUNTI.
!       I PARAMETRI DA CUI FCN(X) DIPENDE SONO PASSATI NELLA FUNCTION
!       FCN(X) CON IL COMMON/PARFCN/
!
!
	IMPLICIT REAL*8 (A-H,O-Z)
!    real*8 :: fcn, area
    COMMON/ESTREMI/A,B
	DIMENSION X(48),W(48)
!
	DATA X /.016276744849602969579, .048812985136049731112,&
            .081297495464425558994, .113695850110665920911,&
            .145973714654896941989, .178096882367618602759,&
            .210031310460567203603, .241743156163840012328,&
            .273198812591049141487, .304364944354496353024,&
            .335208522892625422616, .365696861472313635031,&
            .395797649828908603285, .425478988407300545365,&
            .454709422167743008636, .483457973920596359768,&
            .511694177154667673586, .539388108324357436227,&
            .566510418561397168404, .593032364777572080684,&
            .618925840125468570386, .644163403784967106798,&
            .668718310043916153953, .692564536642171561344,&
            .715676812348967626225, .738030643744400132851,&
            .759602341176647498703, .780369043867433217604,&
            .800308744139140817229, .819400310737931675539,&
            .837623511228187121494, .854959033434601455463,&
            .871388505909296502874, .886894517402420416057,&
            .901460635315852341319, .915071423120898074206,&
            .927712456722308690965, .939370339752755216932,&
            .950032717784437635756, .959688291448742539300,&
            .968326828463264212174, .975939174585136466453,&
            .982517263563014677447, .988054126329623799481,&
            .992543900323762624572, .995981842987209290650,&
            .998364375863181677724, .999689503883230766828/
!
    DATA W /.032550614492363166242, .032516118713868835987,&
            .032447163714064269364, .032343822568575928429,&
            .032206204794030250669, .032034456231992663218,&
            .031828758894411006535, .031589330770727168558,&
            .031316425596861355813, .031010332586313837423,&
            .030671376123669149014, .030299915420827593794,&
            .029896344136328385984, .029461089958167905970,&
            .028994614150555236543, .028497411065085385646, &
            .027970007616848334440, .027412962726029242823,&
            .026826866725591762198, .026212340735672413913,&
            .025570036005349361499, .024900633222483610288,&
            .024204841792364691282, .023483399085926219842,&
            .022737069658329374001, .021966644438744349195,&
            .021172939892191298988, .020356797154333324595,&
            .019519081140145022410, .018660679627411467385,&
            .017782502316045260838, .016885479864245172450,&
            .015970562902562291381, .015038721026994938006,&
            .014090941772314860916, .013128229566961572637,&
            .012151604671088319635, .011162102099838498591,&
            .010160770535008415758, .009148671230783386633,&
            .008126876925698759217, .007096470791153865269,&
            .006058545504235961683, .005014202742927517693,&
            .003964554338444686674, .002910731817934946408,&
            .001853960788946921732, .000796792065552012429/
!
	AREA=0.D0
	DO 1 I=1,48
	XI = (B-A)/2.D0*X(I)+(B+A)/2.0D0
        AREA = AREA + W(I)*FCN(XI)
	XI =-(B-A)/2.D0*X(I)+(B+A)/2.0D0
        AREA = AREA + W(I)*FCN(XI)
 1	CONTINUE
    AREA = (B-A)/2.D0*AREA 
    RETURN
    END

