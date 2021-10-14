C************************************************
C     PROGRAM PLATELET DEPOSITION AND ACTIVATION
C     Plasminogenesis-Included
C     CODED BY YINGMING YI, Modded by J.Jones
C     2021.10
C************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NX=200,NY=50)
      COMMON /BLK1/ RP(-1:NX+1,-1:NY+1),AP(-1:NX+1,-1:NY+1),
     &              APR(-1:NX+1,-1:NY+1),APS(-1:NX+1,-1:NY+1),
     &              PT(-1:NX+1,-1:NY+1),T(-1:NX+1,-1:NY+1),
     &              AT(-1:NX+1,-1:NY+1),FG(-1:NX+1,-1:NY+1),
     &              FB(-1:NX+1,-1:NY+1),TPA(-1:NX+1,-1:NY+1),
     &              PLS(-1:NX+1,-1:NY+1),PLA(-1:NX+1,-1:NY+1),
     &              R_L2AP(-1:NX+1,-1:NY+1)
      COMMON /BLK2/ U(0:NX,0:NY),V(0:NX,0:NY)
      COMMON /BLK3/ RPT(0:NX,0:NY),APT(0:NX,0:NY),
     &              APRT(0:NX,0:NY),APST(0:NX,0:NY) ,PTT(0:NX,0:NY),
     &              TT(0:NX,0:NY) ,ATT(0:NX,0:NY),
     &              FGT(0:NX,0:NY) ,FBT(0:NX,0:NY),
     &              TPAT(0:NX,0:NY),PLST(0:NX,0:NY),PLAT(0:NX,0:NY),
     &              R_L2APT(0:NX,0:NY)
      COMMON /BLK4/ XP(0:NX,0:NY), YP(0:NX,0:NY)
      COMMON /BLK5/ RJ(0:NX,0:NY), ALPHA(0:NX,0:NY),
     &              BETA(0:NX,0:NY), GAMMA(0:NX,0:NY)
      COMMON /BLK6/ ARJYTE(0:NX,0:NY), GRJYTN(0:NX,0:NY),
     &              DLX(0:NX,0:NY), DLY(0:NX,0:NY), YWL(0:NX,0:NY)
      COMMON /BLK7/ UC(0:NX,0:NY), VC(0:NX,0:NY)
      COMMON /BLK8/ XXI(0:NX,0:NY), YXI(0:NX,0:NY),
     &              XET(0:NX,0:NY), YET(0:NX,0:NY)
      COMMON /BLK9/ UP(0:NX,0:NY),VP(0:NX,0:NY)
      COMMON /BLK10/ OMEGA(0:NX,0:NY),RNUT(0:NX,0:NY),
     &              NU(0:NX,0:NY) 
C
      REAL*8    JR(0:NX),JA(0:NX),JPR(0:NX),JPS(0:NX),
     &              JPT(0:NX),JT(0:NX),DYB(0:NX),S(0:NX),
     &              M(0:NX),MTO(0:NX),MAS(0:NX),
     &              MR(0:NX), MAT(0:NX),DS(0:NX),
     &              BGAMMA(0:NX,0:NY),TAUD(0:NX,0:NY),
     &              DSS(0:NX,0:NY),RRJ2(0:NX,0:NY),
     &              KPA(0:NX,0:NY),F(0:NX,0:NY),OMEGA1(0:NX,0:NY),
     &              DT,ALLT,TACT,Drp,Dap,Dapr,Daps,Dpt,Dth,Dat,
     &              K1_T,KAT,KT,ALPHACOF,BETACOF,FEIAT,FEIRT,TIM,
     &              K1_J,Lambdaj,Spj,THETA,H,KRS,KAS,KAA,TAUDD,Re,
     &              UM, DSin(0:NY),TS,N,D1,
     &              Dfg,Dfb,K1,K1_m,JFG(0:NX),JFB(0:NX),
     &              KDEP,RHOS,DELTAPLT,DELTAFB,
C       ---Plasminogenesis-related constants
     &              Dtpa,Dpls,Dpla,Dl2ap,H_PLA,K_TPA,
     &              H1,H1_m

      CHARACTER NAME1*25,NAME2*25,NAME3*25,NAME4*25
      CHARACTER NAME5*25,NAME6*25,NAME7*25,NAME8*25
      CHARACTER NAME9*25,NAME10*25,NAME11*25
      CHARACTER NAME12*25,NAME13*25,NAME14*25
      CHARACTER NAME15*25,NAME16*25,NAME17*25
      CHARACTER NAME18*25,NAME19*25,NAME20*25
      CHARACTER NAME21*25,NAME22*25,NAME23*25,NAME24*25
      CHARACTER NAME25*25,NAME26*25,NAME27*25,NAME28*25
      CHARACTER NAME29*25,NAME30*25,NAME31*25,NAME32*25
      CHARACTER NAME33*25,NAME34*25,NAME35*25,NAME36*25

C-------------------------------------------------------------------------
C
      PARAMETER (D=2.2,Pi=3.14159265359)
      OPEN(20,FILE='../GRID_D/AB',FORM='UNFORMATTED')
      OPEN(21,FILE='../RESULT/VEL.DAT',FORM='UNFORMATTED')
      OPEN(24,FILE='../A-TOCYUU/TAUD.DAT',FORM='UNFORMATTED')
      OPEN(25,FILE='../RESULT/RNUT.DAT',FORM='UNFORMATTED', 
     &     STATUS='UNKNOWN') 
 
       READ(20) XP, YP, XXI, YXI, XET, YET,
     &         RJ, ALPHA, BETA, GAMMA,
     &         ARJYTE, GRJYTN, DLX, DLY,YWL
      READ(21) U, V
      READ(24) TAUD
      READ(25) RNUT
      CLOSE(20)
      CLOSE(21)
      CLOSE(24)
      CLOSE(25)

c      WRITE(6,*) 'PLEASE Input Re Number'
c      READ(5,*) Re
      Re = 5300

        UM = Re/3180.0D0*13.159973D0

      WRITE(6,*) 'RE=',RE,'UM=', UM

c      WRITE(6,*) 'PLEASE Choose Re Number:3180,4240,5300,6360,7420'
c      READ(5,*) Re
      Re = 5300

c      WRITE(6,*) 'PLEASE Choose AB:0 or CB 1'
c      READ(5,*) IFLAG
      IFLAG = 0

      IF(RE.EQ.3180)THEN
        UM = 13.159973
      ENDIF
      IF(RE.EQ.4240)THEN
        UM = 17.54663
      ENDIF
      IF(RE.EQ.5300)THEN
        UM = 21.9
      ENDIF
      IF(RE.EQ.6360)THEN
        UM = 26.319945
      ENDIF
      IF(RE.EQ.7420)THEN
        UM = 30.706603
      ENDIF
      WRITE(6,*) 'RE=',RE,'UM=', UM

      TRP = 0D0
      TAP = 0D0
      TAPP = 0D0
      TAPR = 0D0
      TAPS = 0D0
      TPT = 0D0
      TTRB = 0D0
      TAT = 0D0
      TDEP = 0D0
      TDEPRP = 0D0
      TDEPAP = 0D0
      TDEPAPP = 0D0
      TFG = 0D0
      TFB = 0D0
      TTPA=0D0
      TPLS=0D0
      TPLA=0D0
      TL2AP=0D0


      DO 25 I=0,NX
       DO 26 J=0,NY
       RX = 1D0
       RY = 1D0

       IF(I.EQ.0.OR.I.EQ.NX) THEN
          RX = 0.5D0
       ENDIF
       IF(J.EQ.0.OR.J.EQ.NY) THEN
          RY = 0.5D0
       ENDIF

      DSS(I,J) = RJ(I,J)*RX*RY*D*D
       
  26   ENDDO
  25  ENDDO

      DO 37 I=0,NX 
       DO 38 J=0,NY 

      RNUT(I,J) = RNUT(I,J)*UM*D 
       
  38   ENDDO 
  37  ENDDO 

       DO 950 J=1,NY-1
 
      DSin(J) = (YP(0,J+1)-YP(0,J-1))*0.5*2*Pi*
     &        YP(0,J) *D*D
       
  950   ENDDO

      DSin(0) = (YP(0,1)-YP(0,0))*0.5*2*Pi*
     &        YP(0,0) *D*D
      DSin(NY) = (YP(0,NY)-YP(0,NY-1))*0.5*2*Pi*
     &        YP(0,NY) *D*D

      KTT=0
      ALLT=12.0D0
      H=0.417
      KRS=0.0037
      KAS=0.0046
      KAA=0.0046	
      TACT=1D0
      THETA=0D0

c      WRITE(6,*) 'PLEASE INPUT Threshold Shear Rate Or 0 means no limit'
c      READ(5,*) STH
       STH = 500

        OPEN(70,FILE='Threshold Shear Rate Number',STATUS='UNKNOWN')
        WRITE(70,*) STH

      DO 400 I=0,NX
       DYB(I) = D*((YP(I,NY) - YP(I,NY-1))**2.0+(XP(I,NY)
     &               - XP(I,NY-1))**2.0)**0.5
  400  ENDDO

      DO 401 I=0,NX
       DO 402 J=0,NY

       TAUD(I,J) = TAUD(I,J)*UM/D

       IF(TAUD(I,J).GE.STH) THEN
          F(I,J) = 1D0
       ELSE
          F(I,J) = 0D0
        ENDIF
 
  402   ENDDO
  401  ENDDO

C----------PARAMETER SET ----------------------

      Drp=1.58D-9
      Dap=1.58D-9
      Dapr=2.57D-6
      Daps=2.14D-6
      Dpt=3.32D-7
      Dth=4.16D-7
      Dat=3.49D-7
      Dfg=3.10D-7
      Dfb=2.47D-7
      Dtpa=5.28D-7
      Dpls=4.81D-7
      Dpla=4.93D-7
      Dl2ap=5.25D-7
      
      K1_T=13.333
      KAT=0.100
      KT=3.50D-2
      ALPHACOF=1.0
      BETACOF=9.11D-3
      FEIAT=3.69D-9
      FEIRT=6.50D-10
      K1_J=0.0161
      Lambdaj=2.4D-8
      Spj=9.5D-12
C--Reaction rate K1(/s) and Michaelis-Menten constant K1_m(μM)
      K1=59
      K1_m=3.16
C--Reaction rates for Plasminogenesis (Anand,2003)
      H1=0.1161
      H1_m=2.90
      K_TPA=0.2
C--H_PLA = 160*10^4 (M^-1 s^-1)
C        = 1.6*10^6 (M^-1 s^-1)
C        = 1.6 (10^-6 * M s)^-1
C        = 1.6 (μM^-1 s^-1)
      H_PLA=1.6

      RHOS=1050
      DELTAPLT=2.0D-12
      DELTAFB=66D-5

C---Note that the actual time-step is set to 
C-- 1/10 of the laminar version      
      DT=1.0D-6
    


c      WRITE(6,*) 'PLEASE Confirm whether first cal. YES:0 ,NO:1'
c      READ(5,*) IEG
       IEG = 0


       IF(IEG.EQ.0) THEN

      DO 10 I=-1,NX+1
       DO 11 J=-1,NY+1
       RP(I,J)=0.0
       AP(I,J)=0.0
       APR(I,J)=0.0
       APS(I,J)=0.0
       PT(I,J)=0.0
       T(I,J)=0.0
       AT(I,J)=0.0
       FG(I,J)=0.0
       FB(I,J)=0.0
       TPA(I,J)=0.0
       PLS(I,J)=0.0
       PLA(I,J)=0.0
       R_L2AP(I,J)=0.0
  11   ENDDO
  10  ENDDO
C
C--Initial concentration(micro-M)
C--Reference: Sorensen(1999) and Anand(2003) for FG, FB
      DO 114 I=5,10
       DO 115 J=-1,NY+1
       RP(I,J)=2.0D8
       AP(I,J)=0.0
       APR(I,J)=0.0
       APS(I,J)=0.0
       PT(I,J)=1.1
       T(I,J)=0.0
       AT(I,J)=2.844
       FG(I,J)=7.0
       FB(I,J)=0.0
       TPA(I,J)=8.0D-5
       PLS(I,J)=2.18
       PLA(I,J)=0.0
       R_L2AP(I,J)=0.105
  115   ENDDO
  114  ENDDO

      DO 31 I=5,10
       RP(I,NY+1)=RP(I,NY)
       AP(I,NY+1)=AP(I,NY)
       APR(I,NY+1)=APR(I,NY)
       APS(I,NY+1)=APS(I,NY)
       PT(I,NY+1)=PT(I,NY)
       T(I,NY+1)=T(I,NY)
       AT(I,NY+1)=AT(I,NY)
       FG(I,NY+1)=FG(I,NY)
       FB(I,NY+1)=FB(I,NY)
       TPA(I,NY+1)=TPA(I,NY)
       PLS(I,NY+1)=PLS(I,NY)
       PLA(I,NY+1)=PLA(I,NY)
       R_L2AP(I,NY+1)=R_L2AP(I,NY)
  31  ENDDO

C--Surface-flux model
      TIM=0D0
      KDEP=0D0

      DO 12 I=0, NX
         JR(I) = 0.0
         JA(I) = 0.0
         JPR(I) = 0.0
         JPS(I) = 0.0
         JPT(I) = 0.0
         JT(I) = 0.0
         S(I) = 1.0
         M(I) = 0.0
         MAS(I) = 0.0
         MR(I) = 0.0
         MAT(I) = 0.0
C--Added reactivity wall flux for FG, FB         
         JFG(I)=0.0
         JFB(I)=0.0
  12  CONTINUE

       ELSE
       
      WRITE(6,*) 'PLEASE Input KTT'
      READ(5,*) KTT

      TIM=KTT*0.00001

       K13=KTT
       K12=K13/1000000000
       K11=K13-K12*1000000000
       K10=K11/100000000
       KK9=K11-K10*100000000
       KK8=KK9/10000000
       KK7=KK9-KK8*10000000
       KK6=KK7/1000000
       KK5=KK7-KK6*1000000
       KK4=KK5/100000
       KK3=KK5-KK4*100000
       KK2=KK3/10000
       KK1=KK3-KK2*10000
       KK0=KK1/1000
       KL1=KK1-KK0*1000
       KL2=KL1/100
       KL3=KL1-KL2*100
       KL4=KL3/10
       KL5=KL3-KL4*10
       KL6=KL5/1

       NAME1='RP'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME2='AP'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME3='APR'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME4='APS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME5='PT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME6='T'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME7='AT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME23='FG'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME24='FB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'     


       NAME11='MTO'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME12='M'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME13='MAS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME14='MR'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME15='MAT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

      OPEN(10001,FILE='../fei 100s/'//NAME1,STATUS='UNKNOWN')
      OPEN(10002,FILE='../fei 100s/'//NAME2,STATUS='UNKNOWN')
      OPEN(10003,FILE='../fei 100s/'//NAME3,STATUS='UNKNOWN')
      OPEN(10004,FILE='../fei 100s/'//NAME4,STATUS='UNKNOWN')
      OPEN(10005,FILE='../fei 100s/'//NAME5,STATUS='UNKNOWN')
      OPEN(10006,FILE='../fei 100s/'//NAME6,STATUS='UNKNOWN')
      OPEN(10007,FILE='../fei 100s/'//NAME7,STATUS='UNKNOWN')
      OPEN(10011,FILE='../fei 100s/'//NAME11,STATUS='UNKNOWN')
      OPEN(10012,FILE='../fei 100s/'//NAME12,STATUS='UNKNOWN')
      OPEN(10013,FILE='../fei 100s/'//NAME13,STATUS='UNKNOWN')
      OPEN(10014,FILE='../fei 100s/'//NAME14,STATUS='UNKNOWN')
      OPEN(10015,FILE='../fei 100s/'//NAME15,STATUS='UNKNOWN')
      OPEN(10023,FILE='../fei 100s/'//NAME23,STATUS='UNKNOWN')
      OPEN(10024,FILE='../fei 100s/'//NAME24,STATUS='UNKNOWN')


      READ(10001,*)
      READ(10002,*)
      READ(10003,*)
      READ(10004,*)
      READ(10005,*)
      READ(10006,*)
      READ(10007,*)
      READ(10023,*)
      READ(10024,*)

      DO 420 J = 0, NY
        DO 420 I = 0, NX

      READ(10001,*) RP(I,J),RP(I,J),RP(I,J)
      READ(10002,*) AP(I,J),AP(I,J),AP(I,J)
      READ(10003,*) APR(I,J),APR(I,J),APR(I,J)
      READ(10004,*) APS(I,J),APS(I,J),APS(I,J)
      READ(10005,*) PT(I,J),PT(I,J),PT(I,J)
      READ(10006,*) T(I,J),T(I,J),T(I,J)
      READ(10007,*) AT(I,J),AT(I,J),AT(I,J)
      READ(10023,*) FG(I,J),FG(I,J),FG(I,J)
  420 READ(10024,*) FB(I,J),FB(I,J),FB(I,J)


       DO 450 J=-1,NY+1

       RP(-1,J)=RP(0,J)
       AP(-1,J)=AP(0,J)
       APR(-1,J)=APR(0,J)
       APS(-1,J)=APS(0,J)
       PT(-1,J)=PT(0,J)
       T(-1,J)=T(0,J)
       AT(-1,J)=AT(0,J)
       FG(-1,J)=FG(0,J)
       FB(-1,J)=FB(0,J)

        RP(NX+1,J)=RP(NX,J)
        AP(NX+1,J)=AP(NX,J)
        APR(NX+1,J)=APR(NX,J)
        APS(NX+1,J)=APS(NX,J)
        PT(NX+1,J)=PT(NX,J)
        T(NX+1,J)=T(NX,J)
        AT(NX+1,J)=AT(NX,J)
        FG(NX+1,J)=FG(NX,J)
        FB(NX+1,J)=FB(NX,J)

  450   ENDDO

      DO 460 I = -1, NX+1
        RP(I,NY+1)=RP(I,NY)
        AP(I,NY+1)=AP(I,NY)
        APR(I,NY+1)=APR(I,NY)
        APS(I,NY+1)=APS(I,NY)
        PT(I,NY+1)=PT(I,NY)
        T(I,NY+1)=T(I,NY)
        AT(I,NY+1)=AT(I,NY)
        FG(I,NY+1)=FG(I,NY)
        FB(I,NY+1)=FB(I,NY)

        RP(I,-1)=RP(I,0)
        AP(I,-1)=AP(I,0)
        APR(I,-1)=APR(I,0)
        APS(I,-1)=APS(I,0)
        PT(I,-1)=PT(I,0)
        T(I,-1)=T(I,0)
        AT(I,-1)=AT(I,0)
        FG(I,-1)=FG(I,0)
        FB(I,-1)=FB(I,0)
  460  CONTINUE



      READ(10011,*) MTO
      READ(10012,*) M
      READ(10013,*) MAS
      READ(10014,*) MR
      READ(10015,*) MAT

       DO 410 I = 0, NX    
           S(I) = 1.0- M(I)/7.0D6
  410  CONTINUE

C--Reactive wall flux  
      DO 360 I=0, NX
          JR(I) = S(I)*KRS*RP(I,NY)
          JA(I) = (S(I)*KAS+(MAS(I)/7.0D6)*KAA)*AP(I,NY)
          JPR(I) = Lambdaj*(THETA*S(I)*KRS*RP(I,NY))
          JPS(I) = MAT(I)*Spj
          JPT(I) = BETACOF*(MAT(I)*(3.22D-8)+MR(I)*FEIRT)*PT(I,NY)
          JT(I) = (MAT(I)*(3.22D-8)+MR(I)*FEIRT)*PT(I,NY)
          JFB(I) = ((KDEP*FB(I,NY))/(RHOS-DELTAFB*KDEP*FB(I,NY)))
     &            *(DELTAPLT*JR(I)+DELTAPLT*JA(I))
C--Fibrin
C          JFB = KDEP * FB(I,NY) * (H(I-1,NY)-H(I,NY))/ DT         
  360 CONTINUE

       ENDIF



      DO 213 I = 1, NX-1
       
           DS(I) = SQRT(((XP(I-1,NY)-XP(I+1,NY))*0.5)**2
     &             +((YP(I-1,NY)-YP(I+1,NY))*0.5)**2)*D

  213  CONTINUE

           DS(0) = SQRT(((XP(0,NY)-XP(1,NY))*0.5)**2
     &             +((YP(0,NY)-YP(1,NY))*0.5)**2)*D
           DS(NX) = SQRT(((XP(NX,NY)-XP(NX-1,NY))*0.5)**2
     &             +((YP(NX,NY)-YP(NX-1,NY))*0.5)**2)*D

      DO 30 I=0,NX
       DO 40 J=0,NY
         
          RRJ2(I,J) = 1.0/(RJ(I,J)*RJ(I,J))
  
  40   ENDDO
  30  ENDDO



C----------TIME CHANGE---------------------------

C-----------VELOCITY CALCULATION-----------------
C
      DO 210 J = 0, NY
        DO 200 I = 0, NX
          UP(I,J) = (  YET(I,J)*U(I,J) - XET(I,J)*V(I,J) )*YP(I,J)
          VP(I,J) = (- YXI(I,J)*U(I,J) + XXI(I,J)*V(I,J) )*YP(I,J) 
  200   CONTINUE
  210 CONTINUE
C
      DO 222 J = 0, NY-1
        DO 223 I = 0, NX-1
          UC(I,J) = 0.5D0*( UP(I+1,J) + UP(I,J) )
          VC(I,J) = 0.5D0*( VP(I,J+1) + VP(I,J) )*MIN(1,J,NY-1-J)
  223   CONTINUE
  222 ENDDO

        DO 224 I = 0, NX-1
          UC(I,NY) = 0.0
          VC(I,NY) = 0.0
  224   CONTINUE

  300  DO 430 I=0, NX
        DO 440 J=0, NY

       IF(IFLAG.EQ.0) THEN

       IF(XP(I,J)*D.GE.-2.0.AND.XP(I,J)*D.LE.17.0) THEN

       OMEGA(I,J) = (APR(I,J)/2.0 
     &       + APS(I,J)/0.6 + T(I,J)/0.1)*F(I,J)
       KPA(I,J) = OMEGA(I,J)/TACT

       ELSE
          KPA(I,J) = 0D0
        ENDIF

       ELSE

       IF(XP(I,J)*D.GE.-5.0.AND.XP(I,J)*D.LE.14.0) THEN

       OMEGA(I,J) = (APR(I,J)/2.0 
     &       + APS(I,J)/0.6 + T(I,J)/0.1)*F(I,J)
       KPA(I,J) = OMEGA(I,J)/TACT

       ELSE
          KPA(I,J) = 0D0
        ENDIF

        ENDIF

  440   CONTINUE
  430 CONTINUE

C--------------CONCENTRATION CALCULATION---------------
C-----------------Resting PLT CALCULATION-------------
C
	    D1=Drp 

      DO 610 J = 1, NY-1
        DO 600 I = 1, NX-1

CC -- CONVECTION  --

          CALL CNVUP(RP, I, J, CNVRP)
CC -- DIFFUSION  --
          CALL FDIFF(RP,RPXIXI,RPXIET,RPETXI,RPETET,I,J,D1)
C
          RPXI = 0.5D0*(RP(I+1,J) - RP(I-1,J))
          RPET = 0.5D0*(RP(I,J+1) - RP(I,J-1))
          DIFFRP =( ALPHA(I,J)*RPXIXI -   BETA(I,J)*(RPXIET + RPETXI)
     &             + GAMMA(I,J)*RPETET
     &     + (RNUT(I,J)+Drp)*( DLX(I,J)*( YXI(I,J)*RPET 
     &      - YET(I,J)*RPXI )
     &             + DLY(I,J)*( XET(I,J)*RPXI - XXI(I,J)*RPET ) )
     &           )*RRJ2(I,J)
C
          RPT(I,J) = RP(I,J) + DT*( - UM/D*CNVRP +  1.0/D/D*
     &           DIFFRP-KPA(I,J)*RP(I,J))


  600  CONTINUE
  610 ENDDO

        DO 500 I = 1, NX-1
          RPT(I,NY) = RPT(I,NY-1) - DT*2*JR(I)/DYB(I)
  500  CONTINUE

C
C-------------Activated PLT CALCULATION---------
      D1=Dap

      DO 611 J = 1, NY-1
      DO 601 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(AP, I, J, CNVAP)
CC -- DIFFUSION  --
          CALL FDIFF(AP,APXIXI,APXIET,APETXI,APETET,I,J,D1)
C
          APXI = 0.5D0*(AP(I+1,J) - AP(I-1,J))
          APET = 0.5D0*(AP(I,J+1) - AP(I,J-1))
          DIFFAP =( ALPHA(I,J)*APXIXI -   BETA(I,J)*(APXIET + APETXI)
     &             + GAMMA(I,J)*APETET
     &     + (RNUT(I,J)+Dap)*( DLX(I,J)*( YXI(I,J)*APET
     &     - YET(I,J)*APXI )
     &             + DLY(I,J)*( XET(I,J)*APXI - XXI(I,J)*APET ) )
     &           )*RRJ2(I,J)
C
          APT(I,J) = AP(I,J) + DT*( - UM/D*CNVAP + 1.0/D/D*
     &           DIFFAP +KPA(I,J)*RP(I,J))

  601   CONTINUE
  611 ENDDO

      DO 501 I = 1, NX-1
          APT(I,NY) = APT(I,NY-1) - DT*2*JA(I)/DYB(I)
  501   CONTINUE

C-------------PLT-released agonists CALCULATION----------
      D1=Dapr

      DO 612 J = 1, NY-1
      DO 602 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(APR, I, J, CNVAPR)
CC -- DIFFUSION  --
          CALL FDIFF(APR,APRXIXI,APRXIET,APRETXI,APRETET,I,J,D1)
C
          APRXI = 0.5D0*(APR(I+1,J) - APR(I-1,J))
          APRET = 0.5D0*(APR(I,J+1) - APR(I,J-1))
          DIFFAPR = ( ALPHA(I,J)*APRXIXI -   BETA(I,J)*(APRXIET 
     &             + APRETXI) + GAMMA(I,J)*APRETET
     &     + (RNUT(I,J)+Dapr)*( DLX(I,J)*( YXI(I,J)*APRET 
     &     - YET(I,J)*APRXI )
     &     + DLY(I,J)*( XET(I,J)*APRXI - XXI(I,J)*APRET ) ))*RRJ2(I,J)

          APRT(I,J) = APR(I,J) + DT*( -UM/D*CNVAPR + 1.0/D/D*
     &          DIFFAPR +Lambdaj*KPA(I,J)*RP(I,J)-K1_J*APR(I,J))
  602   CONTINUE
  612 ENDDO 

      DO 502 I = 1, NX-1
          APRT(I,NY) = APRT(I,NY-1) + DT*2*JPR(I)/DYB(I)
  502   CONTINUE

C----------------- PLT-synthesized agonist CALCULATION-------------
C
      D1=Daps

      DO 613 J = 1, NY-1
      DO 603 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(APS, I, J, CNVAPS)
CC -- DIFFUSION  --
          CALL FDIFF(APS,APSXIXI,APSXIET,APSETXI,APSETET,I,J,D1)
C
          APSXI = 0.5D0*(APS(I+1,J) - APS(I-1,J))
          APSET = 0.5D0*(APS(I,J+1) - APS(I,J-1))
          DIFFAPS =( ALPHA(I,J)*APSXIXI -   BETA(I,J)*(APSXIET 
     &             + APSETXI) + GAMMA(I,J)*APSETET
     &     + (RNUT(I,J)+Daps)*( DLX(I,J)*( YXI(I,J)*APSET 
     &     - YET(I,J)*APSXI )
     &             + DLY(I,J)*( XET(I,J)*APSXI - XXI(I,J)*APSET ) )
     &           )*RRJ2(I,J)
C
          APST(I,J) = APS(I,J) + DT*( - UM/D*CNVAPS + 1.0/D/D*
     &                    DIFFAPS +spj*AP(I,J)-k1_j*APS(I,J))       
  603  CONTINUE
  613 ENDDO

      DO 503 I = 1, NX-1
          APST(I,NY) = APST(I,NY-1) + DT*2*JPS(I)/DYB(I)
  503  CONTINUE


C
C----------------- Prothrombin  CALCULATION-------------
C
      D1=Dpt

      DO 614 J = 1, NY-1
      DO 604 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(PT, I, J, CNVPT)
CC -- DIFFUSION  --
          CALL FDIFF(PT,PTXIXI,PTXIET,PTETXI,PTETET,I,J,D1)
C
          PTXI = 0.5D0*(PT(I+1,J) - PT(I-1,J))
          PTET = 0.5D0*(PT(I,J+1) - PT(I,J-1))
          DIFFPT =( ALPHA(I,J)*PTXIXI -   BETA(I,J)*(PTXIET + PTETXI)
     &             + GAMMA(I,J)*PTETET
     &     + (RNUT(I,J)+Dpt)*( DLX(I,J)*( YXI(I,J)*PTET
     &     - YET(I,J)*PTXI )
     &             + DLY(I,J)*( XET(I,J)*PTXI - XXI(I,J)*PTET ) )
     &           )*RRJ2(I,J)
C
          PTT(I,J) = PT(I,J) + DT*( - UM/D*CNVPT + 1.0/D/D*DIFFPT 
     &                -BETACOF*PT(I,J)*(FEIAT*AP(I,J)+FEIRT*RP(I,J)))   
  604  CONTINUE
  614 ENDDO

      DO 504 I = 1, NX-1
          PTT(I,NY) = PTT(I,NY-1) - DT*2*JPT(I)/DYB(I)     
  504  CONTINUE

      DO 6160 I = 0, NX
        DO 6060 J = 0, NY
        BGAMMA(I,J) = K1_T*H*AT(I,J)/(ALPHACOF*KAT*KT
     &                  +ALPHACOF*KAT* BETACOF*T(I,J)              
     &                  +AT(I,J)*BETACOF*T(I,J))              

6060  CONTINUE
6160  ENDDO

C
C----------------- Thrombin  CALCULATION-------------
C
      D1=Dth

      DO 615 J = 1, NY-1
      DO 605 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(T, I, J, CNVT)
CC -- DIFFUSION  --
          CALL FDIFF(T,TXIXI,TXIET,TETXI,TETET,I,J,D1)
C
          TXI = 0.5D0*(T(I+1,J) - T(I-1,J))
          TET = 0.5D0*(T(I,J+1) - T(I,J-1))
          DIFFT =( ALPHA(I,J)*TXIXI -   BETA(I,J)*(TXIET + TETXI)
     &             + GAMMA(I,J)*TETET
     &     +  (RNUT(I,J)+Dth)*( DLX(I,J)*( YXI(I,J)*TET - YET(I,J)*TXI )
     &             + DLY(I,J)*( XET(I,J)*TXI - XXI(I,J)*TET ) )
     &           )*RRJ2(I,J)
C
          TT(I,J) = T(I,J) + DT*( - UM/D*CNVT + 1.0/D/D*DIFFT
     &        -BGAMMA(I,J)*T(I,J)+PT(I,J)*(FEIAT*AP(I,J)+ 
     &           FEIRT*RP(I,J)))  
  605  CONTINUE
  615 ENDDO

      DO 505 I = 1, NX-1
       TT(I,NY) = TT(I,NY-1) + DT*2*JT(I)/DYB(I)
  505  CONTINUE



C
C----------------- ATIII  CALCULATION-------------
C
      D1=Dat

      DO 606 I = 1, NX-1
        DO 616 J = 1, NY-1
CC -- CONVECTION  --
          CALL CNVUP(AT, I, J, CNVAT)
CC -- DIFFUSION  --
          CALL FDIFF(AT,ATXIXI,ATXIET,ATETXI,ATETET,I,J,D1)
C
          ATXI = 0.5D0*(AT(I+1,J) - AT(I-1,J))
          ATET = 0.5D0*(AT(I,J+1) - AT(I,J-1))
          DIFFAT =( ALPHA(I,J)*ATXIXI -   BETA(I,J)*(ATXIET + ATETXI)
     &             + GAMMA(I,J)*ATETET
     &     +  (RNUT(I,J)+Dat) *( DLX(I,J)*( YXI(I,J)*ATET 
     &     - YET(I,J)*ATXI )
     &             + DLY(I,J)*( XET(I,J)*ATXI - XXI(I,J)*ATET ) )
     &           )*RRJ2(I,J)
C
          ATT(I,J) = (AT(I,J) + DT*( - UM/D*CNVAT + 1.0/D/D*DIFFAT
     &                 -BGAMMA(I,J)*BETACOF*T(I,J)))   
  616  CONTINUE
  606 ENDDO
  
      DO 506 I = 1, NX-1
          ATT(I,NY) = ATT(I,NY-1) 
  506 CONTINUE

C-------------FIBRINOGEN CALCULATION---------
      D1=Dfg

      DO 617 J = 1, NY-1
      DO 607 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(FG, I, J, CNVFG)
CC -- DIFFUSION  --
          CALL FDIFF(FG,FGXIXI,FGXIET,FGETXI,FGETET,I,J,D1)
C
          FGXI = 0.5D0*(FG(I+1,J) - FG(I-1,J))
          FGET = 0.5D0*(FG(I,J+1) - FG(I,J-1))
          DIFFFG =( ALPHA(I,J)*FGXIXI -   BETA(I,J)*(FGXIET + FGETXI)
     &             + GAMMA(I,J)*FGETET
     &     +  (RNUT(I,J)+Dfg)*( DLX(I,J)*( YXI(I,J)*FGET 
     &      - YET(I,J)*FGXI )
     &             + DLY(I,J)*( XET(I,J)*FGXI - XXI(I,J)*FGET ) )
     &           )*RRJ2(I,J)
C
      FGT(I,J) = FG(I,J) + DT*( - UM/D*CNVFG + 1.0/D/D*
     &           DIFFFG - (K1*BETACOF*T(I,J)*FG(I,J))/(K1_M + FG(I,J)))

  607   CONTINUE
  617 ENDDO

C-- Reactive wall flux is unknown for FG
      DO 507 I = 1, NX-1
          FGT(I,NY) = FGT(I,NY-1)
  507   CONTINUE

C------------------------------------------------

C-------------FIBRIN CALCULATION--------- 
      D1=Dfb

      DO 618 J = 1, NY-1
      DO 608 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(FB, I, J, CNVFB)
CC -- DIFFUSION  --
          CALL FDIFF(FB,FBXIXI,FBXIET,FBETXI,FBETET,I,J,D1)
C
          FBXI = 0.5D0*(FB(I+1,J) - FB(I-1,J))
          FBET = 0.5D0*(FB(I,J+1) - FB(I,J-1))
          DIFFFB =( ALPHA(I,J)*FBXIXI -   BETA(I,J)*(FBXIET + FBETXI)
     &             + GAMMA(I,J)*FBETET
     &     +  (RNUT(I,J)+Dfb)*( DLX(I,J)*( YXI(I,J)*FBET 
     &    - YET(I,J)*FBXI )
     &             + DLY(I,J)*( XET(I,J)*FBXI - XXI(I,J)*FBET ) )
     &           )*RRJ2(I,J)
C
      FBT(I,J) = FB(I,J) + DT*( -UM/D*CNVFB + 1.0/D/D*
     &           DIFFFB + (K1*BETACOF*T(I,J)*FG(I,J))/(K1_M + FG(I,J))
     &                  - H1*PLA(I,J)*FB(I,J)/(H1_m+FB(I,J)))

  608   CONTINUE
  618 ENDDO

C-- Reactive wall flux for FB is assumed to be 0, as blood clot isn't simulated.
C-- reference may be in Strong(1987)
      DO 508 I = 1, NX-1
          FBT(I,NY) = FBT(I,NY-1) - DT*2*JFB(I)/DYB(I)
  508   ENDDO

C-------------TPA CALCULATION---------
      D1=Dtpa

      DO 619 J = 1, NY-1
      DO 609 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(TPA, I, J, CNVTPA)
CC -- DIFFUSION  --
          CALL FDIFF(TPA,TPAXIXI,TPAXIET,TPAETXI,TPAETET,I,J,D1)
C
          TPAXI = 0.5D0*(TPA(I+1,J) - TPA(I-1,J))
          TPAET = 0.5D0*(TPA(I,J+1) - TPA(I,J-1))
          DIFFTPA =( ALPHA(I,J)*TPAXIXI - BETA(I,J)*(TPAXIET+ TPAETXI)
     &             + GAMMA(I,J)*TPAETET
     &     +  (RNUT(I,J)+Dtpa)*( DLX(I,J)*( YXI(I,J)*TPAET 
     &      - YET(I,J)*TPAXI )
     &             + DLY(I,J)*( XET(I,J)*TPAXI - XXI(I,J)*TPAET ) )
     &           )*RRJ2(I,J)
C
      TPAT(I,J) = TPA(I,J) + DT*( - UM/D*CNVTPA + 1.0/D/D*
     &           DIFFTPA +0)

  609   CONTINUE
  619 ENDDO

C-- Reactive wall flux is unknown for TPA
      DO 509 I = 1, NX-1
          TPAT(I,NY) = TPAT(I,NY-1)
  509   CONTINUE

C-------------PLS CALCULATION---------
      D1=Dpls

      DO 630 J = 1, NY-1
      DO 620 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(PLS, I, J, CNVPLS)
CC -- DIFFUSION  --
          CALL FDIFF(PLS,PLSXIXI,PLSXIET,PLSETXI,PLSETET,I,J,D1)
C
          PLSXI = 0.5D0*(PLS(I+1,J) - PLS(I-1,J))
          PLSET = 0.5D0*(PLS(I,J+1) - PLS(I,J-1))
          DIFFPLS =( ALPHA(I,J)*PLSXIXI-BETA(I,J)*(PLSXIET + PLSETXI)
     &             + GAMMA(I,J)*PLSETET
     &     +  (RNUT(I,J)+Dpls)*( DLX(I,J)*( YXI(I,J)*PLSET 
     &      - YET(I,J)*PLSXI )
     &             + DLY(I,J)*( XET(I,J)*PLSXI - XXI(I,J)*PLSET ) )
     &           )*RRJ2(I,J)
C
      PLST(I,J) = PLS(I,J) + DT*( - UM/D*CNVPLS + 1.0/D/D*
     &           DIFFPLS -K_TPA*TPA(I,J))

  620   CONTINUE
  630 ENDDO

C-- Reactive wall flux is unknown for PLS
      DO 520 I = 1, NX-1
          PLST(I,NY) = PLST(I,NY-1)
  520   CONTINUE

C-------------PLA CALCULATION---------
      D1=Dpla

      DO 631 J = 1, NY-1
      DO 621 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(PLA, I, J, CNVPLA)
CC -- DIFFUSION  --
          CALL FDIFF(PLA,PLAXIXI,PLAXIET,PLAETXI,PLAETET,I,J,D1)
C
          PLAXI = 0.5D0*(PLA(I+1,J) - PLA(I-1,J))
          PLAET = 0.5D0*(PLA(I,J+1) - PLA(I,J-1))
          DIFFPLA =( ALPHA(I,J)*PLAXIXI-BETA(I,J)*(PLAXIET + PLAETXI)
     &             + GAMMA(I,J)*PLAETET
     &     +  (RNUT(I,J)+Dpla)*( DLX(I,J)*( YXI(I,J)*PLAET 
     &      - YET(I,J)*PLAXI )
     &             + DLY(I,J)*( XET(I,J)*PLAXI - XXI(I,J)*PLAET ) )
     &           )*RRJ2(I,J)
C
      PLAT(I,J) = PLA(I,J) + DT*( - UM/D*CNVPLA + 1.0/D/D*
     &        DIFFPLA +K_TPA*TPA(I,J)-H_PLA*PLA(I,J)*R_L2AP(I,J))

  621   CONTINUE
  631 ENDDO

C-- Reactive wall flux is unknown for PLA
      DO 521 I = 1, NX-1
          PLAT(I,NY) = PLAT(I,NY-1)
  521   CONTINUE


C-------------R_L2AP CALCULATION---------
      D1=Dl2ap

      DO 632 J = 1, NY-1
      DO 622 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(R_L2AP, I, J, CNVR_L2AP)
CC -- DIFFUSION  --
          CALL FDIFF(R_L2AP,R_L2APXIXI,R_L2APXIET
     &                  ,R_L2APETXI,R_L2APETET,I,J,D1)
C
          R_L2APXI = 0.5D0*(R_L2AP(I+1,J) - R_L2AP(I-1,J))
          R_L2APET = 0.5D0*(R_L2AP(I,J+1) - R_L2AP(I,J-1))
          DIFFR_L2AP =( ALPHA(I,J)*R_L2APXIXI-BETA(I,J)
     &      *(R_L2APXIET + R_L2APETXI)+ GAMMA(I,J)*R_L2APETET
     &     +  (RNUT(I,J)+Dl2ap)*( DLX(I,J)*( YXI(I,J)*R_L2APET 
     &      - YET(I,J)*R_L2APXI )+ DLY(I,J)*( XET(I,J)
     &      *R_L2APXI - XXI(I,J)*R_L2APET ) ))*RRJ2(I,J)
C
      R_L2APT(I,J) = R_L2AP(I,J) + DT*( - UM/D*CNVR_L2AP + 1.0/D/D*
     &           DIFFR_L2AP -H_PLA*PLA(I,J)*R_L2AP(I,J))

  622   CONTINUE
  632 ENDDO

C-- Reactive wall flux is unknown for R_L2AP
      DO 522 I = 1, NX-1
          R_L2APT(I,NY) = R_L2APT(I,NY-1)
  522   CONTINUE



C------------------------------------------------

C
      TIM=TIM+DT
      KDEP=19.3888692D-4*TIM-9.3982062D-7*(TIM**2)  

       WRITE(6,*) TIM,TPA(63,50),PLS(63,50),PLA(63,50),R_L2AP(63,50)

       DO 201 J=1,NY
       DO 202 I=1,NX-1
         RP(I,J)=RPT(I,J)
         AP(I,J)=APT(I,J)
         APR(I,J)=APRT(I,J)
         APS(I,J)=APST(I,J)
         PT(I,J)=PTT(I,J)
         T(I,J)=TT(I,J)
         AT(I,J)=ATT(I,J)
         FG(I,J)=FGT(I,J)
         FB(I,J)=FBT(I,J)
         TPA(I,J)=TPAT(I,J)
         PLS(I,J)=PLST(I,J)
         PLA(I,J)=PLAT(I,J)
         R_L2AP(I,J)=R_L2APT(I,J)

  202  ENDDO
  201  ENDDO
C
      DO 113 I = 1, NX-1
        RP(I,0)=RP(I,1)
        AP(I,0)=AP(I,1)
        APR(I,0)=APR(I,1)
        APS(I,0)=APS(I,1)
        PT(I,0)=PT(I,1)
        T(I,0)=T(I,1)
        AT(I,0)=AT(I,1)
        FG(I,0)=FG(I,1)
        FB(I,0)=FB(I,1)
        TPA(I,0)=TPA(I,1)
        PLS(I,0)=PLS(I,1)
        PLA(I,0)=PLA(I,1)
        R_L2AP(I,0)=R_L2AP(I,1)

        RP(I,-1)=RP(I,1)
        AP(I,-1)=AP(I,1)
        APR(I,-1)=APR(I,1)
        APS(I,-1)=APS(I,1)
        PT(I,-1)=PT(I,1)
        T(I,-1)=T(I,1)
        AT(I,-1)=AT(I,1)
        FG(I,-1)=FG(I,1)
        FB(I,-1)=FB(I,1)
        TPA(I,-1)=TPA(I,1)
        PLS(I,-1)=PLS(I,1)
        PLA(I,-1)=PLA(I,1)
        R_L2AP(I,-1)=R_L2AP(I,1)

  113 CONTINUE

       DO 256 J=-1,NY

       RP(0,J)=RP(NX-1,J)
       AP(0,J)=AP(NX-1,J)
       APR(0,J)=APR(NX-1,J)
       APS(0,J)=APS(NX-1,J)
       PT(0,J)=PT(NX-1,J)
       T(0,J)=T(NX-1,J)
       AT(0,J)=AT(NX-1,J)
       FG(0,J)=FG(NX-1,J)
       FB(0,J)=FB(NX-1,J)
       TPA(0,J)=TPA(NX-1,J)
       PLS(0,J)=PLS(NX-1,J)
       PLA(0,J)=PLA(NX-1,J)
       R_L2AP(0,J)=R_L2AP(NX-1,J)

       RP(-1,J)=RP(NX-1,J)
       AP(-1,J)=AP(NX-1,J)
       APR(-1,J)=APR(NX-1,J)
       APS(-1,J)=APS(NX-1,J)
       PT(-1,J)=PT(NX-1,J)
       T(-1,J)=T(NX-1,J)
       AT(-1,J)=AT(NX-1,J)
       FG(-1,J)=FG(NX-1,J)
       FB(-1,J)=FB(NX-1,J)
       TPA(-1,J)=TPA(NX-1,J)
       PLS(-1,J)=PLS(NX-1,J)
       PLA(-1,J)=PLA(NX-1,J)
       R_L2AP(-1,J)=R_L2AP(NX-1,J)

        RP(NX,J)=RP(NX-1,J)
        AP(NX,J)=AP(NX-1,J)
        APR(NX,J)=APR(NX-1,J)
        APS(NX,J)=APS(NX-1,J)
        PT(NX,J)=PT(NX-1,J)
        T(NX,J)=T(NX-1,J)
        AT(NX,J)=AT(NX-1,J)
        FG(NX,J)=FG(NX-1,J)
        FB(NX,J)=FB(NX-1,J)
        TPA(NX,J)=TPA(NX-1,J)
        PLS(NX,J)=PLS(NX-1,J)
        PLA(NX,J)=PLA(NX-1,J)
        R_L2AP(NX,J)=R_L2AP(NX-1,J)

        RP(NX+1,J)=RP(NX-1,J)
        AP(NX+1,J)=AP(NX-1,J)
        APR(NX+1,J)=APR(NX-1,J)
        APS(NX+1,J)=APS(NX-1,J)
        PT(NX+1,J)=PT(NX-1,J)
        T(NX+1,J)=T(NX-1,J)
        AT(NX+1,J)=AT(NX-1,J)
        FG(NX+1,J)=FG(NX-1,J)
        FB(NX+1,J)=FB(NX-1,J)
        TPA(NX+1,J)=TPA(NX-1,J)
        PLS(NX+1,J)=PLS(NX-1,J)
        PLA(NX+1,J)=PLA(NX-1,J)
        R_L2AP(NX+1,J)=R_L2AP(NX-1,J)

  256   ENDDO

      DO 214 I = -1, NX+1
        RP(I,NY+1)=RP(I,NY)
        AP(I,NY+1)=AP(I,NY)
        APR(I,NY+1)=APR(I,NY)
        APS(I,NY+1)=APS(I,NY)
        PT(I,NY+1)=PT(I,NY)
        T(I,NY+1)=T(I,NY)
        AT(I,NY+1)=AT(I,NY)
        FG(I,NY+1)=FG(I,NY)
        FB(I,NY+1)=FB(I,NY)
        TPA(I,NY+1)=TPA(I,NY)
        PLS(I,NY+1)=PLS(I,NY)
        PLA(I,NY+1)=PLA(I,NY)
        R_L2AP(I,NY+1)=R_L2AP(I,NY)
  214  CONTINUE

      DO 350 I=0, NX
          JR(I) = S(I)*KRS*RP(I,NY)
          JA(I) = (S(I)*KAS+(MAS(I)/7.0D6)*KAA)*AP(I,NY)
          JPR(I) = Lambdaj*(THETA*S(I)*KRS*RP(I,NY))
          JPS(I) = MAT(I)*Spj
          JPT(I) = BETACOF*(MAT(I)*(3.22D-8)+MR(I)*FEIRT)*PT(I,NY)
          JT(I) = (MAT(I)*(3.22D-8)+MR(I)*FEIRT)*PT(I,NY)
          JFB(I) = ((KDEP*FB(I,NY))/(RHOS-DELTAFB*KDEP*FB(I,NY)))
     &            *(DELTAPLT*JR(I)+DELTAPLT*JA(I))
  350 CONTINUE

       DO 203 I = 0, NX    
           S(I) = 1.0- M(I)/7.0D6
           M(I) = M(I) + DT*(JR(I)+KAS*S(I)*AP(I,NY))
           MAS(I) = MAS(I) + DT*(THETA*JR(I)+KAS*S(I)*AP(I,NY)
     &             )
           MR(I) = MR(I) + DT*((1-THETA)*JR(I))
           MAT(I) = MAT(I) + DT*(THETA*JR(I)+JA(I))
           MTO(I) = MTO(I) + DT*(JR(I)+JA(I))
  203  CONTINUE

       IF(MOD(KTT,20000).EQ.0) THEN
       K13=KTT
       K12=K13/1000000000
       K11=K13-K12*1000000000
       K10=K11/100000000
       KK9=K11-K10*100000000
       KK8=KK9/10000000
       KK7=KK9-KK8*10000000
       KK6=KK7/1000000
       KK5=KK7-KK6*1000000
       KK4=KK5/100000
       KK3=KK5-KK4*100000
       KK2=KK3/10000
       KK1=KK3-KK2*10000
       KK0=KK1/1000
       KL1=KK1-KK0*1000
       KL2=KL1/100
       KL3=KL1-KL2*100
       KL4=KL3/10
       KL5=KL3-KL4*10
       KL6=KL5/1

       NAME1='RP'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME2='AP'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME3='APR'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME4='APS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME5='PT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME6='T'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME7='AT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME8='KPA'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME11='MTO'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME12='M'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME13='MAS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME14='MR'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME15='MAT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME16='RPB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME17='APB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME18='APRB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME19='APSB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME20='PTB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME21='TB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME22='ATB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME23='FG'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME24='FB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME25='FGB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME26='FBB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME27='TPA'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME28='TPAB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME29='PLS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME30='PLSB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME31='PLA'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME32='PLAB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME33='L2AP'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME34='L2APB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

      OPEN(2010,FILE=NAME1,STATUS='UNKNOWN')
      OPEN(2011,FILE=NAME2,STATUS='UNKNOWN')
      OPEN(2012,FILE=NAME3,STATUS='UNKNOWN')
      OPEN(2013,FILE=NAME4,STATUS='UNKNOWN')
      OPEN(2014,FILE=NAME5,STATUS='UNKNOWN')
      OPEN(2015,FILE=NAME6,STATUS='UNKNOWN')
      OPEN(2016,FILE=NAME7,STATUS='UNKNOWN')
      OPEN(2017,FILE=NAME8,STATUS='UNKNOWN')
      OPEN(2018,FILE=NAME23,STATUS='UNKNOWN')
      OPEN(2019,FILE=NAME24,STATUS='UNKNOWN')
      OPEN(2020,FILE=NAME27,STATUS='UNKNOWN')
      OPEN(2021,FILE=NAME29,STATUS='UNKNOWN')
      OPEN(2022,FILE=NAME31,STATUS='UNKNOWN')
      OPEN(2023,FILE=NAME33,STATUS='UNKNOWN')

      WRITE(2010,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2011,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2012,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2013,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2014,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2015,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2016,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2017,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2018,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2019,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2020,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2021,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2022,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2023,*) 'ZONE F=POINT,I=201,J=51'

      DO 800 J = 0, NY
        DO 800 I = 0, NX

      WRITE(2010,*) XP(I,J)*D,YP(I,J)*D,RP(I,J)
      WRITE(2011,*) XP(I,J)*D,YP(I,J)*D,AP(I,J)
      WRITE(2012,*) XP(I,J)*D,YP(I,J)*D,APR(I,J)
      WRITE(2013,*) XP(I,J)*D,YP(I,J)*D,APS(I,J)
      WRITE(2014,*) XP(I,J)*D,YP(I,J)*D,PT(I,J)
      WRITE(2015,*) XP(I,J)*D,YP(I,J)*D,T(I,J)
      WRITE(2016,*) XP(I,J)*D,YP(I,J)*D,AT(I,J)
      WRITE(2017,*) XP(I,J)*D,YP(I,J)*D,KPA(I,J)
      WRITE(2018,*) XP(I,J)*D,YP(I,J)*D,FG(I,J)
      WRITE(2019,*) XP(I,J)*D,YP(I,J)*D,FB(I,J)
      WRITE(2020,*) XP(I,J)*D,YP(I,J)*D,TPA(I,J)
      WRITE(2021,*) XP(I,J)*D,YP(I,J)*D,PLS(I,J)
      WRITE(2022,*) XP(I,J)*D,YP(I,J)*D,PLA(I,J)
  800 WRITE(2023,*) XP(I,J)*D,YP(I,J)*D,R_L2AP(I,J)



      DO 212 J = 1, NY-1
        DO 212 I = 1, NX-1

       IF(XP(I,J)*D.GE.0.0) THEN

      TRP = TRP + DSS(I,J)*RP(I,J)*2*Pi*YP(I,J)*D
      TAP = TAP + DSS(I,J)*AP(I,J)*2*Pi*YP(I,J)*D
      TAPR = TAPR + DSS(I,J)*APR(I,J)*2*Pi*YP(I,J)*D
      TAPS = TAPS + DSS(I,J)*APS(I,J)*2*Pi*YP(I,J)*D
      TPT = TPT + DSS(I,J)*PT(I,J)*2*Pi*YP(I,J)*D

C    TOTAL THROMBIN CALCULATIION

      TTRB = TTRB + DSS(I,J)*T(I,J)*2*Pi*YP(I,J)*D
      TAT = TAT + DSS(I,J)*AT(I,J)*2*Pi*YP(I,J)*D
      
C     Added Fibrinogen and Fibrin calculation

      TFG = TFG + DSS(I,J)*FG(I,J)*2*Pi*YP(I,J)*D
      TFB = TFB + DSS(I,J)*FB(I,J)*2*Pi*YP(I,J)*D

C     Added Plasminogen and Plasmin Calculation
      TTPA = TTPA + DSS(I,J)*TPA(I,J)*2*Pi*YP(I,J)*D
      TPLS = TPLS + DSS(I,J)*PLS(I,J)*2*Pi*YP(I,J)*D
      TPLA = TPLA + DSS(I,J)*PLA(I,J)*2*Pi*YP(I,J)*D
      TL2AP = TL2AP + DSS(I,J)*R_L2AP(I,J)*2*Pi*YP(I,J)*D


        ENDIF

  212 ENDDO

      DO 235 J = 0, NY
  235 TAPP = TAPP + DSin(J)*AP(0,J)*U(0,J)*UM

       DO 370 I=1,NX-1
 
       IF(XP(I,NY)*D.GE.0.0) THEN

        TDEP = TDEP + DS(I)*MTO(I)*2*Pi*YP(I,NY)*D
        TDEPRP = TDEPRP + DS(I)*MR(I)*2*Pi*YP(I,NY)*D
        TDEPAPP = TDEPAPP + DS(I)*MAS(I)*2*Pi*YP(I,NY)*D
        TDEPAP = TDEPAP + DS(I)*MAT(I)*2*Pi*YP(I,NY)*D

        ENDIF
       
  370   ENDDO

      OPEN(3000,FILE='TRP',STATUS='UNKNOWN')
      OPEN(3100,FILE='TAP',STATUS='UNKNOWN')
      OPEN(3200,FILE='TAPR',STATUS='UNKNOWN')
      OPEN(3300,FILE='TAPS',STATUS='UNKNOWN')
      OPEN(3400,FILE='TPT',STATUS='UNKNOWN')
      OPEN(3500,FILE='TTRB',STATUS='UNKNOWN')
      OPEN(3600,FILE='TAT',STATUS='UNKNOWN')
      OPEN(5500,FILE='TFG',STATUS='UNKNOWN')
      OPEN(5600,FILE='TFB',STATUS='UNKNOWN')
      OPEN(5700,FILE='TTPA',STATUS='UNKNOWN')
      OPEN(5800,FILE='TPLS',STATUS='UNKNOWN')
      OPEN(5900,FILE='TPLA',STATUS='UNKNOWN')
      OPEN(6000,FILE='TL2AP',STATUS='UNKNOWN')

      OPEN(3700,FILE='TDEPOSITION',STATUS='UNKNOWN')
      OPEN(3800,FILE='TDEPRP',STATUS='UNKNOWN')
      OPEN(3900,FILE='TDEPAP',STATUS='UNKNOWN')
      OPEN(5000,FILE='TDEPAPP',STATUS='UNKNOWN')
      OPEN(5100,FILE='MAT150',STATUS='UNKNOWN')
      OPEN(5200,FILE='APTNumRate',STATUS='UNKNOWN')
      OPEN(5300,FILE='APpoint',STATUS='UNKNOWN')
      OPEN(5400,FILE='Information Needed',STATUS='UNKNOWN')


      WRITE(3000,*) TRP
      WRITE(3100,*) TAP
      WRITE(3200,*) TAPR
      WRITE(3300,*) TAPS
      WRITE(3400,*) TPT
      WRITE(3500,*) TTRB
      WRITE(3600,*) TAT
      WRITE(5500,*) TFG
      WRITE(5600,*) TFB
      WRITE(5700,*) TTPA
      WRITE(5800,*) TPLS
      WRITE(5900,*) TPLA
      WRITE(6000,*) TL2AP

      WRITE(3700,*) TDEP
      WRITE(3800,*) TDEPRP
      WRITE(3900,*) TDEPAP
      WRITE(5000,*) TDEPAPP
      WRITE(5100,*)  MAT(150)
      WRITE(5200,*) TAPP
      WRITE(5300,*) AP(138,49),AP(138,50)
      WRITE(5400,*) AP(157,49),AP(157,50),MAS(157)


        OPEN(50,FILE=NAME11,STATUS='UNKNOWN')
        WRITE(50,*) MTO

        OPEN(51,FILE=NAME12,STATUS='UNKNOWN')
        WRITE(51,*) M

        OPEN(52,FILE=NAME13,STATUS='UNKNOWN')
        WRITE(52,*) MAS

        OPEN(53,FILE=NAME14,STATUS='UNKNOWN')
        WRITE(53,*) MR

        OPEN(54,FILE=NAME15,STATUS='UNKNOWN')
        WRITE(54,*) MAT

      OPEN(4000,FILE=NAME16,STATUS='UNKNOWN')
      OPEN(4100,FILE=NAME17,STATUS='UNKNOWN')
      OPEN(4200,FILE=NAME18,STATUS='UNKNOWN')
      OPEN(4300,FILE=NAME19,STATUS='UNKNOWN')
      OPEN(4400,FILE=NAME20,STATUS='UNKNOWN')
      OPEN(4500,FILE=NAME21,STATUS='UNKNOWN')
      OPEN(4600,FILE=NAME22,STATUS='UNKNOWN')
      OPEN(4700,FILE=NAME25,STATUS='UNKNOWN')
      OPEN(4800,FILE=NAME26,STATUS='UNKNOWN')
      OPEN(4810,FILE=NAME28,STATUS='UNKNOWN')
      OPEN(4820,FILE=NAME30,STATUS='UNKNOWN')
      OPEN(4830,FILE=NAME32,STATUS='UNKNOWN')
      OPEN(4840,FILE=NAME34,STATUS='UNKNOWN')

        DO 900 I = 0, NX

      WRITE(4000,*) RP(I,NY)
      WRITE(4100,*) AP(I,NY)
      WRITE(4200,*) APR(I,NY)
      WRITE(4300,*) APS(I,NY)
      WRITE(4400,*) PT(I,NY)
      WRITE(4500,*) T(I,NY)
      WRITE(4600,*) AT(I,NY)
      WRITE(4700,*) FG(I,NY)
      WRITE(4800,*) FB(I,NY)
      WRITE(4810,*) TPA(I,NY)
      WRITE(4820,*) PLS(I,NY)
      WRITE(4830,*) PLA(I,NY)
  900 WRITE(4840,*) R_L2AP(I,NY)

       ENDIF

      KTT=KTT+1

      TS = 0D0
      TRP = 0D0
      TAP = 0D0
      TAPP = 0D0
      TAPR = 0D0
      TAPS = 0D0
      TPT = 0D0
      TTRB = 0D0
      TAT = 0D0
      TFG = 0D0
      TFB = 0D0
      TTPA = 0D0
      TPLS = 0D0
      TPLA = 0D0
      TL2AP = 0D0
      TDEP = 0D0
      TDEPRP = 0D0
      TDEPAP = 0D0
      TDEPAPP = 0D0

       IF(TIM.LT.ALLT) THEN
         GOTO 300
        ENDIF

C---------------------------------------------------------

      END
C
C---------------------SUBROUTINE-------------------------
C************************************************
      SUBROUTINE CNVUP(A, I, J, CNVA)
C************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NX=200,NY=50)
      COMMON /BLK1/ RP(-1:NX+1,-1:NY+1),AP(-1:NX+1,-1:NY+1),
     &              APR(-1:NX+1,-1:NY+1),APS(-1:NX+1,-1:NY+1),
     &              PT(-1:NX+1,-1:NY+1),T(-1:NX+1,-1:NY+1),
     &              AT(-1:NX+1,-1:NY+1),FG(-1:NX+1,-1:NY+1),
     &              FB(-1:NX+1,-1:NY+1),TPA(-1:NX+1,-1:NY+1),
     &              PLS(-1:NX+1,-1:NY+1),PLA(-1:NX+1,-1:NY+1),
     &              R_L2AP(-1:NX+1,-1:NY+1)
      COMMON /BLK2/ U(0:NX,0:NY),V(0:NX,0:NY)
      COMMON /BLK3/ RPT(0:NX,0:NY),APT(0:NX,0:NY),
     &              APRT(0:NX,0:NY),APST(0:NX,0:NY) ,PTT(0:NX,0:NY),
     &              TT(0:NX,0:NY) ,ATT(0:NX,0:NY),
     &              FGT(0:NX,0:NY) ,FBT(0:NX,0:NY),
     &              TPAT(0:NX,0:NY),PLST(0:NX,0:NY),PLAT(0:NX,0:NY),
     &              R_L2APT(0:NX,0:NY)
      COMMON /BLK4/ XP(0:NX,0:NY), YP(0:NX,0:NY)
      COMMON /BLK5/ RJ(0:NX,0:NY), ALPHA(0:NX,0:NY),
     &              BETA(0:NX,0:NY), GAMMA(0:NX,0:NY)
      COMMON /BLK6/ ARJYTE(0:NX,0:NY), GRJYTN(0:NX,0:NY),
     &              DLX(0:NX,0:NY), DLY(0:NX,0:NY), YWL(0:NX,0:NY)
      COMMON /BLK7/ UC(0:NX,0:NY), VC(0:NX,0:NY)
      COMMON /BLK8/ XXI(0:NX,0:NY), YXI(0:NX,0:NY),
     &              XET(0:NX,0:NY), YET(0:NX,0:NY)
      COMMON /BLK9/ UP(0:NX,0:NY),VP(0:NX,0:NY)
      COMMON /BLK10/ OMEGA(0:NX,0:NY),RNUT(0:NX,0:NY),
     &              NU(0:NX,0:NY) 


C
C----------------------------------------------------------------------
      DIMENSION A(-1:NX+1,-1:NY+1)
      YPRJ = 0.5D0/(YP(I,J)*RJ(I,J))
      AC = A(I,J)
      AR = A(I+1,J)
      AL = A(I-1,J)
      AT1 = A(I,J+1)
      AB = A(I,J-1)
      FAR = UC(I  ,J)*(AC + AR) + ABS(UC(I  ,J))*(AC - AR)
      FAL = UC(I-1,J)*(AC + AL) - ABS(UC(I-1,J))*(AC - AL)
      FAT = VC(I,J  )*(AC + AT1) + ABS(VC(I,J  ))*(AC - AT1)
      FAB = VC(I,J-1)*(AC + AB) - ABS(VC(I,J-1))*(AC - AB)

      CNVA = ( FAR - FAL + FAT - FAB )*YPRJ

      END
C
C******************************************************
      SUBROUTINE FDIFF(F,FXIXI,FXIET,FETXI,FETET,I,J,C)
C******************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NX=200,NY=50)
      COMMON /BLK1/ RP(-1:NX+1,-1:NY+1),AP(-1:NX+1,-1:NY+1),
     &              APR(-1:NX+1,-1:NY+1),APS(-1:NX+1,-1:NY+1),
     &              PT(-1:NX+1,-1:NY+1),T(-1:NX+1,-1:NY+1),
     &              AT(-1:NX+1,-1:NY+1),FG(-1:NX+1,-1:NY+1),
     &              FB(-1:NX+1,-1:NY+1),TPA(-1:NX+1,-1:NY+1),
     &              PLS(-1:NX+1,-1:NY+1),PLA(-1:NX+1,-1:NY+1),
     &              R_L2AP(-1:NX+1,-1:NY+1)
      COMMON /BLK2/ U(0:NX,0:NY),V(0:NX,0:NY)
      COMMON /BLK3/ RPT(0:NX,0:NY),APT(0:NX,0:NY),
     &              APRT(0:NX,0:NY),APST(0:NX,0:NY) ,PTT(0:NX,0:NY),
     &              TT(0:NX,0:NY) ,ATT(0:NX,0:NY),
     &              FGT(0:NX,0:NY) ,FBT(0:NX,0:NY),
     &              TPAT(0:NX,0:NY),PLST(0:NX,0:NY),PLAT(0:NX,0:NY),
     &              R_L2APT(0:NX,0:NY)
      COMMON /BLK4/ XP(0:NX,0:NY), YP(0:NX,0:NY)
      COMMON /BLK5/ RJ(0:NX,0:NY), ALPHA(0:NX,0:NY),
     &              BETA(0:NX,0:NY), GAMMA(0:NX,0:NY)
      COMMON /BLK6/ ARJYTE(0:NX,0:NY), GRJYTN(0:NX,0:NY),
     &              DLX(0:NX,0:NY), DLY(0:NX,0:NY), YWL(0:NX,0:NY)
      COMMON /BLK7/ UC(0:NX,0:NY), VC(0:NX,0:NY)
      COMMON /BLK8/ XXI(0:NX,0:NY), YXI(0:NX,0:NY),
     &              XET(0:NX,0:NY), YET(0:NX,0:NY)
      COMMON /BLK9/ UP(0:NX,0:NY),VP(0:NX,0:NY)
      COMMON /BLK10/ OMEGA(0:NX,0:NY),RNUT(0:NX,0:NY),
     &              NU(0:NX,0:NY) 

C
      REAL*8    C 

C
      DIMENSION F(-1:NX+1,-1:NY+1)
C----------------------------------------------------------------------
C
      RNUTR= 0.5D0*(RNUT(I,J)+RNUT(I+1,J)) 
      RNUTL= 0.5D0*(RNUT(I,J)+RNUT(I-1,J)) 
      RNUTT= 0.5D0*(RNUT(I,J)+RNUT(I,J+1)) 
      RNUTB= 0.5D0*(RNUT(I,J)+RNUT(I,J-1)) 
C
      FXIR = F(I+1,J) - F(I,J)
      FXIL = F(I,J) - F(I-1,J)


      FXIXI = (RNUTR+C)*FXIR - (RNUTL+C)*FXIL
C
      FXIT = 0.5D0*(F(I+1,J+1) - F(I-1,J+1))
      FXIB = 0.5D0*(F(I+1,J-1) - F(I-1,J-1))


      FXIET = 0.5D0*((RNUTT+C)*FXIT - (RNUTB+C)*FXIB)
C
      FETR = 0.5D0*(F(I+1,J+1) - F(I+1,J-1))
      FETL = 0.5D0*(F(I-1,J+1) - F(I-1,J-1))


      FETXI = 0.5D0*((RNUTR+C)*FETR - (RNUTL+C)*FETL)
C
      FETT = F(I,J+1) - F(I,J)
      FETB = F(I,J) - F(I,J-1)


      FETET = (RNUTT+C)*FETT - (RNUTB+C)*FETB
C
      RETURN

      END
C



