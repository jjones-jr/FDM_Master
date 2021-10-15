C************************************************
C     PROGRAM PLATELET DEPOSITION AND ACTIVATION
C     CODED BY J.E.Jones III
C     2021.10
C************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NX=200,NY=50)
      COMMON /BLK1/ RP(-1:NX+1,-1:NY+1),AP(-1:NX+1,-1:NY+1),
     &              APR(-1:NX+1,-1:NY+1),APS(-1:NX+1,-1:NY+1),
     &              f_I(-1:NX+1,-1:NY+1),f_IA(-1:NX+1,-1:NY+1),
     &              f_II(-1:NX+1,-1:NY+1),f_IIA(-1:NX+1,-1:NY+1),
     &              f_V(-1:NX+1,-1:NY+1),f_VA(-1:NX+1,-1:NY+1),
     &              f_VIII(-1:NX+1,-1:NY+1),f_VIIIA(-1:NX+1,-1:NY+1),
     &              f_IX(-1:NX+1,-1:NY+1),f_IXA(-1:NX+1,-1:NY+1),
     &              f_X(-1:NX+1,-1:NY+1),f_XA(-1:NX+1,-1:NY+1),
     &              f_XI(-1:NX+1,-1:NY+1),f_XIA(-1:NX+1,-1:NY+1),
     &              ATIII(-1:NX+1,-1:NY+1),Z(-1:NX+1,-1:NY+1),
     &              W(-1:NX+1,-1:NY+1),TFPI(-1:NX+1,-1:NY+1),
     &              APC(-1:NX+1,-1:NY+1),PC(-1:NX+1,-1:NY+1),
     &              rL1AT(-1:NX+1,-1:NY+1),TPA(-1:NX+1,-1:NY+1),
     &              PLA(-1:NX+1,-1:NY+1),PLS(-1:NX+1,-1:NY+1),
     &              rL2AP(-1:NX+1,-1:NY+1)
      COMMON /BLK2/ U(0:NX,0:NY),V(0:NX,0:NY)
      COMMON /BLK3/ RPT(0:NX,0:NY),APT(0:NX,0:NY),
     &              APRT(0:NX,0:NY),APST(0:NX,0:NY), 
     &              f_IT(0:NX,0:NY),f_IAT(0:NX,0:NY),
     &              f_IIT(0:NX,0:NY),f_IIAT(0:NX,0:NY),
     &              f_VT(0:NX,0:NY),f_VAT(0:NX,0:NY),
     &              f_VIIIT(0:NX,0:NY),f_VIIIAT(0:NX,0:NY),
     &              f_IXT(0:NX,0:NY),f_IXAT(0:NX,0:NY),
     &              f_XT(0:NX,0:NY),f_XAT(0:NX,0:NY),
     &              f_XIT(0:NX,0:NY),f_XIAT(0:NX,0:NY),
     &              ATIIIT(0:NX,0:NY),ZT(0:NX,0:NY),
     &              WT(0:NX,0:NY),TFPIT(0:NX,0:NY),
     &              APCT(0:NX,0:NY),PCT(0:NX,0:NY),
     &              rL1ATT(0:NX,0:NY),TPAT(0:NX,0:NY),
     &              PLAT(0:NX,0:NY),PLST(0:NX,0:NY),
     &              rL2APT(0:NX,0:NY)
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
        REAL*8        JR(0:NX),JA(0:NX),JPR(0:NX),JPS(0:NX),
     &              Jf_II(0:NX),Jf_IIA(0:NX),Jf_IA(0:NX),
     &              DYB(0:NX),S(0:NX),
     &              M(0:NX),MTO(0:NX),MAS(0:NX),
     &              MR(0:NX), MAT(0:NX),DS(0:NX),
     &              BGAMMA(0:NX,0:NY),TAUD(0:NX,0:NY),
     &              DSS(0:NX,0:NY),RRJ2(0:NX,0:NY),
     &              KPA(0:NX,0:NY),F(0:NX,0:NY),OMEGA1(0:NX,0:NY),
    !  &              K8_9(0:NX,0:NY),H8_9(0:NX,0:NY),
    !  &              K5_10(0:NX,0:NY),H5_10(0:NX,0:NY),
     &              DT,ALLT,TACT,K_9,H_9,K_8,H_8,K_A,K_5,H_5,
     &              K_10,H_10,K_TFPI,K_2,K_2M,H_2,KAPC,HAPC,K_1,K_1M,
     &              H1,H_1M,K_11,H_11,KIIA_AP,KAP_AP,KTPA_IA_PLA,
     &              HPLA,DRP,DAP,Dapr,Daps,Df_I,Df_Ia,
     &              Df_II,Df_IIA,Df_V,Df_Va,Df_VIII,
     &              Df_VIIIa,Df_IX,Df_IXa,Df_X,Df_Xa,
     &              Df_XI,Df_XIa,DPC,DAPC,DATIII,DZ,DW
     &              DTFPI,DrL1AT,DtPA,DPLS,DPLA,DrL2AP,
     &              K1_T,KAT,KT,ALPHACOF,BETACOF,FEIAT,FEIRT,TIM,
     &              K1_J,Lambdaj,Spj,THETA,H,KRS,KAS,KAA,TAUDD,Re,
     &              UM, DSin(0:NY),TS,RP_init,N,D1,
     &              KDEP,RHOS,DELTAPLT,DELTAFB,
     &              Z_bunshi,Z_bunbo_1,Z_bunbo_2

      CHARACTER NAME1*25,NAME2*25,NAME3*25,NAME4*25
      CHARACTER NAME5*25,NAME6*25,NAME7*25,NAME8*25
      CHARACTER NAME9*25,NAME10*25,NAME11*25
      CHARACTER NAME12*25,NAME13*25,NAME14*25
      CHARACTER NAME15*25,NAME16*25,NAME17*25
      CHARACTER NAME18*25,NAME19*25,NAME20*25
      CHARACTER NAME21*25,NAME22*25,NAME23*25,NAME24*25
      CHARACTER NAME25*25,NAME26*25,NAME27*25,NAME28*25
      CHARACTER NAME29*25
      CHARACTER NAME30*25,NAME31*25,NAME32*25,NAME33*25
      CHARACTER NAME34*25,NAME35*25,NAME36*25,NAME37*25
      CHARACTER NAME38*25,NAME39*25
      CHARACTER NAME40*25,NAME41*25,NAME42*25,NAME43*25
      CHARACTER NAME44*25,NAME45*25,NAME46*25,NAME47*25
      CHARACTER NAME48*25,NAME49*25
      CHARACTER NAME50*25,NAME51*25,NAME52*25,NAME53*25
      CHARACTER NAME54*25,NAME55*25,NAME56*25,NAME57*25
      CHARACTER NAME58*25,NAME59*25,NAME60*25,NAME61*25
      CHARACTER NAME62*25,NAME63*25,NAME64*25
   

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

      ! WRITE(6,*) 'PLEASE Input Re Number'
      ! READ(5,*) Re
      Re = 5300

        UM = Re/3180.0D0*13.159973D0

      WRITE(6,*) 'RE=',RE,'UM=', UM

      ! WRITE(6,*) 'PLEASE Choose Re Number:3180,4240,5300,6360,7420'

      Re = 5300

      ! WRITE(6,*) 'PLEASE Choose AB:0 or CB 1'
      ! READ(5,*) IFLAG
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
      Tf_IXA = 0D0
      Tf_IX = 0D0
      Tf_VIIIA = 0D0
      Tf_VIII = 0D0
      Tf_VA = 0D0
      Tf_V = 0D0
      TZ = 0D0
      TW = 0D0
      Tf_XA = 0D0
      Tf_X = 0D0
      Tf_IIA = 0D0
      Tf_II = 0D0
      Tf_IA = 0D0
      Tf_I = 0D0
      Tf_XIA = 0D0
      Tf_XI = 0D0
      TATIII = 0D0
      TTFPI = 0D0
      TAPC = 0D0
      TPC = 0D0
      TrL1AT = 0D0
      TPLA = 0D0
      TPLS = 0D0
      TrL2AP = 0D0

      TDEP = 0D0
      TDEPRP = 0D0
      TDEPAP = 0D0
      TDEPAPP = 0D0

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

      ! WRITE(6,*) 'PLEASE INPUT Threshold Shear Rate Or 0 means no limit'
      ! READ(5,*) STH
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

      DRP=1.58D-9
      DAP=1.58D-9
      Dapr=2.57D-6
      Daps=2.14D-6
      Df_I=3.1D-7
      Df_Ia=2.47D-7
      Df_II=5.21D-7
      Df_IIa=6.47D-7
      Df_V=3.12D-7
      Df_Va=3.82D-7
      Df_VIII=3.12D-7
      Df_VIIIa=3.92D-7
      Df_IX=5.63D-7
      Df_IXa=6.25D-7
      Df_X=5.63D-7
      Df_Xa=7.37D-7
      Df_XI=3.97D-7
      Df_XIa=5.0D-7
      DPC=5.44D-7
      DAPC=5.50D-7
      DATIII=5.57D-7
      DTFPI=6.30D-7
      DrL1AT=5.82D-7
      DtPA=5.28D-7
      DPLS=4.81D-7
      DPLA=4.93D-7
      DrL2AP=5.25D-7
      DZ=0
      DW=0

      K1_T=13.333
      KAT=0.100
      KT=3.50D-2
      K_9=20
      H_9=0.2
      K_8=1.0D-5
      H_8=0.31
      K_A=1.2D0
      K_5=0.17
      H_5=0.31
      K_10=1200
      H_10=1
      K_TFPI=0.96
      K_2=1344
      K_2M=1060
      H_2=1.3
      KAPC=0.0014
      HAPC=0.1
      K_1=3540
      K_1M=3160
      H_1=6.97
      H_1M=2900
      K_11=0.0078
      H_11=0.2
      KIIA_AP=30
      KAP_AP=18
      KTPA_IA_PLA=12
      HPLA=0.096
      RP_init = 2.0D8

      ALPHACOF=1.0
      BETACOF=9.11D-3
      FEIAT=3.69D-9
      FEIRT=6.50D-10
      K1_J=0.0161
      Lambdaj=2.4D-8
      Spj=9.5D-12

C--  Values for Jf_IA(JFB) from Fernolendt(2011)
      RHOS=1050
      DELTAPLT=2.0D-12
      DELTAFB=66D-5

      DT=1.0D-6

c      WRITE(6,*) 'PLEASE Confirm whether first cal. YES:0 ,NO:1'
c      READ(5,*) IEG
      IEG = 0
c       IF(IEG.EQ.0) THEN

      DO 10 I=-1,NX+1
       DO 11 J=-1,NY+1
       RP(I,J)=0.0
       AP(I,J)=0.0
       APR(I,J)=0.0
       APS(I,J)=0.0
       f_IX(I,J)=0.0
       f_IXA(I,J)=0.0
       f_VIIIA(I,J)=0.0
       f_VIII(I,J)=0.0
       f_VA(I,J)=0.0
       f_V(I,J)=0.0
       Z(I,J)=0.0
       W(I,J)=0.0
       f_XA(I,J)=0.0
       f_X(I,J)=0.0
       f_IIA(I,J)=0.0
       f_II(I,J)=0.0
       f_IA(I,J)=0.0
       f_I(I,J)=0.0
       f_XIA(I,J)=0.0
       f_XI(I,J)=0.0
       ATIII(I,J)=0.0
       TFPI(I,J)=0.0
       APC(I,J)=0.0
       PC(I,J)=0.0
       rL1AT(I,J)=0.0
       TPA(I,J)=0.0
       PLA(I,J)=0.0
       PLS(I,J)=0.0
       rL2AP(I,J)=0.0
  11   ENDDO
  10  ENDDO
C
C--Initial concentration(micro-M)
C--Reference: Anand(2003), Table IV

C      DO 114 I=5,10
      DO 114 I=-1,NX+1
       DO 115 J=-1,NY+1
C----RP, AP [PLT/mL]
C---- Conversion between [PLT/mL] to [nM]
C---- 1 = 1.6606*10^-12 nM/(PLT/mL)
       RP(I,J)=RP_init
       AP(I,J)=1.0D7

C----APR,APS [µM]       
       APR(I,J)=0.0
       APS(I,J)=0.0

C----Coagulants [nM]
       f_IXA(I,J)=0.0       
       f_IX(I,J)=90
       f_VIIIA(I,J)=0.0
       f_VIII(I,J)=0.7
       f_VA(I,J)=0.0
       f_V(I,J)=20
       Z(I,J)=0.0
       W(I,J)=0.0
       f_XA(I,J)=0.0
       f_X(I,J)=170
       f_IIA(I,J)=0.0
       f_II(I,J)=1400
       f_IA(I,J)=0.0
       f_I(I,J)=7000
       f_XIA(I,J)=0.0
       f_XI(I,J)=30
       ATIII(I,J)=2410
       TFPI(I,J)=2.5
       APC(I,J)=0.0
       PC(I,J)=60
       rL1AT(I,J)=45000
       TPA(I,J)=0.08
       PLA(I,J)=0.0
       PLS(I,J)=2180
       rL2AP(I,J)=105
  115   ENDDO
  114  ENDDO


      ! DO 31 I=5,10
      DO 31 I=-1,NX+1
       RP(I,NY+1)=RP(I,NY)
       AP(I,NY+1)=AP(I,NY)
       APR(I,NY+1)=APR(I,NY)
       APS(I,NY+1)=APS(I,NY)
       f_IX(I,NY+1)=f_IX(I,NY)
       f_IXA(I,NY+1)=f_IXA(I,NY)
       f_VIIIA(I,NY+1)=f_VIIIA(I,NY)
       f_VIII(I,NY+1)=f_VIII(I,NY)
       f_VA(I,NY+1)=f_VA(I,NY)
       f_V(I,NY+1)=f_V(I,NY)
       Z(I,NY+1)=Z(I,NY)
       W(I,NY+1)=W(I,NY)
       f_XA(I,NY+1)=f_XA(I,NY)
       f_X(I,NY+1)=f_X(I,NY)
       f_IIA(I,NY+1)=f_IIA(I,NY)
       f_II(I,NY+1)=f_II(I,NY)
       f_IA(I,NY+1)=f_IA(I,NY)
       f_I(I,NY+1)=f_I(I,NY)
       f_XIA(I,NY+1)=f_XIA(I,NY)
       f_XI(I,NY+1)=f_XI(I,NY)
       ATIII(I,NY+1)=ATIII(I,NY)
       TFPI(I,NY+1)=TFPI(I,NY)
       APC(I,NY+1)=APC(I,NY)
       PC(I,NY+1)=PC(I,NY)
       rL1AT(I,NY+1)=rL1AT(I,NY)
       TPA(I,NY+1)=TPA(I,NY)
       PLA(I,NY+1)=PLA(I,NY)
       PLS(I,NY+1)=PLS(I,NY)
       rL2AP(I,NY+1)=rL2AP(I,NY)
  31  ENDDO

C--Surface-flux model  
      TIM=0D0
      KDEP=0D0

      DO 12 I=0, NX
         JR(I) = 0.0
         JA(I) = 0.0
         JPR(I) = 0.0
         JPS(I) = 0.0
         Jf_II(I) = 0.0
         Jf_IIA(I) = 0.0
         S(I) = 1.0
         M(I) = 0.0
         MAS(I) = 0.0
         MR(I) = 0.0
         MAT(I) = 0.0
c         Jf_I(I)=0.0
         Jf_IA(I)=0.0

  12  CONTINUE

C---Moved to IEG=ELSE.f

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
     &       + APS(I,J)/0.6 + (f_IIA(I,J)/(1000*BETACOF))/0.1)*F(I,J)
c-----Beta coefficient was used to convert f_IIA(I,J) [nM] to [U/mL]
C-----Note that Factor IIa is Thrombin.

       KPA(I,J) = OMEGA(I,J)/TACT

       ELSE
          KPA(I,J) = 0D0
        ENDIF

       ELSE

       IF(XP(I,J)*D.GE.-5.0.AND.XP(I,J)*D.LE.14.0) THEN

       OMEGA(I,J) = (APR(I,J)/2.0 
     &       + APS(I,J)/0.6 + (f_IIA(I,J)/(1000*BETACOF))/0.1)*F(I,J)
       KPA(I,J) = OMEGA(I,J)/TACT

       ELSE
          KPA(I,J) = 0D0
        ENDIF

        ENDIF

  440   CONTINUE
  430 CONTINUE

C--------------CONCENTRATION CALCULATION---------------


C-----------------Factor IXa CALCULATION-------------
      D1=Df_IXa

      DO 610 J = 1, NY-1
      DO 600 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_IXA, I, J, CNVf_IXA)
CC -- DIFFUSION  --
          CALL FDIFF(f_IXA,f_IXAXIXI,f_IXAXIET,
     &          f_IXAETXI,f_IXAETET,I,J,D1)
C
          f_IXAXI = 0.5D0*(f_IXA(I+1,J) - f_IXA(I-1,J))
          f_IXAET = 0.5D0*(f_IXA(I,J+1) - f_IXA(I,J-1))
          DIFFf_IXA =( ALPHA(I,J)*f_IXAXIXI
     &      -BETA(I,J)*(f_IXAXIET + f_IXAETXI)
     &      +GAMMA(I,J)*f_IXAETET
     &      +(RNUT(I,J)+Df_IXa)*( DLX(I,J)*( YXI(I,J)*f_IXAET 
     &      -YET(I,J)*f_IXAXI )
     &      +DLY(I,J)*( XET(I,J)*f_IXAXI-XXI(I,J)*f_IXAET ))
     &      )*RRJ2(I,J)
C
          f_IXAT(I,J) = f_IXA(I,J) + DT*(-UM/D*CNVf_IXA 
     &      +1.0/D/D*DIFFf_IXA 
     &      +K_9*f_XIA(I,J)-H_9*f_IXA(I,J))

  600   CONTINUE
  610 ENDDO

      DO 500 I = 1, NX-1
          f_IXAT(I,NY) = f_IXAT(I,NY-1)
  500   CONTINUE  

C-----Factor IX (Christmas factor, plasma thromboplastin component, PTC) CALCULATION----
      D1=Df_IX

      DO 611 J = 1, NY-1
      DO 601 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_IX, I, J, CNVf_IX)
CC -- DIFFUSION  --
          CALL FDIFF(f_IX,f_IXXIXI,f_IXXIET,f_IXETXI,f_IXETET,I,J,D1)
C
          f_IXXI = 0.5D0*(f_IX(I+1,J) - f_IX(I-1,J))
          f_IXET = 0.5D0*(f_IX(I,J+1) - f_IX(I,J-1))
          DIFFf_IX =( ALPHA(I,J)*f_IXXIXI
     &      - BETA(I,J)*(f_IXXIET + f_IXETXI)
     &      + GAMMA(I,J)*f_IXETET+ (RNUT(I,J)+Df_IX)
     &      *( DLX(I,J)*( YXI(I,J)*f_IXET
     &     - YET(I,J)*f_IXXI )+ DLY(I,J)*( XET(I,J)*f_IXXI 
     &      - XXI(I,J)*f_IXET )))*RRJ2(I,J)
C
          f_IXT(I,J) = f_IX(I,J) + DT*( - UM/D*CNVf_IX + 1.0/D/D*
     &           DIFFf_IX -K_9*f_XIA(I,J) )

  601   CONTINUE
  611 ENDDO

      DO 501 I = 1, NX-1
          f_IXT(I,NY) = f_IXT(I,NY-1)
  501   CONTINUE

C----------Factor VIIIa CALCULATION------
      D1=Df_VIIIa

      DO 612 J = 1, NY-1
      DO 602 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_VIIIA, I, J, CNVf_VIIIA)
CC -- DIFFUSION  --
          CALL FDIFF(f_VIIIA,f_VIIIAXIXI,f_VIIIAXIET,
     &   f_VIIIAETXI,f_VIIIAETET,I,J,D1)
C
          f_VIIIAXI = 0.5D0*(f_VIIIA(I+1,J) - f_VIIIA(I-1,J))
          f_VIIIAET = 0.5D0*(f_VIIIA(I,J+1) - f_VIIIA(I,J-1))
          DIFFf_VIIIA =( ALPHA(I,J)*f_VIIIAXIXI 
     &      - BETA(I,J)*(f_VIIIAXIET + f_VIIIAETXI)
     &      +GAMMA(I,J)*f_VIIIAETET+ (RNUT(I,J)+Df_VIIIa)
     &      *( DLX(I,J)*( YXI(I,J)*f_VIIIAET 
     &      -YET(I,J)*f_VIIIAXI )+ DLY(I,J)*( XET(I,J)*f_VIIIAXI
     &      -XXI(I,J)*f_VIIIAET)))*RRJ2(I,J)

C
          f_VIIIAT(I,J) = f_VIIIA(I,J) + DT*( 
     &      - UM/D*CNVf_VIIIA + 1.0/D/D
     &      *DIFFf_VIIIA
     &      +K_8*f_IIA(I,J)
     &      -H_8*f_VIIIA(I,J)
     &      -K_A*APC(I,J)*( f_VIIIA(I,J) + Z(I,J)))

  602   CONTINUE
  612 ENDDO

      DO 502 I = 1, NX-1
          f_VIIIAT(I,NY) = f_VIIIAT(I,NY-1)
  502   CONTINUE  


C----------Factor VIII (Antihemophilc globulin, AHG) CALCULATION------
      D1=Df_VIII

      DO 613 J = 1, NY-1
      DO 603 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_VIII, I, J, CNVf_VIII)
CC -- DIFFUSION  --
          CALL FDIFF(f_VIII,f_VIIIXIXI,f_VIIIXIET,
     &  f_VIIIETXI,f_VIIIETET,I,J,D1)
C
          f_VIIIXI = 0.5D0*(f_VIII(I+1,J) - f_VIII(I-1,J))
          f_VIIIET = 0.5D0*(f_VIII(I,J+1) - f_VIII(I,J-1))
          DIFFf_VIII =( ALPHA(I,J)*f_VIIIXIXI 
     &      - BETA(I,J)*(f_VIIIXIET + f_VIIIETXI)
     &      +GAMMA(I,J)*f_VIIIETET+ (RNUT(I,J)+Df_VIII)*
     &      ( DLX(I,J)*( YXI(I,J)*f_VIIIET 
     &      -YET(I,J)*f_VIIIXI )+ DLY(I,J)*( XET(I,J)*f_VIIIXI
     &      -XXI(I,J)*f_VIIIET)))*RRJ2(I,J)

C
          f_VIIIT(I,J) = f_VIII(I,J) + DT*( 
     &           - UM/D*CNVf_VIII + 1.0/D/D*
     &           DIFFf_VIII -K_8*f_IIA(I,J))

  603   CONTINUE
  613 ENDDO

      DO 503 I = 1, NX-1
          f_VIIIT(I,NY) = f_VIIIT(I,NY-1)
  503   CONTINUE

C-----------------Factor Va CALCULATION-------------
      D1=Df_Va

      DO 614 J = 1, NY-1
      DO 604 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_VA, I, J, CNVf_VA)
CC -- DIFFUSION  --
          CALL FDIFF(f_VA,f_VAXIXI,f_VAXIET,f_VAETXI,f_VAETET,I,J,D1)
C
          f_VAXI = 0.5D0*(f_VA(I+1,J) - f_VA(I-1,J))
          f_VAET = 0.5D0*(f_VA(I,J+1) - f_VA(I,J-1))
          DIFFf_VA =( ALPHA(I,J)*f_VAXIXI 
     &     - BETA(I,J)*(f_VAXIET + f_VAETXI)
     &      +GAMMA(I,J)*f_VAETET+ (RNUT(I,J)
     &       +Df_Va)*( DLX(I,J)*( YXI(I,J)*f_VAET 
     &      - YET(I,J)*f_VAXI )+ DLY(I,J)*( XET(I,J)*f_VAXI 
     &      - XXI(I,J)*f_VAET )))*RRJ2(I,J)
C
          f_VAT(I,J) = f_VA(I,J) + DT*(-UM/D*CNVf_VA 
     &      +1.0/D/D*DIFFf_VA +K_5*f_IIA(I,J)
     &      -H_5*f_VA(I,J)-K_A*APC(I,J)*(f_VA(I,J) 
     &      +W(I,J)))

  604   CONTINUE
  614 ENDDO

      DO 504 I = 1, NX-1
          f_VAT(I,NY) = f_VAT(I,NY-1)
  504   CONTINUE

C-----------------Factor fV (Proaccelerin) CALCULATION-------------
      D1=Df_V
      DO 615 J = 1, NY-1
      DO 605 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_V, I, J, CNVf_V)
CC -- DIFFUSION  --
          CALL FDIFF(f_V,f_VXIXI,f_VXIET,f_VETXI,f_VETET,I,J,D1)
C
          f_VXI = 0.5D0*(f_V(I+1,J) - f_V(I-1,J))
          f_VET = 0.5D0*(f_V(I,J+1) - f_V(I,J-1))
          DIFFf_V =( ALPHA(I,J)*f_VXIXI - 
     &      BETA(I,J)*(f_VXIET + f_VETXI)+ GAMMA(I,J)*f_VETET
     &     +(RNUT(I,J)+Df_V)*( DLX(I,J)*
     &      ( YXI(I,J)*f_VET - YET(I,J)*f_VXI )
     &     + DLY(I,J)*( XET(I,J)*f_VXI 
     &     - XXI(I,J)*f_VET ) ))*RRJ2(I,J)
C
          f_VT(I,J) = f_V(I,J) + DT*( - UM/D*CNVf_V + 1.0/D/D*
     &           DIFFf_V -K_5*f_IIA(I,J))

  605   CONTINUE
  615 ENDDO

      DO 505 I = 1, NX-1
          f_VT(I,NY) = f_VT(I,NY-1)
  505   CONTINUE

C-----------------Z-COMPLEX CALCULATION-------------
      D1 = DZ 

      DO 616 J = 1, NY-1
      DO 606 I = 1, NX-1

CC -- CONVECTION  --
          CALL CNVUP(Z, I, J, CNVZ)
CC -- DIFFUSION  --
          CALL FDIFF(Z,ZXIXI,ZXIET,ZETXI,ZETET,I,J,D1)
C
          ZXI = 0.5D0*(Z(I+1,J) - Z(I-1,J))
          ZET = 0.5D0*(Z(I,J+1) - Z(I,J-1))
          DIFFZ =( ALPHA(I,J)*ZXIXI - 
     &      BETA(I,J)*(ZXIET + ZETXI)+ GAMMA(I,J)*ZETET
     &     +(RNUT(I,J)+DZ)*( DLX(I,J)*
     &      ( YXI(I,J)*ZET - YET(I,J)*ZXI )
     &     + DLY(I,J)*( XET(I,J)*ZXI 
     &     - XXI(I,J)*ZET ) ))*RRJ2(I,J)

          ZT(I,J) = Z(I,J) + DT*( - UM/D*CNVZ + 1.0/D/D*
     &           DIFFZ + 1.0D2*(AP(I,J)/RP_init)*
     &      f_VIIIA(I,J)*f_IXA(I,J)/(1.0D2*(AP(I,J)/RP_init)
     &      +K_A*APC(I,J)) )

  606   CONTINUE
  616 ENDDO

      DO 506 I = 1, NX-1
          ZT(I,NY) = ZT(I,NY-1)
  506   CONTINUE  

C-----------------W-COMPLEX CALCULATION-------------
      D1 = DW

      DO 617 J = 1, NY-1
      DO 607 I = 1, NX-1

CC -- CONVECTION  --
          CALL CNVUP(W, I, J, CNVW)
CC -- DIFFUSION  --
          CALL FDIFF(W,WXIXI,WXIET,WETXI,WETET,I,J,D1)
C
          WXI = 0.5D0*(W(I+1,J) - W(I-1,J))
          WET = 0.5D0*(W(I,J+1) - W(I,J-1))
          DIFFW =( ALPHA(I,J)*WXIXI - 
     &      BETA(I,J)*(WXIET + WETXI)+ GAMMA(I,J)*WETET
     &     +(RNUT(I,J)+DW)*( DLX(I,J)*
     &      ( YXI(I,J)*WET - YET(I,J)*WXI )
     &     + DLY(I,J)*( XET(I,J)*WXI 
     &     - XXI(I,J)*WET ) ))*RRJ2(I,J)

          WT(I,J) = W(I,J) + DT*(  - UM/D*CNVW + 1.0/D/D*
     &           DIFFW + 1.0D2*
     &      (AP(I,J)/RP_init)*
     &      f_VA(I,J)*f_XA(I,J)/( 1.0D2*(AP(I,J)/RP_init)
     &      +K_A*APC(I,J)))

  607   CONTINUE
  617 ENDDO

      DO 507 I = 1, NX-1
          WT(I,NY) = WT(I,NY-1)
  507   CONTINUE

C-----------------Factor Xa CALCULATION-------------
      D1=Df_Xa

      DO 618 J = 1, NY-1
      DO 608 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_XA, I, J, CNVf_XA)
CC -- DIFFUSION  --
          CALL FDIFF(f_XA,f_XAXIXI,f_XAXIET,f_XAETXI,f_XAETET,I,J,D1)
C
          f_XAXI = 0.5D0*(f_XA(I+1,J) - f_XA(I-1,J))
          f_XAET = 0.5D0*(f_XA(I,J+1) - f_XA(I,J-1))
          DIFFf_XA =( ALPHA(I,J)*f_XAXIXI 
     &       -BETA(I,J)*(f_XAXIET + f_XAETXI)
     &     + GAMMA(I,J)*f_XAETET+ (RNUT(I,J)+Df_Xa)*
     &    ( DLX(I,J)*( YXI(I,J)*f_XAET 
     &    - YET(I,J)*f_XAXI )+ DLY(I,J)*( XET(I,J)*f_XAXI
     &     - XXI(I,J)*f_XAET )))*RRJ2(I,J)
C
          f_XAT(I,J) = f_XA(I,J) + DT*( - UM/D*CNVf_XA 
     &      +1.0/D/D*DIFFf_XA +K_10*Z(I,J)
     &      -H_10*f_XA(I,J) -K_TFPI*TFPI(I,J)*f_XA(I,J))

  608   CONTINUE
  618 ENDDO

      DO 508 I = 1, NX-1
          f_XAT(I,NY) = f_XAT(I,NY-1)
  508   CONTINUE


C-----------------Factor X (Stuart-Prower factor, 
C                   autothromboplastin III) CALCULATION-------------
      D1=Df_X

      DO 619 J = 1, NY-1
      DO 609 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_X, I, J, CNVf_X)
CC -- DIFFUSION  --
          CALL FDIFF(f_X,f_XXIXI,f_XXIET,f_XETXI,f_XETET,I,J,D1)
C
          f_XXI = 0.5D0*(f_X(I+1,J) - f_X(I-1,J))
          f_XET = 0.5D0*(f_X(I,J+1) - f_X(I,J-1))
          DIFFf_X =( ALPHA(I,J)*f_XXIXI
     &      -BETA(I,J)*(f_XXIET + f_XETXI)
     &      +GAMMA(I,J)*f_XETET
     &      +(RNUT(I,J)+Df_X)*( DLX(I,J)*( YXI(I,J)*f_XET 
     &      -YET(I,J)*f_XXI )
     &      +DLY(I,J)*( XET(I,J)*f_XXI-XXI(I,J)*f_XET ))
     &      )*RRJ2(I,J)
C
          f_XT(I,J) = f_X(I,J) + DT*(-UM/D*CNVf_X
     &      +1.0/D/D*DIFFf_X-K_10*Z(I,J))

  609   CONTINUE
  619 ENDDO

      DO 509 I = 1, NX-1
          f_XT(I,NY) = f_XT(I,NY-1)
  509   CONTINUE

C-----------------Factor IIa (Thrombin) CALCULATION-------------
      D1=Df_IIa

      DO 630 J = 1, NY-1
      DO 620 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_IIA, I, J, CNVf_IIA)
CC -- DIFFUSION  --
          CALL FDIFF(f_IIA,f_IIAXIXI,f_IIAXIET,
     &       f_IIAETXI,f_IIAETET,I,J,D1)
C
          f_IIAXI = 0.5D0*(f_IIA(I+1,J) - f_IIA(I-1,J))
          f_IIAET = 0.5D0*(f_IIA(I,J+1) - f_IIA(I,J-1))
          DIFFf_IIA =( ALPHA(I,J)*f_IIAXIXI
     &      -BETA(I,J)*(f_IIAXIET + f_IIAETXI)
     &      +GAMMA(I,J)*f_IIAETET+(RNUT(I,J)+Df_IIa)
     &      *( DLX(I,J)*( YXI(I,J)*f_IIAET
     &      -YET(I,J)*f_IIAXI )+ DLY(I,J)*( XET(I,J)*f_IIAXI
     &      -XXI(I,J)*f_IIAET ) ))*RRJ2(I,J)
C
          f_IIAT(I,J) = f_IIA(I,J) + DT*(-UM/D*CNVf_IIA 
     &      +1.0/D/D*DIFFf_IIA +K_2*W(I,J)*f_II(I,J)/
     &      (K_2M + f_II(I,J))-H_2*f_IIA(I,J))

  620   CONTINUE
  630 ENDDO

      DO 520 I = 1, NX-1
          f_IIAT(I,NY) = f_IIAT(I,NY-1) 
     &        + (DT*2*Jf_IIa(I)/DYB(I))*1000*BETACOF
  520   CONTINUE  

C-----------------Factor II (Prothrombin) CALCULATION-------------
      D1=Df_II

      DO 631 J = 1, NY-1
      DO 621 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_II, I, J, CNVf_II)
CC -- DIFFUSION  --
          CALL FDIFF(f_II,f_IIXIXI,f_IIXIET,f_IIETXI,f_IIETET,I,J,D1)
C
          f_IIXI = 0.5D0*(f_II(I+1,J) - f_II(I-1,J))
          f_IIET = 0.5D0*(f_II(I,J+1) - f_II(I,J-1))
          DIFFf_II =( ALPHA(I,J)*f_IIXIXI
     &      -BETA(I,J)*(f_IIXIET + f_IIETXI)
     &      +GAMMA(I,J)*f_IIETET+(RNUT(I,J)+Df_II)
     &     *( DLX(I,J)*( YXI(I,J)*f_IIET
     &      -YET(I,J)*f_IIXI )+ DLY(I,J)*( XET(I,J)*f_IIXI
     &      -XXI(I,J)*f_IIET ) ))*RRJ2(I,J)
C
          f_IIT(I,J) = f_II(I,J) + DT*( - UM/D*CNVf_II 
     &      + 1.0/D/D*DIFFf_II -K_2*W(I,J)*f_II(I,J)/(K_2M + f_II(I,J)))

  621   CONTINUE
  631 ENDDO

      DO 521 I = 1, NX-1
          f_IIT(I,NY) = f_IIT(I,NY-1) - (DT*2*Jf_II(I)/DYB(I))*1000
  521   CONTINUE  

C-----------------Factor Ia (Fibrin) CALCULATION-------------
      D1=Df_Ia

      DO 632 J = 1, NY-1
      DO 622 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_IA, I, J, CNVf_IA)
CC -- DIFFUSION  --
          CALL FDIFF(f_IA,f_IAXIXI,f_IAXIET,f_IAETXI,f_IAETET,I,J,D1)
C
          f_IAXI = 0.5D0*(f_IA(I+1,J) - f_IA(I-1,J))
          f_IAET = 0.5D0*(f_IA(I,J+1) - f_IA(I,J-1))
          DIFFf_IA =( ALPHA(I,J)*f_IAXIXI -BETA(I,J)*
     &     (f_IAXIET + f_IAETXI) + GAMMA(I,J)*f_IAETET
     &     + (RNUT(I,J)+Df_Ia)*( DLX(I,J)*( YXI(I,J)*f_IAET
     &     - YET(I,J)*f_IAXI )+ DLY(I,J)*
     &     (XET(I,J)*f_IAXI - XXI(I,J)*f_IAET)))*RRJ2(I,J)
C
          f_IAT(I,J) = f_IA(I,J) + DT*( - UM/D*CNVf_IA
     &       + 1.0/D/D*DIFFf_IA
     &      +K_1*f_IIA(I,J)*f_I(I,J)/(K_1M + f_I(I,J))
     &      -H_1*PLA(I,J)*f_IA(I,J)/(H_1M + f_IA(I,J)))

  622   CONTINUE
  632 ENDDO

        DO 522 I = 1, NX-1
          f_IAT(I,NY) = f_IAT(I,NY-1) - (DT*2*Jf_IA(I)/DYB(I))*1000
  522   CONTINUE


C-----------------Factor I (Fibrinogen) CALCULATION-------------
      D1=Df_I

      DO 633 J = 1, NY-1
      DO 623 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_I, I, J, CNVf_I)
CC -- DIFFUSION  --
          CALL FDIFF(f_I,f_IXIXI,f_IXIET,f_IETXI,f_IETET,I,J,D1)
C
          f_IXI = 0.5D0*(f_I(I+1,J) - f_I(I-1,J))
          f_IET = 0.5D0*(f_I(I,J+1) - f_I(I,J-1))
          DIFFf_I =( ALPHA(I,J)*f_IXIXI-BETA(I,J)
     &     *(f_IXIET + f_IETXI)+ GAMMA(I,J)*f_IETET
     &   + (RNUT(I,J)+Df_I)*( DLX(I,J)*
     &   ( YXI(I,J)*f_IET - YET(I,J)*f_IXI)
     &   + DLY(I,J)*( XET(I,J)*f_IXI - XXI(I,J)
     &   *f_IET )))*RRJ2(I,J)
C
          f_IT(I,J) = f_I(I,J) + DT*( - UM/D*CNVf_I + 1.0/D/D*
     &      DIFFf_I -K_1*f_IIA(I,J)*f_I(I,J)/(K_1M + f_I(I,J)))

  623   CONTINUE
  633 ENDDO

      DO 523 I = 1, NX-1
          f_IT(I,NY) = f_IT(I,NY-1)
  523   CONTINUE

C-----------------f_actor XIa CALCULATION-------------
      D1=Df_XIa

      DO 634 J = 1, NY-1
      DO 624 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_XIA, I, J, CNVf_XIA)
CC -- DIFFUSION  --
          CALL FDIFF(f_XIA,f_XIAXIXI,f_XIAXIET,
     &     f_XIAETXI,f_XIAETET,I,J,D1)
C
          f_XIAXI = 0.5D0*(f_XIA(I+1,J) - f_XIA(I-1,J))
          f_XIAET = 0.5D0*(f_XIA(I,J+1) - f_XIA(I,J-1))
          DIFFf_XIA =( ALPHA(I,J)*f_XIAXIXI
     &      -BETA(I,J)*(f_XIAXIET + f_XIAETXI)
     &      +GAMMA(I,J)*f_XIAETET+(RNUT(I,J)+Df_XIa)*
     &      ( DLX(I,J)*( YXI(I,J)*f_XIAET
     &      -YET(I,J)*f_XIAXI )+DLY(I,J)*( XET(I,J)*f_XIAXI
     &      -XXI(I,J)*f_XIAET )))*RRJ2(I,J)
C
          f_XIAT(I,J) = f_XIA(I,J) + DT*(-UM/D*CNVf_XIA+1.0/D/D*DIFFf_XIA
     &      +K_11*f_IIA(I,J)-H_11*f_XIA(I,J))

  624   CONTINUE
  634 ENDDO

      DO 524 I = 1, NX-1
          f_XIAT(I,NY) = f_XIAT(I,NY-1)
  524   CONTINUE


C-------Factor XI (Plasma thromboplastin antecedent, PTA) CALCULATION--------
      D1=Df_XI

      DO 635 J = 1, NY-1
      DO 625 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(f_XI, I, J, CNVf_XI)
CC -- DIFFUSION  --
          CALL FDIFF(f_XI,f_XIXIXI,f_XIXIET,f_XIETXI,f_XIETET,I,J,D1)
C
          f_XIXI = 0.5D0*(f_XI(I+1,J) - f_XI(I-1,J))
          f_XIET = 0.5D0*(f_XI(I,J+1) - f_XI(I,J-1))
          DIFFf_XI =( ALPHA(I,J)*f_XIXIXI-BETA(I,J)*(f_XIXIET+f_XIETXI)
     &             + GAMMA(I,J)*f_XIETET
     &     + (RNUT(I,J)+Df_XI)*( DLX(I,J)*( YXI(I,J)*f_XIET
     &     - YET(I,J)*f_XIXI )
     &    + DLY(I,J)*( XET(I,J)*f_XIXI - XXI(I,J)*f_XIET)))*RRJ2(I,J)
C
          f_XIT(I,J) = f_XI(I,J) + DT*( - UM/D*CNVf_XI + 1.0/D/D*
     &           DIFFf_XI -K_11*f_IIA(I,J) )

  625   CONTINUE
  635 ENDDO

      DO 525 I = 1, NX-1
          f_XIT(I,NY) = f_XIT(I,NY-1)
  525   CONTINUE  


C-----------------Antithrombin III CALCULATION-------------
      D1=DATIII

      DO 636 J = 1, NY-1
      DO 626 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(ATIII, I, J, CNVATIII)
CC -- DIFFUSION  --
          CALL FDIFF(ATIII,ATIIIXIXI,ATIIIXIET,
     &    ATIIIETXI,ATIIIETET,I,J,D1)
C
          ATIIIXI = 0.5D0*(ATIII(I+1,J) - ATIII(I-1,J))
          ATIIIET = 0.5D0*(ATIII(I,J+1) - ATIII(I,J-1))
          DIFFATIII =( ALPHA(I,J)*ATIIIXIXI
     &      -BETA(I,J)*(ATIIIXIET + ATIIIETXI)
     &      +GAMMA(I,J)*ATIIIETET
     &      +(RNUT(I,J)+DATIII)*( DLX(I,J)*( YXI(I,J)*ATIIIET
     &      -YET(I,J)*ATIIIXI )+DLY(I,J)*( XET(I,J)*ATIIIXI
     &      -XXI(I,J)*ATIIIET )))*RRJ2(I,J)
C
          ATIIIT(I,J) = ATIII(I,J) + DT*(-UM/D*CNVATIII
     &      +1.0/D/D*DIFFATIII-H_9*f_IXA(I,J)
     &      -H_10*f_XA(I,J)-H_2*f_IIA(I,J))

  626   CONTINUE
  636 ENDDO

      DO 526 I = 1, NX-1
          ATIIIT(I,NY) = ATIIIT(I,NY-1)
  526   CONTINUE

C-----------------TFPI CALCULATION-------------
      D1=DTFPI

      DO 637 J = 1, NY-1
      DO 627 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(TFPI, I, J, CNVTFPI)
CC -- DIFFUSION  --
          CALL FDIFF(TFPI,TFPIXIXI,TFPIXIET,TFPIETXI,TFPIETET,I,J,D1)
C
          TFPIXI = 0.5D0*(TFPI(I+1,J) - TFPI(I-1,J))
          TFPIET = 0.5D0*(TFPI(I,J+1) - TFPI(I,J-1))
          DIFFTFPI =( ALPHA(I,J)*TFPIXIXI
     &      -BETA(I,J)*(TFPIXIET+TFPIETXI)
     &      +GAMMA(I,J)*TFPIETET+(RNUT(I,J)+DTFPI)*
     &      (DLX(I,J)*(YXI(I,J)*TFPIET
     &      -YET(I,J)*TFPIXI)+DLY(I,J)
     &      *(XET(I,J)*TFPIXI
     &      -XXI(I,J)*TFPIET)))*RRJ2(I,J)
C
          TFPIT(I,J) = TFPI(I,J)+DT*(-UM/D*CNVTFPI
     &      +1.0/D/D*DIFFTFPI-K_TFPI*TFPI(I,J)*f_XA(I,J))

  627   CONTINUE
  637 ENDDO

      DO 527 I = 1, NX-1
          TFPIT(I,NY) = TFPIT(I,NY-1)
  527   CONTINUE


C-----------------APC CALCULATION-------------
      D1=DAPC

      DO 638 J = 1, NY-1
      DO 628 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(APC, I, J, CNVAPC)
CC -- DIFFUSION  --
          CALL FDIFF(APC,APCXIXI,APCXIET,APCETXI,APCETET,I,J,D1)
C
          APCXI = 0.5D0*(APC(I+1,J) - APC(I-1,J))
          APCET = 0.5D0*(APC(I,J+1) - APC(I,J-1))
          DIFFAPC =(ALPHA(I,J)*APCXIXI
     &      -BETA(I,J)*(APCXIET +APCETXI)
     &      +GAMMA(I,J)*APCETET+(RNUT(I,J)+DAPC)*( DLX(I,J)*
     &      ( YXI(I,J)*APCET - YET(I,J)*APCXI )
     &      +DLY(I,J)*( XET(I,J)*APCXI
     &      -XXI(I,J)*APCET)))*RRJ2(I,J)
C
          APCT(I,J) = APC(I,J) + DT*(-UM/D*CNVAPC 
     &  +1.0/D/D*DIFFAPC +KAPC*f_IIA(I,J)-HAPC*APC(I,J))

  628   CONTINUE
  638 ENDDO

      DO 528 I = 1, NX-1
          APCT(I,NY) = APCT(I,NY-1)
  528   CONTINUE

C-----------------PC CALCULATION-------------
      D1=DPC

      DO 639 J = 1, NY-1
      DO 629 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(PC, I, J, CNVPC)
CC -- DIFFUSION  --
          CALL FDIFF(PC,PCXIXI,PCXIET,PCETXI,PCETET,I,J,D1)
C
          PCXI = 0.5D0*(PC(I+1,J) - PC(I-1,J))
          PCET = 0.5D0*(PC(I,J+1) - PC(I,J-1))
          DIFFPC =( ALPHA(I,J)*PCXIXI - BETA(I,J)*(PCXIET + PCETXI)
     &             + GAMMA(I,J)*PCETET
     &     + (RNUT(I,J)+DPC)*( DLX(I,J)*( YXI(I,J)*PCET
     &      - YET(I,J)*PCXI )+ DLY(I,J)*( XET(I,J)*PCXI 
     &      - XXI(I,J)*PCET )))*RRJ2(I,J)
C
          PCT(I,J) = PC(I,J) + DT*( - UM/D*CNVPC + 1.0/D/D*
     &           DIFFPC -KAPC*f_IIA(I,J))

 629   CONTINUE
 639  ENDDO

      DO 529 I = 1, NX-1
          PCT(I,NY) = PCT(I,NY-1)
  529   CONTINUE  

C-----------------rL1AT CALCULATION-------------
      D1=DrL1AT

      DO 650 J = 1, NY-1
      DO 640 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(rL1AT, I, J, CNVrL1AT)
CC -- DIFFUSION  --
          CALL FDIFF(rL1AT,rL1ATXIXI,rL1ATXIET,
     &            rL1ATETXI,rL1ATETET,I,J,D1)
C
          rL1ATXI = 0.5D0*(rL1AT(I+1,J) - rL1AT(I-1,J))
          rL1ATET = 0.5D0*(rL1AT(I,J+1) - rL1AT(I,J-1))
          DIFFrL1AT =(ALPHA(I,J)*rL1ATXIXI-BETA(I,J)*
     &      (rL1ATXIET +rL1ATETXI)
     &      +GAMMA(I,J)*rL1ATETET
     &      +(RNUT(I,J)+DrL1AT)*( DLX(I,J)*( YXI(I,J)*rL1ATET
     &      -YET(I,J)*rL1ATXI )
     &      +DLY(I,J)*( XET(I,J)*rL1ATXI
     &      -XXI(I,J)*rL1ATET)))*RRJ2(I,J)
C
          rL1ATT(I,J) = rL1AT(I,J) + DT*( - UM/D*CNVrL1AT + 1.0/D/D*
     &           DIFFrL1AT -HAPC*APC(I,J))

  640   CONTINUE
  650 ENDDO

      DO 540 I = 1, NX-1
          rL1ATT(I,NY) = rL1ATT(I,NY-1)
  540   CONTINUE

C-----------------tPA CALCULATION-------------
      D1=DTPA

      DO 651 J = 1, NY-1
      DO 641 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(TPA, I, J, CNVTPA)
CC -- DIFFUSION  --
          CALL FDIFF(TPA,TPAXIXI,TPAXIET,TPAETXI,TPAETET,I,J,D1)
C
          TPAXI = 0.5D0*(TPA(I+1,J) - TPA(I-1,J))
          TPAET = 0.5D0*(TPA(I,J+1) - TPA(I,J-1))
          DIFFTPA =( ALPHA(I,J)*TPAXIXI
     &      -BETA(I,J)*(TPAXIET
     &      +TPAETXI)
     &      +GAMMA(I,J)*TPAETET
     &      +(RNUT(I,J)+DTPA)*( DLX(I,J)*( YXI(I,J)*TPAET
     &      -YET(I,J)*TPAXI)
     &      +DLY(I,J)*( XET(I,J)*TPAXI
     &      -XXI(I,J)*TPAET)))*RRJ2(I,J)
C
          TPAT(I,J) = TPA(I,J) + DT*( - UM/D*CNVTPA + 1.0/D/D*
     &           DIFFTPA +0)

  641   CONTINUE
  651 ENDDO

      DO 541 I = 1, NX-1
          TPAT(I,NY) = TPAT(I,NY-1)
  541   CONTINUE



C-----------------PLA CALCULATION-------------
      D1=DPLA

      DO 652 J = 1, NY-1
      DO 642 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(PLA, I, J, CNVPLA)
CC -- DIFFUSION  --
          CALL FDIFF(PLA,PLAXIXI,PLAXIET,PLAETXI,PLAETET,I,J,D1)
C
          PLAXI = 0.5D0*(PLA(I+1,J) - PLA(I-1,J))
          PLAET = 0.5D0*(PLA(I,J+1) - PLA(I,J-1))
          DIFFPLA =( ALPHA(I,J)*PLAXIXI
     &      -BETA(I,J)*(PLAXIET + PLAETXI)
     &      + GAMMA(I,J)*PLAETET
     &      +(RNUT(I,J)+DPLA)*( DLX(I,J)*( YXI(I,J)*PLAET
     &      -YET(I,J)*PLAXI )+ DLY(I,J)*( XET(I,J)*PLAXI
     &      -XXI(I,J)*PLAET )))*RRJ2(I,J)
C
          PLAT(I,J) = PLA(I,J) + DT*( - UM/D*CNVPLA + 1.0/D/D*
     &           DIFFPLA +KTPA_IA_PLA*TPA(I,J)
     &      -HPLA*PLA(I,J)*rL2AP(I,J))

  642   CONTINUE
  652 ENDDO

      DO 542 I = 1, NX-1
          PLAT(I,NY) = PLAT(I,NY-1)
  542   CONTINUE

C-----------------PLS CALCULATION-------------
      D1=DPLS

      DO 653 J = 1, NY-1
      DO 643 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(PLS, I, J, CNVPLS)
CC -- DIFFUSION  --
          CALL FDIFF(PLS,PLSXIXI,PLSXIET,PLSETXI,PLSETET,I,J,D1)
C
          PLSXI = 0.5D0*(PLS(I+1,J) - PLS(I-1,J))
          PLSET = 0.5D0*(PLS(I,J+1) - PLS(I,J-1))
          DIFFPLS =( ALPHA(I,J)*PLSXIXI-
     &      BETA(I,J)*(PLSXIET + PLSETXI)
     &      + GAMMA(I,J)*PLSETET
     &      + (RNUT(I,J)+DPLS)*( DLX(I,J)*( YXI(I,J)*PLSET
     &      -YET(I,J)*PLSXI )+DLY(I,J)*( XET(I,J)*PLSXI
     &      -XXI(I,J)*PLSET ) )
     &           )*RRJ2(I,J)
C
          PLST(I,J) = PLS(I,J) + DT*( - UM/D*CNVPLS + 1.0/D/D*
     &           DIFFPLS -KTPA_IA_PLA*TPA(I,J))

  643   CONTINUE
  653 ENDDO

      DO 543 I = 1, NX-1
          PLST(I,NY) = PLST(I,NY-1)
  543   CONTINUE

C-----------------rL2AP CALCULATION-------------
      D1=DrL2AP

      DO 654 J = 1, NY-1
      DO 644 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(rL2AP, I, J, CNVrL2AP)
CC -- DIFFUSION  --
          CALL FDIFF(rL2AP,rL2APXIXI,rL2APXIET,
     &              rL2APETXI,rL2APETET,I,J,D1)
C
          rL2APXI = 0.5D0*(rL2AP(I+1,J) - rL2AP(I-1,J))
          rL2APET = 0.5D0*(rL2AP(I,J+1) - rL2AP(I,J-1))
          DIFFrL2AP =( ALPHA(I,J)*rL2APXIXI
     &      -BETA(I,J)*(rL2APXIET + rL2APETXI)
     &      +GAMMA(I,J)*rL2APETET
     &      +(RNUT(I,J)+DrL2AP)*( DLX(I,J)*
     &      ( YXI(I,J)*rL2APET - YET(I,J)*rL2APXI )
     &      +DLY(I,J)*(XET(I,J)*rL2APXI - XXI(I,J)*rL2APET))
     &      )*RRJ2(I,J)
C
          rL2APT(I,J) = rL2AP(I,J) + DT*( - UM/D*CNVrL2AP + 1.0/D/D*
     &           DIFFrL2AP -HPLA*PLA(I,J)*rL2AP(I,J))

  644   CONTINUE
  654 ENDDO

      DO 544 I = 1, NX-1
          rL2APT(I,NY) = rL2APT(I,NY-1)
  544   ENDDO

C-----------------Resting PLT CALCULATION-------------
C
      D1=DRP

      DO 655 J = 1, NY-1
        DO 645 I = 1, NX-1

CC -- CONVECTION  --

          CALL CNVUP(RP, I, J, CNVRP)
CC -- DIFFUSION  --
          CALL FDIFF(RP,RPXIXI,RPXIET,RPETXI,RPETET,I,J,D1)
C
          RPXI = 0.5D0*(RP(I+1,J) - RP(I-1,J))
          RPET = 0.5D0*(RP(I,J+1) - RP(I,J-1))
          DIFFRP =( ALPHA(I,J)*RPXIXI -   BETA(I,J)*(RPXIET + RPETXI)
     &             + GAMMA(I,J)*RPETET
     &    + (RNUT(I,J)+DRP)*( DLX(I,J)*( YXI(I,J)*RPET - YET(I,J)*RPXI )
     &             + DLY(I,J)*( XET(I,J)*RPXI - XXI(I,J)*RPET ) )
     &           )*RRJ2(I,J)
C
          RPT(I,J) = RP(I,J) + DT*( - UM/D*CNVRP +  1.0/D/D*
     &           DIFFRP-KPA(I,J)*RP(I,J))


  645  CONTINUE
  655 ENDDO

        DO 545 I = 1, NX-1
          RPT(I,NY) = RPT(I,NY-1) - DT*2*JR(I)/DYB(I)
  545  CONTINUE

C
C-------------Activated PLT CALCULATION---------
      D1=DAP

      DO 656 J = 1, NY-1
      DO 646 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(AP, I, J, CNVAP)
CC -- DIFFUSION  --
          CALL FDIFF(AP,APXIXI,APXIET,APETXI,APETET,I,J,D1)
C
          APXI = 0.5D0*(AP(I+1,J) - AP(I-1,J))
          APET = 0.5D0*(AP(I,J+1) - AP(I,J-1))
          DIFFAP =( ALPHA(I,J)*APXIXI -   BETA(I,J)*(APXIET + APETXI)
     &             + GAMMA(I,J)*APETET
     &    + (RNUT(I,J)+Dap)*( DLX(I,J)*( YXI(I,J)*APET - YET(I,J)*APXI )
     &             + DLY(I,J)*( XET(I,J)*APXI - XXI(I,J)*APET ) )
     &           )*RRJ2(I,J)
C
          APT(I,J) = AP(I,J) + DT*( - UM/D*CNVAP + 1.0/D/D*
     &           DIFFAP +KPA(I,J)*RP(I,J))

  646   CONTINUE
  656 ENDDO

      DO 546 I = 1, NX-1
          APT(I,NY) = APT(I,NY-1) - DT*2*JA(I)/DYB(I)
  546   CONTINUE

C-------------PLT-released agonists CALCULATION----------
      D1=Dapr

      DO 657 J = 1, NY-1
      DO 647 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(APR, I, J, CNVAPR)
CC -- DIFFUSION  --
          CALL FDIFF(APR,APRXIXI,APRXIET,APRETXI,APRETET,I,J,D1)
C
          APRXI = 0.5D0*(APR(I+1,J) - APR(I-1,J))
          APRET = 0.5D0*(APR(I,J+1) - APR(I,J-1))
          DIFFAPR = ( ALPHA(I,J)*APRXIXI -   BETA(I,J)*(APRXIET 
     &             + APRETXI) + GAMMA(I,J)*APRETET + (RNUT(I,J)+Dapr)
     &   *( DLX(I,J)*( YXI(I,J)*APRET - YET(I,J)*APRXI )
     &     + DLY(I,J)*( XET(I,J)*APRXI - XXI(I,J)*APRET ) ))*RRJ2(I,J)

          APRT(I,J) = APR(I,J) + DT*( -UM/D*CNVAPR + 1.0/D/D*
     &          DIFFAPR +Lambdaj*KPA(I,J)*RP(I,J)-K1_J*APR(I,J))
  647   CONTINUE
  657 ENDDO 

      DO 547 I = 1, NX-1
          APRT(I,NY) = APRT(I,NY-1) + DT*2*JPR(I)/DYB(I)
  547   CONTINUE

C----------------- PLT-synthesized agonist CALCULATION-------------
C
      D1=Daps

      DO 658 J = 1, NY-1
      DO 648 I = 1, NX-1
CC -- CONVECTION  --
          CALL CNVUP(APS, I, J, CNVAPS)
CC -- DIFFUSION  --
          CALL FDIFF(APS,APSXIXI,APSXIET,APSETXI,APSETET,I,J,D1)
C
          APSXI = 0.5D0*(APS(I+1,J) - APS(I-1,J))
          APSET = 0.5D0*(APS(I,J+1) - APS(I,J-1))
          DIFFAPS =( ALPHA(I,J)*APSXIXI -   BETA(I,J)*(APSXIET 
     &             + APSETXI) + GAMMA(I,J)*APSETET+ (RNUT(I,J)+Daps)
     &   *( DLX(I,J)*( YXI(I,J)*APSET - YET(I,J)*APSXI )
     &             + DLY(I,J)*( XET(I,J)*APSXI - XXI(I,J)*APSET ) )
     &           )*RRJ2(I,J)
C
          APST(I,J) = APS(I,J) + DT*( - UM/D*CNVAPS + 1.0/D/D*
     &                    DIFFAPS +spj*AP(I,J)-k1_j*APS(I,J))       
  648  CONTINUE
  658 ENDDO

      DO 548 I = 1, NX-1
          APST(I,NY) = APST(I,NY-1) + DT*2*JPS(I)/DYB(I)
  548  CONTINUE

C------------------------------------------------
C
      TIM=TIM+DT
      KDEP=19.3888692D-4*TIM-9.3982062D-7*(TIM**2)  

      ! Z_bunshi = 1.0D2*(AP(150,5)/RP_init)*f_VIIIA(150,5)*f_IXA(150,5)
      ! Z_bunbo_1 = 1.0D2*(AP(150,5)/RP_init)
      ! Z_bunbo_2 = K_A*APC(150,5)

       WRITE(6,*) TIM, ZT(150,5),WT(150,5),APCT(150,5),f_IIAT(150,5)
      !  WRITE(6,*) TIM, ZT(150,5),Z_bunshi,Z_bunbo_1,Z_bunbo_2,APC(150,5)


       DO 201 J=1,NY
       DO 202 I=1,NX-1
         RP(I,J)=RPT(I,J)
         AP(I,J)=APT(I,J)
         APR(I,J)=APRT(I,J)
         APS(I,J)=APST(I,J)
         f_I(I,J)=f_IT(I,J)
         f_IA(I,J)=f_IAT(I,J)
         f_II(I,J)=f_IIT(I,J)
         f_IIA(I,J)=f_IIAT(I,J)
         f_V(I,J)=f_VT(I,J)
         f_VA(I,J)=f_VAT(I,J)
         f_VIII(I,J)=f_VIIIT(I,J)
         f_VIIIA(I,J)=f_VIIIAT(I,J)
         f_IX(I,J)=f_IXT(I,J)
         f_IXA(I,J)=f_IXAT(I,J)
         f_X(I,J)=f_XT(I,J)
         f_XA(I,J)=f_XAT(I,J)
         f_XI(I,J)=f_XIT(I,J)
         f_XIA(I,J)=f_XIAT(I,J)
         ATIII(I,J)=ATIIIT(I,J)
         TFPI(I,J)=TFPIT(I,J)
         PC(I,J)=PCT(I,J)
         APC(I,J)=APCT(I,J)
         rL1AT(I,J)=rL1ATT(I,J)
         TPA(I,J)=TPAT(I,J)
         PLS(I,J)=PLST(I,J)
         PLA(I,J)=PLAT(I,J)
         rL2AP(I,J)=rL2APT(I,J)
         W(I,J)=WT(I,J)
         Z(I,J)=ZT(I,J)
  202  ENDDO
  201  ENDDO
C
      DO 113 I = 1, NX-1
        RP(I,0)=RP(I,1)
        AP(I,0)=AP(I,1)
        APR(I,0)=APR(I,1)
        APS(I,0)=APS(I,1)
        f_I(I,0)=f_I(I,1)
        f_IA(I,0)=f_IA(I,1)
        f_II(I,0)=f_II(I,1)
        f_IIA(I,0)=f_IIA(I,1)
        f_V(I,0)=f_V(I,1)
        f_VA(I,0)=f_VA(I,1)
        f_VIII(I,0)=f_VIII(I,1)
        f_VIIIA(I,0)=f_VIIIA(I,1)
        f_IX(I,0)=f_IX(I,1)
        f_IXA(I,0)=f_IXA(I,1)
        f_X(I,0)=f_X(I,1)
        f_XA(I,0)=f_XA(I,1)
        f_XI(I,0)=f_XI(I,1)
        f_XIA(I,0)=f_XIA(I,1)
        ATIII(I,0)=ATIII(I,1)
        TFPI(I,0)=TFPI(I,1)
        PC(I,0)=PC(I,1)
        APC(I,0)=APC(I,1)
        rL1AT(I,0)=rL1AT(I,1)
        TPA(I,0)=TPA(I,1)
        PLS(I,0)=PLS(I,1)
        PLA(I,0)=PLA(I,1)
        rL2AP(I,0)=rL2AP(I,1)
        W(I,0)=W(I,1)
        Z(I,0)=Z(I,1)

        RP(I,-1)=RP(I,1)
        AP(I,-1)=AP(I,1)
        APR(I,-1)=APR(I,1)
        APS(I,-1)=APS(I,1)
         f_I(I,-1)=f_I(I,1)
         f_IA(I,-1)=f_IA(I,1)
         f_II(I,-1)=f_II(I,1)
         f_IIA(I,-1)=f_IIA(I,1)
         f_V(I,-1)=f_V(I,1)
         f_VA(I,-1)=f_VA(I,1)
         f_VIII(I,-1)=f_VIII(I,1)
         f_VIIIA(I,-1)=f_VIIIA(I,1)
         f_IX(I,-1)=f_IX(I,1)
         f_IXA(I,-1)=f_IXA(I,1)
         f_X(I,-1)=f_X(I,1)
         f_XA(I,-1)=f_XA(I,1)
         f_XI(I,-1)=f_XI(I,1)
         f_XIA(I,-1)=f_XIA(I,1)
         ATIII(I,-1)=ATIII(I,1)
         TFPI(I,-1)=TFPI(I,1)
         PC(I,-1)=PC(I,1)
         APC(I,-1)=APC(I,1)
         rL1AT(I,-1)=rL1AT(I,1)
         TPA(I,-1)=TPA(I,1)
         PLS(I,-1)=PLS(I,1)
         PLA(I,-1)=PLA(I,1)
         rL2AP(I,-1)=rL2AP(I,1)
         W(I,-1)=W(I,1)
         Z(I,-1)=Z(I,1)
  113 CONTINUE

       DO 256 J=-1,NY

       RP(0,J)=RP(NX-1,J)
       AP(0,J)=AP(NX-1,J)
       APR(0,J)=APR(NX-1,J)
       APS(0,J)=APS(NX-1,J)
       f_I(0,J)=f_I(NX-1,J)
       f_IA(0,J)=f_IA(NX-1,J)
       f_II(0,J)=f_II(NX-1,J)
       f_IIA(0,J)=f_IIA(NX-1,J)
       f_V(0,J)=f_V(NX-1,J)
       f_VA(0,J)=f_VA(NX-1,J)
       f_VIII(0,J)=f_VIII(NX-1,J)
       f_VIIIA(0,J)=f_VIIIA(NX-1,J)
       f_IX(0,J)=f_IX(NX-1,J)
       f_IXA(0,J)=f_IXA(NX-1,J)
       f_X(0,J)=f_X(NX-1,J)
       f_XA(0,J)=f_XA(NX-1,J)
       f_XI(0,J)=f_XI(NX-1,J)
       f_XIA(0,J)=f_XIA(NX-1,J)
       ATIII(0,J)=ATIII(NX-1,J)
       TFPI(0,J)=TFPI(NX-1,J)
       PC(0,J)=PC(NX-1,J)
       APC(0,J)=APC(NX-1,J)
       rL1AT(0,J)=rL1AT(NX-1,J)
       TPA(0,J)=TPA(NX-1,J)
       PLS(0,J)=PLS(NX-1,J)
       PLA(0,J)=PLA(NX-1,J)
       rL2AP(0,J)=rL2AP(NX-1,J)
       W(0,J)=W(NX-1,J)
       Z(0,J)=Z(NX-1,J)

       RP(-1,J)=RP(NX-1,J)
       AP(-1,J)=AP(NX-1,J)
       APR(-1,J)=APR(NX-1,J)
       APS(-1,J)=APS(NX-1,J)
       f_I(-1,J)=f_I(NX-1,J)
       f_IA(-1,J)=f_IA(NX-1,J)
       f_II(-1,J)=f_II(NX-1,J)
       f_IIA(-1,J)=f_IIA(NX-1,J)
       f_V(-1,J)=f_V(NX-1,J)
       f_VA(-1,J)=f_VA(NX-1,J)
       f_VIII(-1,J)=f_VIII(NX-1,J)
       f_VIIIA(-1,J)=f_VIIIA(NX-1,J)
       f_IX(-1,J)=f_IX(NX-1,J)
       f_IXA(-1,J)=f_IXA(NX-1,J)
       f_X(-1,J)=f_X(NX-1,J)
       f_XA(-1,J)=f_XA(NX-1,J)
       f_XI(-1,J)=f_XI(NX-1,J)
       f_XIA(-1,J)=f_XIA(NX-1,J)
       ATIII(-1,J)=ATIII(NX-1,J)
       TFPI(-1,J)=TFPI(NX-1,J)
       PC(-1,J)=PC(NX-1,J)
       APC(-1,J)=APC(NX-1,J)
       rL1AT(-1,J)=rL1AT(NX-1,J)
       TPA(-1,J)=TPA(NX-1,J)
       PLS(-1,J)=PLS(NX-1,J)
       PLA(-1,J)=PLA(NX-1,J)
       rL2AP(-1,J)=rL2AP(NX-1,J)
       W(-1,J)=W(NX-1,J)
       Z(-1,J)=Z(NX-1,J)

        RP(NX,J)=RP(NX-1,J)
        AP(NX,J)=AP(NX-1,J)
        APR(NX,J)=APR(NX-1,J)
        APS(NX,J)=APS(NX-1,J)
        f_I(NX,J)=f_I(NX-1,J)
        f_IA(NX,J)=f_IA(NX-1,J)
        f_II(NX,J)=f_II(NX-1,J)
        f_IIA(NX,J)=f_IIA(NX-1,J)
        f_V(NX,J)=f_V(NX-1,J)
        f_VA(NX,J)=f_VA(NX-1,J)
        f_VIII(NX,J)=f_VIII(NX-1,J)
        f_VIIIA(NX,J)=f_VIIIA(NX-1,J)
        f_IX(NX,J)=f_IX(NX-1,J)
        f_IXA(NX,J)=f_IXA(NX-1,J)
        f_X(NX,J)=f_X(NX-1,J)
        f_XA(NX,J)=f_XA(NX-1,J)
        f_XI(NX,J)=f_XI(NX-1,J)
        f_XIA(NX,J)=f_XIA(NX-1,J)
        ATIII(NX,J)=ATIII(NX-1,J)
        TFPI(NX,J)=TFPI(NX-1,J)
        PC(NX,J)=PC(NX-1,J)
        APC(NX,J)=APC(NX-1,J)
        rL1AT(NX,J)=rL1AT(NX-1,J)
        TPA(NX,J)=TPA(NX-1,J)
        PLS(NX,J)=PLS(NX-1,J)
        PLA(NX,J)=PLA(NX-1,J)
        rL2AP(NX,J)=rL2AP(NX-1,J)
        W(NX,J)=W(NX-1,J)
        Z(NX,J)=Z(NX-1,J)


        RP(NX+1,J)=RP(NX-1,J)
        AP(NX+1,J)=AP(NX-1,J)
        APR(NX+1,J)=APR(NX-1,J)
        APS(NX+1,J)=APS(NX-1,J)
        f_I(NX+1,J)=f_I(NX-1,J)
        f_IA(NX+1,J)=f_IA(NX-1,J)
        f_II(NX+1,J)=f_II(NX-1,J)
        f_IIA(NX+1,J)=f_IIA(NX-1,J)
        f_V(NX+1,J)=f_V(NX-1,J)
        f_VA(NX+1,J)=f_VA(NX-1,J)
        f_VIII(NX+1,J)=f_VIII(NX-1,J)
        f_VIIIA(NX+1,J)=f_VIIIA(NX-1,J)
        f_IX(NX+1,J)=f_IX(NX-1,J)
        f_IXA(NX+1,J)=f_IXA(NX-1,J)
        f_X(NX+1,J)=f_X(NX-1,J)
        f_XA(NX+1,J)=f_XA(NX-1,J)
        f_XI(NX+1,J)=f_XI(NX-1,J)
        f_XIA(NX+1,J)=f_XIA(NX-1,J)
        ATIII(NX+1,J)=ATIII(NX-1,J)
        TFPI(NX+1,J)=TFPI(NX-1,J)
        PC(NX+1,J)=PC(NX-1,J)
        APC(NX+1,J)=APC(NX-1,J)
        rL1AT(NX+1,J)=rL1AT(NX-1,J)
        TPA(NX+1,J)=TPA(NX-1,J)
        PLS(NX+1,J)=PLS(NX-1,J)
        PLA(NX+1,J)=PLA(NX-1,J)
        rL2AP(NX+1,J)=rL2AP(NX-1,J)
        W(NX+1,J)=W(NX-1,J)
        Z(NX+1,J)=Z(NX-1,J)

  256   ENDDO

      DO 214 I = -1, NX+1
        RP(I,NY+1)=RP(I,NY)
        AP(I,NY+1)=AP(I,NY)
        APR(I,NY+1)=APR(I,NY)
        APS(I,NY+1)=APS(I,NY)
        f_I(I,NY+1)=f_I(I,NY)
        f_IA(I,NY+1)=f_IA(I,NY)
        f_II(I,NY+1)=f_II(I,NY)
        f_IIA(I,NY+1)=f_IIA(I,NY)
        f_V(I,NY+1)=f_V(I,NY)
        f_VA(I,NY+1)=f_VA(I,NY)
        f_VIII(I,NY+1)=f_VIII(I,NY)
        f_VIIIA(I,NY+1)=f_VIIIA(I,NY)
        f_IX(I,NY+1)=f_IX(I,NY)
        f_IXA(I,NY+1)=f_IXA(I,NY)
        f_X(I,NY+1)=f_X(I,NY)
        f_XA(I,NY+1)=f_XA(I,NY)
        f_XI(I,NY+1)=f_XI(I,NY)
        f_XIA(I,NY+1)=f_XIA(I,NY)
        ATIII(I,NY+1)=ATIII(I,NY)
        TFPI(I,NY+1)=TFPI(I,NY)
        PC(I,NY+1)=PC(I,NY)
        APC(I,NY+1)=APC(I,NY)
        rL1AT(I,NY+1)=rL1AT(I,NY)
        TPA(I,NY+1)=TPA(I,NY)
        PLS(I,NY+1)=PLS(I,NY)
        PLA(I,NY+1)=PLA(I,NY)
        rL2AP(I,NY+1)=rL2AP(I,NY)
        W(I,NY+1)=W(I,NY)
        Z(I,NY+1)=Z(I,NY)
  214  CONTINUE

      DO 350 I=0, NX
          JR(I) = S(I)*KRS*RP(I,NY)
          JA(I) = (S(I)*KAS+(MAS(I)/7.0D6)*KAA)*AP(I,NY)
          JPR(I) = Lambdaj*(THETA*S(I)*KRS*RP(I,NY))
          JPS(I) = MAT(I)*Spj
          Jf_II(I) = BETACOF*(MAT(I)*(3.22D-8)+MR(I)*FEIRT)*f_II(I,NY)/1000
          Jf_IIa(I) = (MAT(I)*(3.22D-8)+MR(I)*FEIRT)*f_II(I,NY)/1000
          Jf_Ia(I) = ((KDEP*f_IA(I,NY)/1000)/(RHOS-DELTAFB*KDEP*
     &           f_IA(I,NY)/1000))*(DELTAPLT*JR(I)+DELTAPLT*JA(I))
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
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME2='AP'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME3='APR'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME4='APS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'


       NAME5='I'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME6='Ia'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME7='II'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME8='IIa'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'


       NAME9='V'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME10='Va'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME11='VIII'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME12='VIIIa'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME13='IX'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME14='IXa'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME15='X'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME16='Xa'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME17='XI'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME18='XIa'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME19='ATIII'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME20='TFPI'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME21='PC'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME22='APC'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME23='rL1AT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME24='TPA'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME25='PLS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME26='PLA'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME27='L2AP'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME28='W'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME29='Z'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'

       NAME30='KPA'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME31='MTO'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME32='M'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME33='MAS'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME34='MR'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME35='MAT'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'

       NAME36='RPB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'
       NAME37='APB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'
       NAME38='APRB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'
       NAME39='APSB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.TXT'
      NAME40='IB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME41='IaB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME42='IIB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME43='IIaB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME44='VB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME45='VaB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME46='VIIIB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME47='VIIIaB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME48='IXB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME49='IXaB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME50='XB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME51='XaB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME52='XIB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME53='XIaB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME54='ATIIIB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME55='TFPIB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME56='PCB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME57='APCB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME58='rL1ATB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME59='TPAB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME60='PLSB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME61='PLAB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME62='L2APB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME63='WB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)
     &           //CHAR(KL4+48)//CHAR(KL6+48)//'.plt'
       NAME64='ZB'//CHAR(K12+48)
     &           //CHAR(K10+48)
     &           //CHAR(KK8+48)
     &           //CHAR(KK6+48)
     &           //CHAR(KK4+48)
     &           //CHAR(KK2+48)
     &           //CHAR(KK0+48)//CHAR(KL2+48)



      OPEN(2000,FILE=NAME1,STATUS='UNKNOWN')
      OPEN(2010,FILE=NAME2,STATUS='UNKNOWN')
      OPEN(2020,FILE=NAME3,STATUS='UNKNOWN')
      OPEN(2030,FILE=NAME4,STATUS='UNKNOWN')
      OPEN(2040,FILE=NAME5,STATUS='UNKNOWN')
      OPEN(2050,FILE=NAME6,STATUS='UNKNOWN')
      OPEN(2060,FILE=NAME7,STATUS='UNKNOWN')
      OPEN(2070,FILE=NAME8,STATUS='UNKNOWN')
      OPEN(2080,FILE=NAME9,STATUS='UNKNOWN')
      OPEN(2090,FILE=NAME10,STATUS='UNKNOWN')
      OPEN(2100,FILE=NAME11,STATUS='UNKNOWN')
      OPEN(2110,FILE=NAME12,STATUS='UNKNOWN')
      OPEN(2120,FILE=NAME13,STATUS='UNKNOWN')
      OPEN(2130,FILE=NAME14,STATUS='UNKNOWN')
      OPEN(2140,FILE=NAME15,STATUS='UNKNOWN')
      OPEN(2150,FILE=NAME16,STATUS='UNKNOWN')
      OPEN(2160,FILE=NAME17,STATUS='UNKNOWN')
      OPEN(2170,FILE=NAME18,STATUS='UNKNOWN')
      OPEN(2180,FILE=NAME19,STATUS='UNKNOWN')
      OPEN(2190,FILE=NAME20,STATUS='UNKNOWN')
      OPEN(2200,FILE=NAME21,STATUS='UNKNOWN')
      OPEN(2210,FILE=NAME22,STATUS='UNKNOWN')
      OPEN(2220,FILE=NAME23,STATUS='UNKNOWN')
      OPEN(2230,FILE=NAME24,STATUS='UNKNOWN')
      OPEN(2240,FILE=NAME25,STATUS='UNKNOWN')
      OPEN(2250,FILE=NAME26,STATUS='UNKNOWN')
      OPEN(2260,FILE=NAME27,STATUS='UNKNOWN')
      OPEN(2270,FILE=NAME28,STATUS='UNKNOWN')
      OPEN(2280,FILE=NAME29,STATUS='UNKNOWN')
      OPEN(2290,FILE=NAME30,STATUS='UNKNOWN')
      

      WRITE(2000,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2010,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2020,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2030,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2040,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2050,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2060,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2070,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2080,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2090,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2100,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2110,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2120,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2130,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2140,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2150,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2160,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2170,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2180,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2190,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2200,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2210,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2220,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2230,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2240,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2250,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2260,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2270,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2280,*) 'ZONE F=POINT,I=201,J=51'
      WRITE(2290,*) 'ZONE F=POINT,I=201,J=51'

      DO 800 J = 0, NY
        DO 800 I = 0, NX

      WRITE(2000,*) XP(I,J)*D,YP(I,J)*D,RP(I,J)
      WRITE(2010,*) XP(I,J)*D,YP(I,J)*D,AP(I,J)
      WRITE(2020,*) XP(I,J)*D,YP(I,J)*D,APR(I,J)
      WRITE(2030,*) XP(I,J)*D,YP(I,J)*D,APS(I,J)
      WRITE(2040,*) XP(I,J)*D,YP(I,J)*D,f_I(I,J)
      WRITE(2050,*) XP(I,J)*D,YP(I,J)*D,f_IA(I,J)
      WRITE(2060,*) XP(I,J)*D,YP(I,J)*D,f_II(I,J)
      WRITE(2070,*) XP(I,J)*D,YP(I,J)*D,f_IIA(I,J)
      WRITE(2080,*) XP(I,J)*D,YP(I,J)*D,f_V(I,J)
      WRITE(2090,*) XP(I,J)*D,YP(I,J)*D,f_VA(I,J)
      WRITE(2100,*) XP(I,J)*D,YP(I,J)*D,f_VIII(I,J)
      WRITE(2110,*) XP(I,J)*D,YP(I,J)*D,f_VIIIA(I,J)
      WRITE(2120,*) XP(I,J)*D,YP(I,J)*D,f_IX(I,J)
      WRITE(2130,*) XP(I,J)*D,YP(I,J)*D,f_IXA(I,J)
      WRITE(2140,*) XP(I,J)*D,YP(I,J)*D,f_X(I,J)
      WRITE(2150,*) XP(I,J)*D,YP(I,J)*D,f_XA(I,J)
      WRITE(2160,*) XP(I,J)*D,YP(I,J)*D,f_XI(I,J)
      WRITE(2170,*) XP(I,J)*D,YP(I,J)*D,f_XIA(I,J)
      WRITE(2180,*) XP(I,J)*D,YP(I,J)*D,ATIII(I,J)
      WRITE(2190,*) XP(I,J)*D,YP(I,J)*D,TFPI(I,J)
      WRITE(2200,*) XP(I,J)*D,YP(I,J)*D,PC(I,J)
      WRITE(2210,*) XP(I,J)*D,YP(I,J)*D,APC(I,J)
      WRITE(2220,*) XP(I,J)*D,YP(I,J)*D,rL1AT(I,J)
      WRITE(2230,*) XP(I,J)*D,YP(I,J)*D,TPA(I,J)
      WRITE(2240,*) XP(I,J)*D,YP(I,J)*D,PLS(I,J)
      WRITE(2250,*) XP(I,J)*D,YP(I,J)*D,PLA(I,J)
      WRITE(2260,*) XP(I,J)*D,YP(I,J)*D,rL2AP(I,J)
      WRITE(2270,*) XP(I,J)*D,YP(I,J)*D,W(I,J)
      WRITE(2280,*) XP(I,J)*D,YP(I,J)*D,Z(I,J)
  800 WRITE(2290,*) XP(I,J)*D,YP(I,J)*D,KPA(I,J)  

      DO 212 J = 1, NY-1
        DO 212 I = 1, NX-1

       IF(XP(I,J)*D.GE.0.0) THEN

      TRP = TRP + DSS(I,J)*RP(I,J)*2*Pi*YP(I,J)*D
      TAP = TAP + DSS(I,J)*AP(I,J)*2*Pi*YP(I,J)*D
      TAPR = TAPR + DSS(I,J)*APR(I,J)*2*Pi*YP(I,J)*D
      TAPS = TAPS + DSS(I,J)*APS(I,J)*2*Pi*YP(I,J)*D
      Tf_I = Tf_I + DSS(I,J)*f_I(I,J)*2*Pi*YP(I,J)*D
      Tf_IA = Tf_IA + DSS(I,J)*f_IA(I,J)*2*Pi*YP(I,J)*D
      Tf_II = Tf_II + DSS(I,J)*f_II(I,J)*2*Pi*YP(I,J)*D
      Tf_IIA = Tf_IIA + DSS(I,J)*f_IIA(I,J)*2*Pi*YP(I,J)*D
      Tf_V = Tf_V + DSS(I,J)*f_V(I,J)*2*Pi*YP(I,J)*D
      Tf_VA = Tf_VA + DSS(I,J)*f_VA(I,J)*2*Pi*YP(I,J)*D
      Tf_VIII = Tf_VIII + DSS(I,J)*f_VIII(I,J)*2*Pi*YP(I,J)*D
      Tf_VIIIA = Tf_VIIIA + DSS(I,J)*f_VIIIA(I,J)*2*Pi*YP(I,J)*D
      Tf_IX = Tf_IX + DSS(I,J)*f_IX(I,J)*2*Pi*YP(I,J)*D
      Tf_IXA = Tf_IXA + DSS(I,J)*f_IXA(I,J)*2*Pi*YP(I,J)*D
      Tf_X = Tf_X + DSS(I,J)*f_X(I,J)*2*Pi*YP(I,J)*D
      Tf_XA = Tf_XA + DSS(I,J)*f_XA(I,J)*2*Pi*YP(I,J)*D
      Tf_XI = Tf_XI + DSS(I,J)*f_XI(I,J)*2*Pi*YP(I,J)*D
      Tf_XIA = Tf_XIA + DSS(I,J)*f_XIA(I,J)*2*Pi*YP(I,J)*D
      TW = TW + DSS(I,J)*W(I,J)*2*Pi*YP(I,J)*D
      TZ = TZ + DSS(I,J)*Z(I,J)*2*Pi*YP(I,J)*D
      TPC = TPC + DSS(I,J)*PC(I,J)*2*Pi*YP(I,J)*D
      TAPC = TAPC + DSS(I,J)*APC(I,J)*2*Pi*YP(I,J)*D
      TATIII = TATIII + DSS(I,J)*ATIII(I,J)*2*Pi*YP(I,J)*D
      TTFPI = TTFPI + DSS(I,J)*TFPI(I,J)*2*Pi*YP(I,J)*D
      TrL1AT = TrL1AT + DSS(I,J)*rL1AT(I,J)*2*Pi*YP(I,J)*D
      TTPA = TTPA + DSS(I,J)*TPA(I,J)*2*Pi*YP(I,J)*D
      TPLS = TPLS + DSS(I,J)*PLS(I,J)*2*Pi*YP(I,J)*D
      TPLA = TPLA + DSS(I,J)*PLA(I,J)*2*Pi*YP(I,J)*D
      TrL2AP = TrL2AP + DSS(I,J)*rL2AP(I,J)*2*Pi*YP(I,J)*D

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
      OPEN(3010,FILE='TAP',STATUS='UNKNOWN')
      OPEN(3020,FILE='TAPR',STATUS='UNKNOWN')
      OPEN(3030,FILE='TAPS',STATUS='UNKNOWN')
      OPEN(3040,FILE='TI',STATUS='UNKNOWN')
      OPEN(3050,FILE='TIA',STATUS='UNKNOWN')
      OPEN(3060,FILE='TII',STATUS='UNKNOWN')
      OPEN(3070,FILE='TIIA',STATUS='UNKNOWN')
      OPEN(3080,FILE='TV',STATUS='UNKNOWN')
      OPEN(3090,FILE='TVA',STATUS='UNKNOWN')
      OPEN(3100,FILE='TVIII',STATUS='UNKNOWN')
      OPEN(3110,FILE='TVIIIA',STATUS='UNKNOWN')
      OPEN(3120,FILE='TIX',STATUS='UNKNOWN')
      OPEN(3130,FILE='TIXA',STATUS='UNKNOWN')
      OPEN(3140,FILE='TX',STATUS='UNKNOWN')
      OPEN(3150,FILE='TXA',STATUS='UNKNOWN')
      OPEN(3160,FILE='TXI',STATUS='UNKNOWN')
      OPEN(3170,FILE='TXIA',STATUS='UNKNOWN')
      OPEN(3180,FILE='TW',STATUS='UNKNOWN')
      OPEN(3190,FILE='TZ',STATUS='UNKNOWN')
      OPEN(3200,FILE='TPC',STATUS='UNKNOWN')
      OPEN(3210,FILE='TAPC',STATUS='UNKNOWN')
      OPEN(3220,FILE='TATIII',STATUS='UNKNOWN')
      OPEN(3230,FILE='TTFPI',STATUS='UNKNOWN')
      OPEN(3240,FILE='TrL1AT',STATUS='UNKNOWN')
      OPEN(3250,FILE='TTPA',STATUS='UNKNOWN')
      OPEN(3260,FILE='TPLS',STATUS='UNKNOWN')
      OPEN(3270,FILE='TPLA',STATUS='UNKNOWN')
      OPEN(3280,FILE='TrL2AP',STATUS='UNKNOWN')
      OPEN(3290,FILE='TDEPOSITION',STATUS='UNKNOWN')
      OPEN(3300,FILE='TDEPRP',STATUS='UNKNOWN')
      OPEN(3310,FILE='TDEPAP',STATUS='UNKNOWN')
      OPEN(5000,FILE='TDEPAPP',STATUS='UNKNOWN')
      OPEN(5100,FILE='MAT150',STATUS='UNKNOWN')
      OPEN(5200,FILE='APTNumRate',STATUS='UNKNOWN')
      OPEN(5300,FILE='APpoint',STATUS='UNKNOWN')
      OPEN(5400,FILE='Information Needed',STATUS='UNKNOWN')

      WRITE(3000,*) TRP
      WRITE(3010,*) TAP
      WRITE(3020,*) TAPR
      WRITE(3030,*) TAPS
      WRITE(3040,*) Tf_I
      WRITE(3050,*) Tf_IA
      WRITE(3060,*) Tf_II
      WRITE(3070,*) Tf_IIA
      WRITE(3080,*) Tf_V
      WRITE(3090,*) Tf_VA
      WRITE(3100,*) Tf_VIII
      WRITE(3110,*) Tf_VIIIA
      WRITE(3120,*) Tf_IX
      WRITE(3130,*) Tf_IXA
      WRITE(3140,*) Tf_X
      WRITE(3150,*) Tf_XA
      WRITE(3160,*) Tf_XI
      WRITE(3170,*) Tf_XIA
      WRITE(3180,*) TW
      WRITE(3190,*) TZ
      WRITE(3200,*) TPC
      WRITE(3210,*) TAPC
      WRITE(3220,*) TATIII
      WRITE(3230,*) TTFPI
      WRITE(3240,*) TrL1AT
      WRITE(3250,*) TTPA
      WRITE(3260,*) TPLS
      WRITE(3270,*) TPLA
      WRITE(3280,*) TrL2AP
      WRITE(3290,*) TDEP
      WRITE(3300,*) TDEPRP
      WRITE(3310,*) TDEPAP
      WRITE(5000,*) TDEPAPP
      WRITE(5100,*)  MAT(150)
      WRITE(5200,*) TAPP
      WRITE(5300,*) AP(138,49),AP(138,50)
      WRITE(5400,*) AP(157,49),AP(157,50),MAS(157)


        OPEN(50,FILE=NAME31,STATUS='UNKNOWN')
        WRITE(50,*) MTO

        OPEN(51,FILE=NAME32,STATUS='UNKNOWN')
        WRITE(51,*) M

        OPEN(52,FILE=NAME33,STATUS='UNKNOWN')
        WRITE(52,*) MAS

        OPEN(53,FILE=NAME34,STATUS='UNKNOWN')
        WRITE(53,*) MR

        OPEN(54,FILE=NAME35,STATUS='UNKNOWN')
        WRITE(54,*) MAT

        OPEN(4000,FILE=NAME36,STATUS='UNKNOWN')
        OPEN(4010,FILE=NAME37,STATUS='UNKNOWN')
        OPEN(4020,FILE=NAME38,STATUS='UNKNOWN')
        OPEN(4030,FILE=NAME39,STATUS='UNKNOWN')
        OPEN(4040,FILE=NAME40,STATUS='UNKNOWN')
        OPEN(4050,FILE=NAME41,STATUS='UNKNOWN')
        OPEN(4060,FILE=NAME42,STATUS='UNKNOWN')
        OPEN(4070,FILE=NAME43,STATUS='UNKNOWN')
        OPEN(4080,FILE=NAME44,STATUS='UNKNOWN')
        OPEN(4090,FILE=NAME45,STATUS='UNKNOWN')
        OPEN(4100,FILE=NAME46,STATUS='UNKNOWN')
        OPEN(4110,FILE=NAME47,STATUS='UNKNOWN')
        OPEN(4120,FILE=NAME48,STATUS='UNKNOWN')
        OPEN(4130,FILE=NAME49,STATUS='UNKNOWN')
        OPEN(4140,FILE=NAME50,STATUS='UNKNOWN')
        OPEN(4150,FILE=NAME51,STATUS='UNKNOWN')
        OPEN(4160,FILE=NAME52,STATUS='UNKNOWN')
        OPEN(4170,FILE=NAME53,STATUS='UNKNOWN')
        OPEN(4180,FILE=NAME54,STATUS='UNKNOWN')
        OPEN(4190,FILE=NAME55,STATUS='UNKNOWN')
        OPEN(4200,FILE=NAME56,STATUS='UNKNOWN')
        OPEN(4210,FILE=NAME57,STATUS='UNKNOWN')
        OPEN(4220,FILE=NAME58,STATUS='UNKNOWN')
        OPEN(4230,FILE=NAME59,STATUS='UNKNOWN')
        OPEN(4240,FILE=NAME60,STATUS='UNKNOWN')
        OPEN(4250,FILE=NAME61,STATUS='UNKNOWN')
        OPEN(4260,FILE=NAME62,STATUS='UNKNOWN')
        OPEN(4270,FILE=NAME63,STATUS='UNKNOWN')
        OPEN(4280,FILE=NAME64,STATUS='UNKNOWN')

        DO 900 I = 0, NX

      WRITE(4000,*) RP(I,NY)
      WRITE(4010,*) AP(I,NY)
      WRITE(4020,*) APR(I,NY)
      WRITE(4030,*) APS(I,NY)
      WRITE(4040,*) f_I(I,NY)
      WRITE(4050,*) f_IA(I,NY)
      WRITE(4060,*) f_II(I,NY)
      WRITE(4070,*) f_IIA(I,NY)
      WRITE(4080,*) f_V(I,NY)
      WRITE(4090,*) f_VA(I,NY)
      WRITE(4100,*) f_VIII(I,NY)
      WRITE(4110,*) f_VIIIA(I,NY)
      WRITE(4120,*) f_IX(I,NY)
      WRITE(4130,*) f_IXA(I,NY)
      WRITE(4140,*) f_X(I,NY)
      WRITE(4150,*) f_XA(I,NY)
      WRITE(4160,*) f_XI(I,NY)
      WRITE(4170,*) f_XIA(I,NY)
      WRITE(4180,*) ATIII(I,NY)
      WRITE(4190,*) TFPI(I,NY)
      WRITE(4200,*) PC(I,NY)
      WRITE(4210,*) APC(I,NY)
      WRITE(4220,*) rL1AT(I,NY)
      WRITE(4230,*) TPA(I,NY)
      WRITE(4240,*) PLS(I,NY)
      WRITE(4250,*) PLA(I,NY)
      WRITE(4260,*) rL2AP(I,NY)
      WRITE(4270,*) W(I,NY)
  900 WRITE(4280,*) Z(I,NY)

       ENDIF

      KTT=KTT+1

      TS = 0D0
      TRP = 0D0
      TAP = 0D0
      TAPP = 0D0
      TAPR = 0D0
      TAPS = 0D0
      Tf_I = 0D0
      Tf_IA = 0D0
      Tf_II = 0D0
      Tf_IIA = 0D0
      Tf_V = 0D0
      Tf_VA = 0D0
      Tf_VIII = 0D0
      Tf_VIIIA = 0D0
      Tf_IX = 0D0
      Tf_IXA = 0D0
      Tf_X = 0D0
      Tf_XA = 0D0
      Tf_XI = 0D0
      Tf_XIA = 0D0
      TW = 0D0
      TZ = 0D0
      TPC = 0D0
      TAPC = 0D0
      TATIII = 0D0
      TTFPI = 0D0
      TrL1AT = 0D0
      TTPA = 0D0
      TPLS = 0D0
      TPLA = 0D0
      TrL2AP = 0D0
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
     &              f_I(-1:NX+1,-1:NY+1),f_IA(-1:NX+1,-1:NY+1),
     &              f_II(-1:NX+1,-1:NY+1),f_IIA(-1:NX+1,-1:NY+1),
     &              f_V(-1:NX+1,-1:NY+1),f_VA(-1:NX+1,-1:NY+1),
     &              f_VIII(-1:NX+1,-1:NY+1),f_VIIIA(-1:NX+1,-1:NY+1),
     &              f_IX(-1:NX+1,-1:NY+1),f_IXA(-1:NX+1,-1:NY+1),
     &              f_X(-1:NX+1,-1:NY+1),f_XA(-1:NX+1,-1:NY+1),
     &              f_XI(-1:NX+1,-1:NY+1),f_XIA(-1:NX+1,-1:NY+1),
     &              ATIII(-1:NX+1,-1:NY+1),Z(-1:NX+1,-1:NY+1),
     &              W(-1:NX+1,-1:NY+1),TFPI(-1:NX+1,-1:NY+1),
     &              APC(-1:NX+1,-1:NY+1),PC(-1:NX+1,-1:NY+1),
     &              rL1AT(-1:NX+1,-1:NY+1),TPA(-1:NX+1,-1:NY+1),
     &              PLA(-1:NX+1,-1:NY+1),PLS(-1:NX+1,-1:NY+1),
     &              rL2AP(-1:NX+1,-1:NY+1)
      COMMON /BLK2/ U(0:NX,0:NY),V(0:NX,0:NY)
      COMMON /BLK3/ RPT(0:NX,0:NY),APT(0:NX,0:NY),
     &              APRT(0:NX,0:NY),APST(0:NX,0:NY), 
     &              f_IT(0:NX,0:NY),f_IAT(0:NX,0:NY),
     &              f_IIT(0:NX,0:NY),f_IIAT(0:NX,0:NY),
     &              f_VT(0:NX,0:NY),f_VAT(0:NX,0:NY),
     &              f_VIIIT(0:NX,0:NY),f_VIIIAT(0:NX,0:NY),
     &              f_IXT(0:NX,0:NY),f_IXAT(0:NX,0:NY),
     &              f_XT(0:NX,0:NY),f_XAT(0:NX,0:NY),
     &              f_XIT(0:NX,0:NY),f_XIAT(0:NX,0:NY),
     &              ATIIIT(0:NX,0:NY),ZT(0:NX,0:NY),
     &              WT(0:NX,0:NY),TFPIT(0:NX,0:NY),
     &              APCT(0:NX,0:NY),PCT(0:NX,0:NY),
     &              rL1ATT(0:NX,0:NY),TPAT(0:NX,0:NY),
     &              PLAT(0:NX,0:NY),PLST(0:NX,0:NY),
     &              rL2APT(0:NX,0:NY)
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
     &              f_I(-1:NX+1,-1:NY+1),f_IA(-1:NX+1,-1:NY+1),
     &              f_II(-1:NX+1,-1:NY+1),f_IIA(-1:NX+1,-1:NY+1),
     &              f_V(-1:NX+1,-1:NY+1),f_VA(-1:NX+1,-1:NY+1),
     &              f_VIII(-1:NX+1,-1:NY+1),f_VIIIA(-1:NX+1,-1:NY+1),
     &              f_IX(-1:NX+1,-1:NY+1),f_IXA(-1:NX+1,-1:NY+1),
     &              f_X(-1:NX+1,-1:NY+1),f_XA(-1:NX+1,-1:NY+1),
     &              f_XI(-1:NX+1,-1:NY+1),f_XIA(-1:NX+1,-1:NY+1),
     &              ATIII(-1:NX+1,-1:NY+1),Z(-1:NX+1,-1:NY+1),
     &              W(-1:NX+1,-1:NY+1),TFPI(-1:NX+1,-1:NY+1),
     &              APC(-1:NX+1,-1:NY+1),PC(-1:NX+1,-1:NY+1),
     &              rL1AT(-1:NX+1,-1:NY+1),TPA(-1:NX+1,-1:NY+1),
     &              PLA(-1:NX+1,-1:NY+1),PLS(-1:NX+1,-1:NY+1),
     &              rL2AP(-1:NX+1,-1:NY+1)
      COMMON /BLK2/ U(0:NX,0:NY),V(0:NX,0:NY)
      COMMON /BLK3/ RPT(0:NX,0:NY),APT(0:NX,0:NY),
     &              APRT(0:NX,0:NY),APST(0:NX,0:NY), 
     &              f_IT(0:NX,0:NY),f_IAT(0:NX,0:NY),
     &              f_IIT(0:NX,0:NY),f_IIAT(0:NX,0:NY),
     &              f_VT(0:NX,0:NY),f_VAT(0:NX,0:NY),
     &              f_VIIIT(0:NX,0:NY),f_VIIIAT(0:NX,0:NY),
     &              f_IXT(0:NX,0:NY),f_IXAT(0:NX,0:NY),
     &              f_XT(0:NX,0:NY),f_XAT(0:NX,0:NY),
     &              f_XIT(0:NX,0:NY),f_XIAT(0:NX,0:NY),
     &              ATIIIT(0:NX,0:NY),ZT(0:NX,0:NY),
     &              WT(0:NX,0:NY),TFPIT(0:NX,0:NY),
     &              APCT(0:NX,0:NY),PCT(0:NX,0:NY),
     &              rL1ATT(0:NX,0:NY),TPAT(0:NX,0:NY),
     &              PLAT(0:NX,0:NY),PLST(0:NX,0:NY),
     &              rL2APT(0:NX,0:NY)
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


