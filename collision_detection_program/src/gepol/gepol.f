C
C       =============================================================
C
C       GEPOL/87 (GEometría POlihedro,1987)
C
C       =============================================================
C
C*** Written by
C
C       J.L. Pascual-Ahuir and E. Silla
C       Departamento de Química Física.
C       Facultad de Química.
C       Universidad de Valencia.
C       C/Dr.Moliner 50.
C       Burjassot (Valencia) 46100
C       SPAIN
C       
C       Phone number: 6-3630011  ext. 341
C
C       J. Tomasi
C       Dipartamento di Chimica.
C       Università di Pisa.
C       Via Risorgimento, 35-I 
C       56100 Pisa
C       ITALY
C
C       R. Bonaccorsi.
C       Istituto di Chimica Quantistica e Energetica Moleculare (CNR)
C       Via Risorgimento, 35-I 
C       56100 Pisa
C       ITALY
C       
C*** Grants
C
C       It has been granted by Comision Asesora Cientifica Y Tecnica
C       (CACYT) of Spain, proyect 714/84.
C
C
C*** References
C
C    A- J.L.Pascual-Ahuir, E.Silla, J.Tomasi y R.Bonacorsi.
C       Electrostatic Interaction of Solute with a Continuum.
C       Improved Description of the Cavity and of the Surface Cavity
C       Bound Charge Distribution.
C       Journal of Computational Chemistry, Vol.8, No.6,778-787(1987)
C
C    B- GEPOL a program to calculate the envelop molecular surface.
C       I. The van der Waals Surface.
C       To be published. 
C
C   C- GEPOL a program to calculate the envelop molecular surface.
C       II. The Molecular Surface.
C       To be published. 
C
C
C*** Aim
C
C       This program calculates the envelop surface for a molecule, as a
C       points distribution, and computes its  correspondent area and volume.
C       As data the coordinates of atoms, the radii assigned to them and
C       the solvent radius (Probe radius) are needed.
C       
C       It calculates the three kind of envelop surfaces:
C
C         - THE VAN DER WAALS MOLECULAR SURFACE. It is the surface 
C       resulting from a set of intersecting spheres with radii of van
C       der Waals, centred on selected nuclei of the molecule.
C
C         - THE ACCESSIBLE MOLECULAR SURFACE. It was defined by
C       B.Lee and F.M.Richards (J.MOL.BIOL.55(1971)379-400. It is the
C       surface defined by the center of the solvent sphere (probe sphere)
C       when it rolls around the van der waals surface. 
C
C         - THE MOLECULAR SURFACE. It was defined by F.M.Richards 
C       (Ann. Rev. Biophys. Bioeng.,6 (1977)151-176).
C       The present program calculates the Molecular Surface creating 
C       a set of new spheres locating among the original spheres 
C       defined at the input geomtry.
C       At this proces we call 'smoothing' of the surface.
C 
C       
C       
C*** Computacional especifications
C
C       This program is written in FORTRAN 77 and runed over VAX VMS.
C              
C*** Input
C       
C       File    FOR005
C
C       The input is composed of the next records:
C     -- 1st record.-- Format (80A1)
C            
C            Title
C
C  
C     -- 2nd record.-- Format(4(3A1,1X))
C         
C         It is composed of four character variables. The options
C         must be written with capital letters.
C      
C         Variable         Options
C     
C         SUAVE             CON   Smoothing.
C                           SIN   No smoothing.
C
C         VECTOR            VES   Print the points and the corresponding vectors
C                                 on the surface.
C                           VEN   No print.
C           
C         SUMA              SUM   Adding the probe radius plus the van der
C                                 Waals radi of atoms.
C                           NSU   No adding.
C
C         OPESF             SFY   Print information about spheres.
C                           SFN   No print.
C
C         If the van der Waals surface is requerided, SUAVE will be
C         SIN and SUMA will be NSU.
C         If you want the accesible surface, SUAVE=SIN and SUMA=SUM
C         If you want the molecular surface, SUAVE=CON and SUMA=NSU
C
C
C     -- 3rd record -- Format(I2,2F10.5)
C                
C         Variable
C           
C           NDIV     It specifies the division level for the triangles on the
C                    surface. 
C                    The accuracy improve when NDIV raise.
C                    We think that NDIV equal to 3 is a good value (reference B
C                    is concerned with study about this point).
C
C           RD       It is the probe radious.
C                    It is only necessary if SUAVE=CON or SUMA=SUM.
C           
C           FRADIO   This parameter controls the new spheres creation
C                    (smoothing).One may introduced values in the range from
C                    0.0 to 1.0. It is only need if SUAVE=CON.
C                    The accuracy improves when the FRADIO value decreases.
C                    We recomend values are about 0.55 (see reference C).
C                           
C
C
C     -- 4th record -- Format(I8)                 
C          
C         Variable   
C         
C           NATOM        Number of atoms.
C
C
C     -- 5th record -- Format(4F10.5,I2)
C        
C           One card by atom.
C
C         Variable
C
C           XE        Coordinate X of the atom. At the same units that RD.
C           YE        Coordinate Y of the atom. "   "   "    "     "   "
C           ZE        Coordinate Z of the atom. "   "   "    "     "   "
C           RE        Radious of the sphere centred on atom.
C                     At the same units that RD.
C           USE    0  This atom is not used for the surface.
C                  1  This atom is used for the surface.
C
C
C*** Output files
C
C           FOR006 Abstract.
C           FOR007 Spheres information.
C           FOR008 Vectors information.
C
C
C*** Example command file for VAX 
C           $ ASSIGN EXAMPLE.INP FOR005
C           $ ASSIGN EXAMPLE.OUT FOR006
C           $ ASSIGN EXAMPLE.SFE FOR007
C           $ ASSIGN EXAMPLE.VEC FOR008
C           $ RUN GEPOL
C
C*** For more information see the user manual.
C
C
C
C    ******************************************************************
C    *  Version for Personal Iris (Silicon Graphics)                   * 
C    *  Some modifications has been made to prepare gepol outputs      * 
C    *  for the TOM Graphics Program:                                  * 
C    *             1.A new output file called EXAMPLE.DIB has been     *
C    *               added. It has been assigned to FOR015. It         *
C    *               conteins information about triangles vertex       * 
C    *               and you can obtein it giving the option DIB       *
C    *               to the variable DIBUJO (in the same record        * 
C    *               and format that variables SUAVE.....)             * 
C    *             2.Information about the kind of spheres is          *
C    *               given in .VEC and .DIB output files with          *
C    *               the variable spheretype (STY) in the output.      *
C    *               STY=1 spheres belonging to van der Waals surface  *
C    *               STY=2 spheres of first generation                 *
C    *               STY=3 spheres of second generation                *
C    *               STY=4 .......                                     *
C    *                                                                 *
C    * By  Fernando Villar & Inyaki Tunyon                             *
C    *     Department of Physical Chemistry                            *
C    *     University of Valencia (Spain)                              *
C    *     C/ Doctor Moliner 50                                        *
C    *        46100-Burjassot (Valencia)                               *
C    *                         SPAIN                                   *
C    *                                                                 *  
C    *******************************************************************
C
C
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*4       (I-N) 
      CHARACTER*3 SUAVE,vector,suma,OPESF,DIBUJO
      CHARACTER*80 TITULO
      LOGICAL NOUSA  
      REAL*8 CV
C        PARAMETER nst=100000
C      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(NST)
C      COMMON/CSFE/XE(NST),YE(NST),ZE(NST),RE(NST),NPEC(NST,2)
      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(100000)
      COMMON/CSFE/XE(100000),YE(100000),ZE(100000),RE(100000)
     &,NPEC(100000,2)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/ESC/TITULO,SUAVE,vector,SUMA,OPESF,FRADIO
     &,DIBUJO
      COMMON/TYPE/ISPHERETYPE(100000)
      DATA DVEC/-1.0/
      DATA OMEGA/50.0/
      NST=100000
C*****Record 1******
      READ(5,'(A)')TITULO
C*****Record 2******
      READ(5,'(6(A,1X))')SUAVE,vector,suma,OPESF,
     &DIBUJO
C*****Record 3*******
      READ(5,'(I2,2F10.5)')NDIV,RD,FRADIO
      CALL GWRITE(6)
      IF(OPESF.EQ.'SFY') CALL GWRITE(7)
      IF(VECTOR.EQ.'VES')CALL GWRITE(8)
      IF(DIBUJO.EQ.'DIB')CALL GWRITE(15)
C*****Records 4 Y 5********
      CALL STAND
C*****Begining ******
      CALL TES
      CALL DIVIDE(NDIV)
      NESFI=NESF
      IF(SUAVE.EQ.'CON') CALL CREA
      NESFNU=NESF-NESFI
      WRITE(6,'(A,I6)')' Number of created spheres=',NESFNU
      CALL GEOCAV
      STOP
      END
C
      SUBROUTINE STAND
C**************************************************************
      IMPLICIT REAL*4          (A-H,O-Z)
      IMPLICIT INTEGER*4       (I-N)
      INTEGER*4 USE
      LOGICAL NOUSA
C      PARAMETER nst=100000
      CHARACTER*3 SUAVE,vector,suma,OPESF,DIBUJO
      CHARACTER*80 TITULO
      COMMON/ESC/TITULO,SUAVE,vector,SUMA,OPESF,FRADIO
     &,DIBUJO
C      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(NST)
C      COMMON/CSFE/XE(NST),YE(NST),ZE(NST),RE(NST),NPEC(NST,2)      
C      COMMON/TYPE/ISPHERETYPE(NST)
      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(100000)
      COMMON/CSFE/XE(100000),YE(100000),ZE(100000),RE(100000)
     &,NPEC(100000,2)
      COMMON/TYPE/ISPHERETYPE(100000)

      NST=100000
C*****Record 4******
      READ (5,'(I8)') NATOM
      WRITE(6,'(A,I8)')' NATOM=', NATOM
      NESF=NATOM
C*****Record 5******
      INDICE=0
      RTOTAL=0.0
      DO 10 I = 1,NATOM
      READ (5,'(4F10.5,I2,I2)') XE(I),YE(I),ZE(I),RE(I),
     &USE,ISPHERETYPE(I)
      IF(SUMA.EQ.'SUM')RE(I)=RE(I)+RD
      IF(USE.EQ.0)RE(I)=0.0
      RTOTAL=RTOTAL+RE(I)
      IF(I.LT.2) WRITE(6,'(4F10.5,I2)') XE(I),YE(I),ZE(I),RE(I),USE
      NOUSA(I)=.FALSE.
      IF(USE.EQ.1)GO TO 10
      NOUSA(I)=.TRUE.
      INDICE=INDICE+1
   10 CONTINUE
      NATESF=NESF-INDICE
      RMEDIO=RTOTAL/NATESF
      RET=RMEDIO*FRADIO
      FRO=FRADIO
      WRITE(6,'(2(A,I6/))')
     & ' Number of atoms used to create the surface',NATESF,
     & '   "    "    " not  "   "   "     "     "   ',INDICE
      DO I=1,NATOM
          ISPHERETYPE(I)=1
      END DO
      RETURN
      END
C
      SUBROUTINE CREA
C*********************************************************************
C     This subroutine creates the new spheres
C*********************************************************************
      IMPLICIT INTEGER*4 (I,N)
      IMPLICIT REAL*4 (A-H,O-Z)
      INTEGER*4 SITUA
      CHARACTER*3 SUAVE,vector,suma,OPESF,DIBUJO
      CHARACTER*80 TITULO
      LOGICAL NOUSA
C      PARAMETER nst=100000
C      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(NST)
      COMMON/ESC/TITULO,SUAVE,vector,SUMA,OPESF,FRADIO
     &,DIBUJO
C      COMMON/CSFE/XE(NST),YE(NST),ZE(NST),RE(NST),NPEC(NST,2)
C      COMMON/TYPE/ISPHERETYPE(NST)
      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(100000)
      COMMON/CSFE/XE(100000),YE(100000),ZE(100000),RE(100000)
     &,NPEC(100000,2)
      COMMON/TYPE/ISPHERETYPE(100000)


C      DIMENSION SITUA(NST)
      DIMENSION SITUA(100000)
      DATA PI/3.1415926535897932D0/
      NST=100000
      FIRST=PI/180.0D0  
      DO 2 I=1,NESF
      NPEC(I,1)=0
      NPEC(I,2)=0
    2 CONTINUE
      OMEGA=OMEGA*FIRST
      COSOM=COS(OMEGA)
      SENOM2=(SIN(OMEGA))**2
      RTDD=RET+RD
      RTDD2=RTDD*RTDD
      NET=NESF
      NN=2   
      NE=NESF
      NEV=NESF
      ITYPE=1
      GO TO 600
  601 continue
      NN=NE+1
      ITYPE=ITYPE+1
      NE=NET
C Loop to select the first sphere of the pair.
  600 continue
      DO 602 I=NN,NE
      NES=I-1
      ISPHERETYPE(I)=ITYPE
      IF(NOUSA(I)) GO TO 602
C
C It fixes how are the rest of spheres respect to I.
C Situa=0 the solvent can go between them.
C Situa=1 the solvent can not go between them.
C Situa=2 if they are linked.
      DO 720 J=1,NEV
      SITUA(J)=0
      IF(NOUSA(J)) GO TO 720
      IF(J.EQ.I)GO TO 720
      CALL DISTAN(XE(I),YE(I),ZE(I),XE(J),YE(J),ZE(J),D2)
      RJD=RE(I)+RE(J)
      RJD2=RJD*RJD
      IF(D2.LE.RJD2) GO TO 722
      TEST1=RD+RD+RJD
      TEST1=TEST1*TEST1
      IF(D2.LT.TEST1)SITUA(J)=1
      GO TO 720
  722 CONTINUE
      SITUA(J)=2
  720 CONTINUE
C
C It selects the second sphere of the pair.
      DO 603 J=1,NES
      IF(SITUA(J).EQ.0)GO TO 603
      CALL DISTAN(XE(I),YE(I),ZE(I),XE(J),YE(J),ZE(J),RIJ2)
      RIJ=SQRT(RIJ2)
  701 continue
      REG=AMAX1(RE(I),RE(J))
      REP=AMIN1(RE(I),RE(J))
      REG2=REG*REG
      REP2=REP*REP
      IF(SITUA(J).NE.2) GO TO 705
      TEST2=REP*COSOM+SQRT(REG2-REP2*SENOM2)                                   
      IF(RIJ.LE.TEST2) GO TO 603
  705 continue
      REGD2=(REG+RD)*(REG+RD)
      TEST3=(REGD2+REG2-RTDD2)/REG
      IF(RIJ.GE.TEST3) GO TO 603
C
C Visibility text
      DO 604 K=1,NEV
      IF(K.EQ.J) GO TO 604
      IF(SITUA(K).EQ.0) GO TO 604
      CALL DISTAN(XE(J),YE(J),ZE(J),XE(K),YE(K),ZE(K),RJK2)
      IF(RJK2.GE.RIJ2) GO TO 604
      CALL DISTAN(XE(I),YE(I),ZE(I),XE(K),YE(K),ZE(K),RIK2)
      IF(RIK2.GE.RIJ2) GO TO 604
       RJK=SQRT(RJK2)
       RIK=SQRT(RIK2)
       SP=(RIJ+RJK+RIK)/2.0D0
       HH=4*(SP*(SP-RIJ)*(SP-RIK)*(SP-RJK))/RIJ2
       REO=RE(K)*FRO
      IF(K.GE.NE)REO=0.0002D0
      REO2=REO*REO
      IF(HH.LT.REO2) GO TO 603
  604 CONTINUE
      REPD2=(REP+RD)**2
      IF(REP.LE.RET) GO TO 700
      IF(SITUA(J).EQ.2) GO TO 605
      TEST8=SQRT(REPD2-RTDD2)+SQRT(REGD2-RTDD2)
      IF(RIJ.LE.TEST8)GO TO 605
C
C Sphere of kind C
  700 continue
      REND2=REGD2+REG2-(REG/RIJ)*(REGD2+RIJ2-REPD2)
      IF(REND2.LE.RTDD2) GO TO 603
      REN=SQRT(REND2)-RD
      FC=REG/(RIJ-REG)
      TEST7=REG-RE(I)
      KG=I
      KP=J
      IF(TEST7.LE.0.000000001D0) GO TO 606
      KG=J
      KP=I
  606 continue
      FC1=FC+1.0
      XEN=(XE(KG)+FC*XE(KP))/FC1
      YEN=(YE(KG)+FC*YE(KP))/FC1
      ZEN=(ZE(KG)+FC*ZE(KP))/FC1
      GO TO 607
C
C Spheres of kind A or B
  605 continue
      R2GN=RIJ-REP+REG
      RGN=R2GN/2.0D0
      REN=SQRT(REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ))-RD
      IF(REN.LT.RET) GO TO 603
      FC=R2GN/(RIJ+REP-REG)
      FC1=FC+1.0D0
      TEST7=REG-RE(I)
      KG=I
      KP=J
      IF(TEST7.LE.0.000000001D0) GO TO 610
      KG=J
      KP=I
  610 continue
      XEN=(XE(KG)+FC*XE(KP))/FC1
      YEN=(YE(KG)+FC*YE(KP))/FC1
      ZEN=(ZE(KG)+FC*ZE(KP))/FC1
  607 continue
      NET=NET+1
      IF(NET.GT.NST)GO TO 21
      XE(NET)=XEN
      YE(NET)=YEN
      ZE(NET)=ZEN
      RE(NET)=REN
      SITUA(NET)=1
      NPEC(NET,1)=KG
      NPEC(NET,2)=KP
      NOUSA(NET)=.FALSE.
  603 CONTINUE
      NEV=NET
  602 CONTINUE
      IF(NET.NE.NE) GO TO 601
    8 continue
      NESF=NET
      RETURN
   21 continue
      WRITE(6,'(A)')'%-ERROR-% DIMENSION. TOO MUCH SPHERES CREATES'
      STOP
      END
C
      SUBROUTINE DISTAN(X,Y,Z,XX,YY,ZZ,D)
      IMPLICIT NONE
      REAL*4 X,Y,Z,XX,YY,ZZ,D,RX,RY,RZ
      RX=X-XX
      RY=Y-YY
      RZ=Z-ZZ
      D=RX*RX+RY*RY+RZ*RZ
      RETURN
      END
C
      SUBROUTINE TES
C*********************************************************************
C     It computes the triangles vertix coordinates for a sphere of radious
C     one, proyecting the pentakisdodecahedro on it.
C*********************************************************************
      IMPLICIT INTEGER*4 (I,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 XC1,YC1,ZC1
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      DIMENSION THEV(6),FIV(6)
      DATA THEV/0.6523581397843682D0,1.1071487177940905D0,
     $          1.3820857960113345D0,1.7595068575784587D0,
     $          2.0344439357957027D0,2.4892345138054251D0/
      DATA FIV/0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              /
      DATA FIR/1.2566370614359173 D0/
      CV(1,1)=0.D0
      CV(1,2)=0.D0
      CV(1,3)=1.D0
      CV(32,1)=0.D0
      CV(32,2)=0.D0
      CV(32,3)=-1.D0
      II=1
      DO 520 I=1,6
      TH=THEV(I)
      FI=FIV(I)
      CTH=DCOS(TH)
      STH=DSIN(TH)
      DO 521 J=1,5
      FI=FI+FIR
      IF(J.EQ.1) FI=FIV(I)
      II=II+1
      CV(II,1)=STH*DCOS(FI)
      CV(II,2)=STH*DSIN(FI)
      CV(II,3)=CTH
  521 CONTINUE
  520 CONTINUE
      RETURN
      END
C
      SUBROUTINE DIVIDE(NDIV)
C*********************************************
C     It divides the spherical triangles
C*********************************************
      IMPLICIT INTEGER*4 (I,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 XC1,YC1,ZC1
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/DIB/JVT1(3,60)
      DIMENSION JVT2(3,4),CVN2(6,3),
     &          CVN3(6,3),CVN4(6,3),CVN5(6,3),CC(3)
      DATA ZERO/0.0D0/,FOUR/4.0D0/
      DATA PI/3.1415926535897932D0/
      DATA JVT1 /  1,   6,   2,  1,   2,   3,  1,   3,   4,
     &             1,   4,   5,  1,   5,   6,  7,   2,   6,
     &             8,   3,   2,  9,   4,   3, 10,   5,   4,
     &            11,   6,   5,  8,   2,  12,  9,   3,  13,
     &            10,   4,  14, 11,   5,  15,  7,   6,  16,
     &             7,  12,   2,  8,  13,   3,  9,  14,   4,
     &            10,  15,   5, 11,  16,   6,  8,  12,  18,
     &             9,  13,  19, 10,  14,  20, 11,  15,  21,
     &             7,  16,  17,  7,  17,  12,  8,  18,  13,
     &             9,  19,  14, 10,  20,  15, 11,  21,  16,
     &            22,  12,  17, 23,  13,  18, 24,  14,  19,
     &            25,  15,  20, 26,  16,  21, 22,  18,  12,
     &            23,  19,  13, 24,  20,  14, 25,  21,  15,
     &            26,  17,  16, 22,  17,  27, 23,  18,  28,
     &            24,  19,  29, 25,  20,  30, 26,  21,  31,
     &            22,  28,  18, 23,  29,  19, 24,  30,  20,
     &            25,  31,  21, 26,  27,  17, 22,  27,  28,
     &            23,  28,  29, 24,  29,  30, 25,  30,  31,
     &            26,  31,  27, 32,  28,  27, 32,  29,  28,
     &            32,  30,  29, 32,  31,  30, 32,  27,  31 /
      DATA JVT2 /  1,   5,   4,
     &             5,   2,   6,
     &             4,   6,   3,
     &             6,   4,   5 /
      IJ=0
C*****Level 1****************************
      DO 10 J=1,60
      NV1=JVT1(1,J)
      NV2=JVT1(2,J)
      NV3=JVT1(3,J)
      XV1=CV(NV1,1)
      YV1=CV(NV1,2)
      ZV1=CV(NV1,3)
      XV2=CV(NV2,1)
      YV2=CV(NV2,2)
      ZV2=CV(NV2,3)
      XV3=CV(NV3,1)
      YV3=CV(NV3,2)
      ZV3=CV(NV3,3)
      IF(NDIV.GT.1) GO TO 20
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 10
C*****Level 2**********************
   20 continue
      CVN2(1,1)=XV1
      CVN2(1,2)=YV1
      CVN2(1,3)=ZV1
      CVN2(2,1)=XV2
      CVN2(2,2)=YV2
      CVN2(2,3)=ZV2
      CVN2(3,1)=XV3
      CVN2(3,2)=YV3
      CVN2(3,3)=ZV3
      CALL CALVER(CVN2)
      DO 21 J2=1,4
      NV21=JVT2(1,J2)
      NV22=JVT2(2,J2)
      NV23=JVT2(3,J2)
      XV1=CVN2(NV21,1)
      YV1=CVN2(NV21,2)
      ZV1=CVN2(NV21,3)
      XV2=CVN2(NV22,1)
      YV2=CVN2(NV22,2)
      ZV2=CVN2(NV22,3)
      XV3=CVN2(NV23,1)
      YV3=CVN2(NV23,2)
      ZV3=CVN2(NV23,3)
      IF(NDIV.GT.2) GO TO 30
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 21
C*****Level 3**********************************
   30 continue
      CVN3(1,1)=XV1
      CVN3(1,2)=YV1
      CVN3(1,3)=ZV1
      CVN3(2,1)=XV2
      CVN3(2,2)=YV2
      CVN3(2,3)=ZV2
      CVN3(3,1)=XV3
      CVN3(3,2)=YV3
      CVN3(3,3)=ZV3
      CALL CALVER(CVN3)
      DO 31 J3=1,4
      NV31=JVT2(1,J3)
      NV32=JVT2(2,J3)
      NV33=JVT2(3,J3)
      XV1=CVN3(NV31,1)
      YV1=CVN3(NV31,2)
      ZV1=CVN3(NV31,3)
      XV2=CVN3(NV32,1)
      YV2=CVN3(NV32,2)
      ZV2=CVN3(NV32,3)
      XV3=CVN3(NV33,1)
      YV3=CVN3(NV33,2)
      ZV3=CVN3(NV33,3)
      IF(NDIV.GT.3) GO TO 40
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 31
C*****Level 4******************************
   40 continue
      CVN4(1,1)=XV1
      CVN4(1,2)=YV1
      CVN4(1,3)=ZV1
      CVN4(2,1)=XV2
      CVN4(2,2)=YV2
      CVN4(2,3)=ZV2
      CVN4(3,1)=XV3
      CVN4(3,2)=YV3
      CVN4(3,3)=ZV3
      CALL CALVER(CVN4)
      DO 41 J4=1,4
      NV41=JVT2(1,J4)
      NV42=JVT2(2,J4)
      NV43=JVT2(3,J4)
      XV1=CVN4(NV41,1)
      YV1=CVN4(NV41,2)
      ZV1=CVN4(NV41,3)
      XV2=CVN4(NV42,1)
      YV2=CVN4(NV42,2)
      ZV2=CVN4(NV42,3)
      XV3=CVN4(NV43,1)
      YV3=CVN4(NV43,2)
      ZV3=CVN4(NV43,3)
      IF(NDIV.GT.4) GO TO 50
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 41
C*****Level 5*************************************
   50 continue
      CVN5(1,1)=XV1
      CVN5(1,2)=YV1
      CVN5(1,3)=ZV1
      CVN5(2,1)=XV2
      CVN5(2,2)=YV2
      CVN5(2,3)=ZV2
      CVN5(3,1)=XV3
      CVN5(3,2)=YV3
      CVN5(3,3)=ZV3
      CALL CALVER(CVN5)
      DO 51 J5=1,4
      NV51=JVT2(1,J5)
      NV52=JVT2(2,J5)
      NV53=JVT2(3,J5)
      XV1=CVN5(NV51,1)
      YV1=CVN5(NV51,2)
      ZV1=CVN5(NV51,3)
      XV2=CVN5(NV52,1)
      YV2=CVN5(NV52,2)
      ZV2=CVN5(NV52,3)
      XV3=CVN5(NV53,1)
      YV3=CVN5(NV53,2)
      ZV3=CVN5(NV53,3)
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
   51 CONTINUE
   41 CONTINUE
   31 CONTINUE
   21 CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALVER(CVN)
C---------------------------------------------------------------------
C  It divides one triangle in four.
C---------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 CVN,XXX,YYY,ZZZ,RRR,FC
      INTEGER*4 N,N1,N2
      DIMENSION CVN(6,3)
      DO 7 N=1,3
      N2=N+3
      N1=N-1
      IF(N.EQ.1)N1=3
      XXX=(CVN(N,1)+CVN(N1,1))/2.0D0
      YYY=(CVN(N,2)+CVN(N1,2))/2.0D0
      ZZZ=(CVN(N,3)+CVN(N1,3))/2.0D0
      RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      FC=1.0D0/RRR
      CVN(N2,1)=XXX*FC
      CVN(N2,2)=YYY*FC
      CVN(N2,3)=ZZZ*FC
    7 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
C---------------------------------------------------------------------
C   It computes the center of a spherical triangle.
C---------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC,XXX,YYY,ZZZ,RRR,FC
      DIMENSION CC(3)
      XXX=(XV1+XV2+XV3)/3.0D0
      YYY=(YV1+YV2+YV3)/3.0D0
      ZZZ=(ZV1+ZV2+ZV3)/3.0D0
      RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      FC=1.0D0/RRR
      CC(1)=XXX*FC
      CC(2)=YYY*FC
      CC(3)=ZZZ*FC
      RETURN
      END

C
      SUBROUTINE GEOCAV
C*********************************************************************
      IMPLICIT INTEGER*4 (I,N)
      IMPLICIT REAL*4 (A-H,O-Z)
      CHARACTER*3 SUAVE,vector,suma,OPESF,DIBUJO
      CHARACTER*6 MARCA
      CHARACTER*80 TITULO
      LOGICAL NOUSA
      REAL*8 CV
      INTEGER*4 IJ
C      PARAMETER nst=100000
C      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(NST)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
C      COMMON/CSFE/XE(NST),YE(NST),ZE(NST),RE(NST),NPEC(NST,2)
      COMMON/ESC/TITULO,SUAVE,vector,SUMA,OPESF,FRADIO
     &,DIBUJO
      COMMON/DIB/JVT1(3,60)
C      COMMON/TYPE/ISPHERETYPE(NST)
      DIMENSION VERTIX(3,3)
C      DIMENSION SRE(NST),SSFE(NST),NPE(NST),IJE(NST)
      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(100000)
      COMMON/CSFE/XE(100000),YE(100000),ZE(100000),RE(100000)
     &,NPEC(100000,2)
      COMMON/TYPE/ISPHERETYPE(100000)
      DIMENSION SRE(100000),SSFE(100000),NPE(100000),IJE(100000)
      
      DATA ZERO/0.0D0/,FOUR/4.0D0/
      DATA PI/3.1415926535897932D0/
      NST=100000
      STOT=ZERO
      VOL=ZERO
      IJ=0
      NES0=0
      NESI=0
C begin
      NTRIAN=4**(NDIV-1)
      FNDIV=FLOAT(60*NTRIAN)
      IF(VECTOR.EQ.'VES')WRITE(8,'(A,1X,A,5X,A,2(10X,A),8X,A,3(9X,A))')
     &'   SN','STY','X','Y','Z','AREA','CX','CY','CZ'
      IF(DIBUJO.EQ.'DIB')WRITE(15,'(/,1X,A,1X,A,6X,3(A,15X),/
     &,10X,9(A,7X))')
     &'ESFERA',' STY','VERTICE 1','VERTICE 2','VERTICE 3',
     &'x','y','z','x','y','z','x','y','z'
C It selects one sphere
      DO 1 I=1,NESF
      NPE(I)=0
      SSFE(I)=ZERO
      IF(NOUSA(I))GO TO 1
      REI=RE(I)
      SRE(I)=FOUR*PI*REI*REI
      ATS=SRE(I)/FNDIV
      REIM=REI-DVEC
      FDR=REIM/REI
      XEI=XE(I)
      YEI=YE(I)
      ZEI=ZE(I)
      DDCO=XEI*XEI+YEI*YEI+ZEI*ZEI
      VCTEN=REI*REI-DDCO
      VCTED=6.0D0*REI
      IIJ=0
C
C  It fixes which spheres are linked to sphere I
      DO 500 J=1,NESF
      IF(NOUSA(J))GO TO 500
      IF(I.EQ.J) GO TO 500
      DIJ2=(XEI-XE(J))*(XEI-XE(J))+
     &     (YEI-YE(J))*(YEI-YE(J))+
     &     (ZEI-ZE(J))*(ZEI-ZE(J))
      SRE2=(REI+RE(J))*(REI+RE(J))
      IF(DIJ2.GT.SRE2) GO TO 500
      SRE2=(RE(J)-REI)*(RE(J)-REI)
      IF(DIJ2.GT.SRE2)GO TO 501
      sre1=RE(J)-REI
      SRE2=sre1*sre1
      IF(DIJ2.GT.SRE2)GO TO 501
      if(sre1.lt.0.0d0)go to 502
      NOUSA(I)=.TRUE. 
      NESI=NESI+1
      GO TO 1
  502 continue
      nousa(j)=.true.
      nesi=nesi+1
      go to 500
  501 continue
      IIJ=IIJ+1
      IJE(IIJ)=J
      NEJCI=IIJ
  500 CONTINUE 
C 
C It selects one main triangle.
      NSUP=0
      DO 2 J=1,60
      XPL=0.0D0
      YPL=0.0D0
      ZPL=0.0D0
      NTS=0
      NINF=NSUP+1
      NSUP=NINF+NTRIAN-1
C
C It selects one secondary triangle.
      DO 3 K=NINF,NSUP
      XSL=XC1(K)*REI
      YSL=YC1(K)*REI
      ZSL=ZC1(K)*REI
      XSM=XSL+XEI
      YSM=YSL+YEI
      ZSM=ZSL+ZEI
C
C It fixes if the secundary triangle is inside or outside.
      DO 8 N3=1,NEJCI
      N4=IJE(N3)
      DD=(XSM-XE(N4))*(XSM-XE(N4))+
     &   (YSM-YE(N4))*(YSM-YE(N4))+
     &   (ZSM-ZE(N4))*(ZSM-ZE(N4))
      RREJ=RE(N4)*RE(N4)
      IF(DD.LT.RREJ) GO TO 3
    8 CONTINUE
C
C It prepares the coordinates to the main triangle
      XPL=XPL+XSL
      YPL=YPL+YSL
      ZPL=ZPL+ZSL
      NTS=NTS+1
    3 CONTINUE
C
C It reduces the secondary triangles to the main triangle.
      IF(NTS.EQ.0)GO TO 2
      ATP=ATS*NTS
      XPL=XPL/FLOAT(NTS)
      YPL=YPL/FLOAT(NTS)
      ZPL=ZPL/FLOAT(NTS)
      RRR=SQRT(XPL*XPL+YPL*YPL+ZPL*ZPL)
      FC=REI/RRR
      XPL=XPL*FC
      YPL=YPL*FC
      ZPL=ZPL*FC
      XPM=XPL+XEI
      YPM=YPL+YEI
      ZPM=ZPL+ZEI
      IJ=IJ+1
      SSFE(I)=SSFE(I)+ATP
      NPE(I)=NPE(I)+1
      DD=XPM*XPM+YPM*YPM+ZPM*ZPM
      VOL=VOL+ATP*(DD+VCTEN)/(VCTED)
      IF(VECTOR.EQ.'VEN') GO TO 2
C
C   It computes the vector whith module equal to DVEC
C   
      XV1=XPL*FDR+XEI
      YV1=YPL*FDR+YEI
      ZV1=ZPL*FDR+ZEI
      XVEC=XV1-XPM
      YVEC=YV1-YPM
      ZVEC=ZV1-ZPM
      WRITE(8,'(I6,I3,7F11.6)')I,ISPHERETYPE(I),
     &XPM,YPM,ZPM,ATP,XVEC,YVEC,ZVEC
C
C     	AQUI CALCULAMOS LOS VERTICES DE LOS TRIANGULOS PRINCIPALES
C	DE CADA ESFERA QUE PERTENECEN A LA SUPERFICIE MOLECULAR
C
      IF(DIBUJO.NE.'DIB') GOTO 2
      NVX1=JVT1(1,J)
      NVX2=JVT1(2,J)
      NVX3=JVT1(3,J)
      DO 77 FF=1,3
	  VERTIX(1,FF)=CV(NVX1,FF)*REI
	  VERTIX(2,FF)=CV(NVX2,FF)*REI
	  VERTIX(3,FF)=CV(NVX3,FF)*REI
   77 CONTINUE
      DO 78 GG=1,3
          VERTIX(GG,1)=VERTIX(GG,1)+XEI
          VERTIX(GG,2)=VERTIX(GG,2)+YEI
          VERTIX(GG,3)=VERTIX(GG,3)+ZEI
   78 CONTINUE
      WRITE(15,'(I6,I3,9(F7.3,X))')I,ISPHERETYPE(I) 
     &,((VERTIX(SS,TT),TT=1,3),SS=1,3)   

    2 CONTINUE
      STOT=STOT+SSFE(I)
      IF(NPE(I).EQ.0)NES0=NES0+1
    1 CONTINUE
C***** Print *******
      WRITE(6,'(A,I8/)')
     &' Number of spheres with area equal to 0 =',NES0
      WRITE(6,'(2(/A,F17.6)/A,I8)')' Area=',STOT,' Volume=',VOL,
     &      ' Number of points=',IJ
      IF(OPESF.EQ.'SFN') GO TO 3333
      WRITE(7,'(//A,3(11X,A),3(7X,A))')
     & '     SN    AREA','X','Y','Z','R','AREA(R)','NP'
      DO 601 I=1,NESF
      MARCA='      '
      IF(NPE(I).EQ.0)MARCA='*CERO*'
      IF(NOUSA(I))MARCA='NOUSED'
      WRITE(7,'(I7,F12.7,3F12.4,F8.4,F12.7,3I7,A)')
     &      I,SSFE(I),XE(I),YE(I),ZE(I),RE(I),SRE(I),
     &      NPE(I),NPEC(I,1),NPEC(I,2),MARCA
  601 CONTINUE
 3333 CONTINUE
      RETURN
      END
C
      SUBROUTINE GWRITE(IUNIT)
C --------------------------------------------------------------
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*4       (I-N) 
      CHARACTER*3 SUAVE,vector,suma,OPESF,DIBUJO
      CHARACTER*80 TITULO
      LOGICAL NOUSA
C      PARAMETER nst=100000
C      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(NST)
      COMMON/SFE1/OMEGA,RD,RET,FRO,NESF,NDIV,DVEC,NESFI,NOUSA(100000)
      COMMON/ESC/TITULO,SUAVE,vector,SUMA,OPESF,FRADIO
     &,DIBUJO
      NST=100000
      WRITE(IUNIT,'(A)') ' GEPOL/87'
      WRITE(IUNIT,'(A,A)')' REMARKS ',TITULO
      WRITE(IUNIT,'(6(A,2X)/A,2X,A,I1,A,F10.5/A,3(2X,A,F10.5))')
     & ' REMARKS',SUAVE,vector,suma,OPESF,DIBUJO,
     & ' REMARKS',' NDIV=',NDIV,'   RD=',RD,
     & ' REMARKS',' OMEGA=',OMEGA,'FRADIO=',FRADIO
      RETURN
      END

