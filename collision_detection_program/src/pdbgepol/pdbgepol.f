CCCCCC It creates a gepol input from PDB files
CCCCCC F.Villar & I.Tunon 26-XI-1989
CCCCCC Departament of Physical-Chemistry
CCCCCC University of Valencia (Spain)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 TITU,SOBRA*4,NAME*4,NOMBRE*1
      character*80 filename,fileoutput*80 
C      PARAMETER NDIV=3,RD=1.5,FRADIO=.7
        NDIV=3
        RD=1.5
        FRADIO=.7
      ZERO=0.0 
      RC=1.0

CCCCCC It asks for specifications
    
      write (6,*)'natoms?'
      read(5,*) natoms
      write (6,*)' primer residuo?'
      read (5,*) imin
      write (6,*)' ultimo residuo?'
      read (5,*) imax 
      write(6,*)'filename input?'
      read(5,'(A)') filename 
      write (6,*)'filename output?'
      read (5,'(A)') fileoutput
      write (6,*)'constante para radio?'
      read (5,*) RC
      if (RC.eq.ZERO) RC=1.0


      open(15,status='old',file=filename)
      open(16,status='new',file=fileoutput) 


      READ(15,'(A)')TITU
      WRITE(16,'(A)')TITU
      WRITE(16,'(A)')'CON,VES,NSU,SFY,DIB'
      WRITE(16,'(I2,2F10.5)')NDIV,RD,FRADIO
      WRITE (16,'(I8)') natoms
      INUM=0 
      ius=0


      DO I=1,natoms
        ius=0
        READ(15,'(A,7X,A,9X,I3,4X,3F8.4)')SOBRA,NAME,IRES,X,Y,Z
        IF (SOBRA.NE.'ATOM') GOTO 10
        IF ((IRES.GE.IMIN).AND.(IRES.LE.IMAX)) IUS=1
        RADIO=1.8000*RC
        itipo=4
        do J=1,4
          IF (NAME(J:J).EQ.' ') GOTO 55
          NOMBRE=NAME(J:J)
          GOTO 65
  55      CONTINUE
        END DO   
  65    CONTINUE
   
        IF (NOMBRE.EQ.'C') RADIO=1.6000*RC
        IF (NOMBRE.EQ.'N') RADIO=1.5000*RC
        IF (NOMBRE.EQ.'H') RADIO=1.2000*RC
        IF (NOMBRE.EQ.'O') RADIO=1.4000*RC
        IF (NOMBRE.EQ.'S') RADIO=1.6000*RC
        IF (NOMBRE.EQ.'C') itipo=1
        IF (NOMBRE.EQ.'N') itipo=2
        IF (NOMBRE.EQ.'H') itipo=4
        IF (NOMBRE.EQ.'O') itipo=3
        IF (NOMBRE.EQ.'S') itipo=4
        WRITE(16,'(4F10.5,I2,I2)')X,Y,Z,RADIO,IUS,itipo
        INUMBER=INUMBER+1
  10    CONTINUE
      END DO
      
        
      close(15)
      close(16) 
      
        
      STOP 
      END
