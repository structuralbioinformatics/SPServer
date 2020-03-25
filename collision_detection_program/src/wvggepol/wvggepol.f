        program wvggepol

C*****************************************************
C
C     (c)  F.Villar & I.Tunon 20-II-1990
C      Departament of Physical-Chemistry
C      University of Valencia (Spain)
C
C*****************************************************
C      It creates a gepol input from WVG files
C ----------------------------------------------------
C      
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 TITU,SOBRA*4,NAME*4,NOMBRE*1
      character*80 filename,fileinp*80 ,filevel*80 
      PARAMETER NDIV=3,RD=0.15,FRADIO=.6, IUS=1,ZERO=0.0
       
      RC=1.0

C  It asks for specifications

  33  format ( /A, $)  
      write (6,33)' first residue number ?  '
      read (5,*) imin
      write (6,33)' last residue number ?  '
      read (5,*) imax
      write (6,33)' number of atoms in those residues ?  '
      read (5, *) numero 
      write(6,33)'input filename ?  '
      read(5,'(A)') filename 
      write (6,33)'enter name of output file .inp :  '
      read (5,'(A)') fileinp
      write (6,33)'enter name of output file .vel :  '
      read (5,'(A)') filevel
      write (6,33)'atomic radio scale factor (1 if nm) ?  '
      read (5,*) RC
      if (RC.eq.ZERO) RC=1.0

      open(15,form='formatted',status='OLD',file=filename)
      open(16,form='formatted',status='NEW',file=fileinp) 
      open(17,form='formatted',status='NEW',file=filevel) 


      READ(15,'(A)')TITU
      read(15,'(I5)') natoms
      WRITE(16,'(A)')TITU
      WRITE(16,'(A)')'CON VEN NSU SFN DIB'
      
      WRITE(17,'(A)')TITU
      WRITE(17,'(A)')'CON VEN NSU SFN DIB'

      k=0
      DO I=1,natoms
        READ(15,'(1X,I4,6X,A,5X,3F8.3,3F8.4)')IRES,NAME,
     &                        X, Y, Z, velx, vely, velz    
        IF ((IRES.GE.IMIN).AND.(IRES.LE.IMAX)) THEN
          
          k=k+1
          if (k.eq.1) then
             inicio=I
             WRITE(16,'(I2,2F10.5,1X,I1,1X,I5)')NDIV,RD*RC
     &            ,FRADIO,1,inicio
             WRITE (16,'(I8)') numero
             WRITE(17,'(I2,2F10.5,1X,I1,1X,I5)')NDIV,RD*RC
     &            ,FRADIO,1,inicio
             WRITE (17,'(I8)') numero
           end if   
             
          RADIO=1.8000*RC
          do J=1,4
            IF (NAME(J:J).EQ.' ') GOTO 55
            NOMBRE=NAME(J:J)
            GOTO 65
  55        CONTINUE
          END DO   
  65      CONTINUE
   
          IF (NOMBRE.EQ.'C') then
              RADIO=0.18000*RC
              natomico=6
              goto 100
          end if    
          IF (NOMBRE.EQ.'N') then
              RADIO=0.15000*RC
              natomico=7
              goto 100
          end if    
          IF (NOMBRE.EQ.'H') then
              RADIO=0.12000*RC
              natomico=1
              goto 100
          end if    
          IF (NOMBRE.EQ.'O') then
              RADIO=0.14000*RC
              natomico=8
              goto 100
          end if    
          IF (NOMBRE.EQ.'S') then 
              RADIO=0.16000*RC
              natomico=16
              goto 100
          end if    
          IF (NOMBRE.EQ.'P') then 
              RADIO=0.17000*RC
              natomico=15
              goto 100
          end if
          goto 200
          
  100     continue       

          WRITE(16,'(4F10.5,I2,6X,I4)') X,Y,Z,RADIO,IUS,natomico
          WRITE(17,'(1X,I5,3(2X,F8.4))') i, velx, vely, velz
        end if  
      END DO
      goto 300
  200 continue    
      write (6, *)' ## Sorry, radio data not available for atom'
     &,' type "',NAME,'"'
  300 continue    
        
      close(15)
      close(16) 
      close(17)
        
      STOP 
      END
