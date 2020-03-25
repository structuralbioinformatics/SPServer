/***  REMOVE ANY LINE ABOVE  *THIS*  LINE USING A STANDARD TEXT EDITOR !!! 

      DO NOT EXTRACT PART OF THIS FILE, E.G. PUTATIVE INCLUDE FILES - 
      YOU HAVE TO DEAL WITH *ONE* FILE ONLY.

      SAVE THIS ENTIRE FILE UNDER THE NAME "dssp.c"  !!!

      COMPILE THE PROGRAM ON SUN's WITH "cc dssp.c -o dssp -lm"
      AND ON IRIS'es WITH "cc -cckr dssp.c -o dssp -lm"

      USAGE: dssp [-na] pdb_file dssp_file
          the -na flag disables the calculation of accessible surface

***/





/************************************************************
*                                                           *
*          dssp.c                                           *
*                                                           *
************************************************************/


/* Output from p2c, the Pascal-to-C translator */
/* From input file "dssp.p" */



/* ------------------------------------------------------------------ */
/*

DSSP version October 1988.
This copy for Baldomero Miguel at Univ_Barcelona
who have agreed to the following software license agreement:


An academic license for the DSSP program
((c) W. Kabsch, C. Sander and MPI-MF, 1983, 1985, 1988)
is granted to in exchange for the following commitments:

I hereby certify that

        (1) I am an academic user at an academic research institution. In
            using the software, I will respect the interests of the authors
            and their institutions.

        (2) I will not use the software in commercial activities without
            a written commercial license agreement; commercial activities
            include, in particular, work under contract from a commercial
            company.

        (3) I will not redistribute the software to others outside of my
            immediate research group. I will suggest to other interested
            research groups to contact the authors directly.

        (4) I will not alter or suppress the run-time copyright message.

        (5) I will acknowledge the program authors on any publication of
            scientific results based in part on use of the program and
            cite the article in which the program was described.

        (6) I will report evidence of program bugs to the authors.

        (7) I will send the source code of any bug corrections and program
            extensions, major or minor, to the original authors, for free
            academic use. If I have made major extensions which are incor-
            porated by the authors, I reserve the right to be appropriately
            included in any future commercial license agreement.

        (8) I will not extract part of the software, e.g. modules or sub-
            routines, for use in other contexts without permission by the
            authors.

        (9) I will not use the program in the context of classified research.
*/
/* PREVIOUS RELEASE: VERSION OCTOBER 1985                             */
/* PREVIOUS RELEASE: VERSION JUNE 1983                                */
/* LANGUAGE: STANDARD PASCAL WITH 128 CHARACTER ASCII SET             */
/* AUTHORS AND COPYRIGHT (1983,1985,1988):
   Wolfgang Kabsch and Chris Sander, Max Planck Institut
   fuer Medizinische Forschung, Jahnstr. 29, 6900 Heidelberg, Germany
   Telephone: +49-6221-486 276  Telex: 461505 mpimf d
   Bitnet:    KABSCH@EMBL
   Current address for Chris Sander:
   Biocomputing, EMBL, 6900 Heidelberg, Germany
   Telephone: +49-6221-387 361 Telex: 461613 embl d
   Telefax:   +49-6221-387 306
   Bitnet:    SANDER@EMBL
   Do report errors if you find any.
   Reference: Kabsch,W. and Sander,C. (1983) Biopolymers 22, 2577-2637*/
/*--------------------------------------------------------------------*/
/* DEFINES SECONDARY STRUCTURE AND SOLVENT EXPOSURE OF PROTEINS FROM
   ATOMIC COORDINATES AS GIVEN BY THE BROOKHAVEN PROTEIN DATA BANK.   */
/*--------------------------------------------------------------------*/
/* This program including sample input and output files for dataset 1PPT
   is available from the authors in exchange for an academic or
   commercial license agreement. The program is no longer available
   from the Brookhaven Protein Data Bank */
/*--------------------------------------------------------------------*/
/* CORRECTION AND MODIFICATION LOG SINCE JUNE 1983 */
/* (1) MODIFICATIONS THAT AFFECT OUTPUT ON FILE TAPEOUT FOR AT LEAST ONE
       OF THE 62 PROTEIN DATA SETS IN THE 1983 BIOPOLYMERS PAPER:
   - SIDECHAIN ATOMS MORE THAN MAXDIST ANGSTROM DISTANT FROM ATOM CA ARE
     DECLARED ILLEGAL AND IGNORED. OUTPUT CHANGE: ACCESSIBILITY VALUES
     FOR ASN 76 OF 1SBT (ILLEGAL ATOM OD1) AND PRO 49 OF 156B (ILLEGAL
     ATOM UNK).
   - ANY RESIDUE WITH INCOMPLETE BACKBONE IS IGNORED. OUTPUT CHANGE:
     CHAIN BREAK BETWEEN RESIDUE SER 11 AND ILE 16 IN 2GCH
     (DUE TO INCOMPLETE COORDINATES FOR SER 11) IS NOW CORRECT.
   (2) MODIFICATIONS THAT DO NOT AFFECT OUTPUT ON FILE TAPEOUT FOR ANY
       OF THE 62 PROTEIN DATA SETS IN THE 1983 BIOPOLYMERS PAPER:
   - SPELLING OF FLAGCHIRALITY AND TESTSSBOND CORRECTED.
   - WARNING MESSAGE FOR RESIDUES WITH NON-STANDARD NUMBER OF
     SIDECHAIN ATOMS HAS BEEN ADDED. FOR EXAMPLE, THIS ALERTS THE USER
     TO BAD DATA FOR RESIDUES 8,12,21,24 AND 44 OF DATA SET 2RXN.
   - WARNING MESSAGE FOR RESIDUES IGNORED DUE TO NON-STANDARD RESIDUE
     NAME SUCH AS 'ACE' AND 'FOR' HAS BEEN ADDED.
   - WARNING MESSAGE FOR ALTERNATE ATOM LOCATION IDENTIFIER HAS BEEN
     ADDED. FOR EXAMPLE, THE USER IS NOW WARNED THAT ATOM CH2 IN ANY
     TRP OF DATA SET 1APP IS IGNORED DUE TO BAD ALTERNATE LOCATION
     IDENTIFIER '2'.
   WE THANK STEVEN SHERIFF, FRANCOIS COLONNA AND JOHN MOULT FOR
   REPORTING PROBLEMS AND STEVEN SHERIF, JANET THORNTON AND
   WILLIE TAYLOR FOR RESULTS OF TEST RUNS ON VAX COMPUTERS.

   Changes after 1985:

   - program speeded up by a factor of two or three by avoiding use
     of square root.
   - hydrogen atoms in data set ignored on input (H of NH of backbone
     is built as before)
   - 18-AUG-1988: CADIST=9.0, replacing CADIST=8.0. Has affected output
     for 63/300 proteins in a minor way. Thanks to Jean Richelle (Bruxelles)
     for pointing out this bug.

     Output changes due to change in parameter CADIST (8 to 9 Angstrom) :
     additional backbone-backbone Hbonds found with slight
     adjustments in secondary structure summary. In about 300 protein
     data sets from the Fall 1988 PDB release, 63 additional
     Hbonds were found, i.e. 0.2 Hbonds per protein (29 of type
     i,i+5;  16 of type i,i+4; 6 of type i,i+3; 10 in antiparallel beta
     bridges and 2 in a parallel beta bridge). These additional
     Hbonds affected the secondary structure summary of 26 of these
     protein data sets in a minor way, typically joining a 4-turn to
     an alpha-helix, changing a geometrical turn to a hydrogen-
     bonded turn or adding an extra residue pair to a beta ladder.
     The changes are (using _ for blank):

     [protein id, old secstruc summary > corrected summary]

     1FC2_E > EE   and  _T > ET
     1GP1GGG > HHH
     1HBSS > T
     1HDSS > T and  GGGGGG > TTHHHH
     1HFM__ > TT
     1HKGSSS > TTT
     1IG2S_ > TT
     1LDX GGG > HTT
     1MEV__ > TT  and  _BS > TBS  and  SSS > TTS
     1PFCSSS > TTS
     1PP2_E > EE  and  _S > ES
     1RN3E_SS_E > EEEEEE  and _E > EE  (>3-res beta bulge)
     1RNSsame as 1RN3
     2ATCHH > TT
     2CABB_ > EE
     2CPPSS > TT  and  GGGGGG > HHHHTT
     2LYZT > H
     2MDHSSS > TTT
     3CPA TTT > HHH
     4CATTTT > HHH
     4SBVS > T
     5API_ > B
     5CPATTT > HHH
     7LYZS > H
     8CAT_ > B  and  _ > B
     8LYZT > H

     Note that this bugfix results in a small variation in the total
     number of Hbonds, compared to the variation which would
     result, say, from changing the (somewhat arbitrary) cutoff of
     -0.5 kcal/mol for the Hbond electrostatic potential energy. We
     cannot here solve the fundamental difficulty of arbitrary
     cutoffs involved in extracting binary classifications (an Hbond
     exists, yes/no) from real numbers (coordinates, energies).
     However, for most purposes the secondary structure summary agrees
     will with anyone's intuitive definition, especially for well-refined and
     high resolution structures. For a more clearcut assignment of protein
     substructure, we recommend using the detailed H-bond and other assignments
     in the columns following the summary column, i.e. columns 19-38 (xxx):

     ....;....1....;....2....;....3....;....4....;....5....;....6....;....7..
                       xxxxxxxxxxxxxxxxxxxx
                       .-- 3-turns/helix
                       |.-- 4-turns/helix
                       ||.-- 5-turns/helix
                       |||.-- geometrical bend
                       ||||.-- chirality
                       |||||.-- beta bridge label
                       ||||||.-- beta bridge label
                       |||||||   .-- beta bridge partner resnum
                       |||||||   |   .-- beta bridge partner resnum
                       |||||||   |   |.-- beta sheet label
                       |||||||   |   ||   .-- solvent accessibility
                       |||||||   |   ||   |
        35   47   I  E     +     0   0    2
        36   48   R  E >  S- K   0  39C  97
        37   49   Q  T 3  S+     0   0   86    (example from 1EST)
        38   50   N  T 3  S+     0   0   34
        39   51   W  E <   -KL  36  98C   6
                                                                           */
/*--------------------------------------------------------------------*/
/* GENERAL PROGRAM INSTALLATION GUIDE. */
/* (1) THE PROGRAM REQUIRES THE FULL STANDARD ASCII 128 CHARACTER SET,
       IN PARTICULAR LOWER CASE LETTERS 'abcdefg....'.
   (2) STANDARD PASCAL MAY NOT RECOGNIZE REAL NUMBERS SUCH AS .1, +.1,
       -.1 ON INPUT. CHANGE TO 0.1,+0.1,-0.1.
   (3) THE NON-STANDARD PROCEDURE 'DATE' RETURNS THE CURRENT DAY, MONTH,
       AND YEAR. IF THE PROCEDURE IS NOT CALLED (LINE COMMENTED OUT)
       THE PSYCHEDELIC DATE DEC 24, 2001 IS RETURNED. YOU MAY  REPLACE
       'DATE' BY THE CORRESPONDING PROCEDURE FROM YOUR PASCAL
       IMPLEMENTATION. THE EXAMPLE GIVEN WORKS IN DEC VAX VMS 5.0.
   (4) DUE TO INCOMPATIBLE ASCII CODES, SQUARE BRACKETS '[' AND ']'
       MAY APPEAR AS '!','?' ETC. USE YOUR EDITOR TO CONVERT THESE.   */
/* INSTALLATION GUIDE FOR VAX/VMS USERS. */
/* (1) THE /OPTIMIZE OPTION OF THE PASCAL COMPILER PRODUCED
       INCORRECT CODE ON THE VAX 8600 AT EMBL RUNNING UNDER VMS V4.2.
       LATER VERSIONS OF VMS (E.G. VMS 5.0) PRODUCED CORRECT CODE.
       IF IN DOUBT, COMPILE USING PASCAL /NOOPTIMIZE.
   (2) COPY BROOKHAVEN DATA BANK COORDINATE INPUT TO A FILE NAMED
       TAPEIN.DAT . OUTPUT WILL BE IN A FILE NAMED TAPEOUT.DAT        */
/* IMPLEMENTATION ON OTHER COMPUTERS */
/* (1) NORD-500. EXECUTION TIME COMPARABLE TO VAX 780.
   (2) SUN-3.    EXECUTION TIME COMPARABLE TO VAX 780.
                 Compile using: pc -L
                 in ORDER to map upper case letters in keywords
                 and identifiers to lower case.
   (3) ATARI 520 ST. RUNS FACTOR 60 SLOWER THAN NORD-500 DUE TO
       SOFTWARE-EMULATED FLOATING POINT OPERATIONS ON MC68000.        */
/*--------------------------------------------------------------------*/
/* INPUT/OUTPUT FILES. */
/* INPUT:   DEFAULT  INPUT UNIT, E.G. YOUR TERMINAL
   OUTPUT:  DEFAULT OUTPUT UNIT, E.G. YOUR TERMINAL,
            USED FOR RUN-TIME MESSAGES. WARNINGS AND ERRORS LOOK
            LIKE THIS: !!! TEXT !!!
   TAPEIN:  FILE WITH PROTEIN DATA BANK COORDINATES, E.G. PDB3PTI.COO
   TAPEOUT: DSSP OUTPUT OF LINE LENGTH 128, E.G. PAPER PRINTER        */
/*--------------------------------------------------------------------*/
/* DESCRIPTION OF OUTPUT ON FILE TAPEOUT:
   LINE LENGTH OF OUTPUT IS 128 CHARCTERS.
   FOR DEFINITONS, SEE ABOVE BIOPOLYMERS ARTICLE.
   IN ADDITION NOTE:
   HISTOGRAMS - E.G. 2 UNDER COLUMN '8' IN LINE 'RESIDUES PER ALPHA
            HELIX' MEANS: THERE ARE 2 ALPHA HELICES OF LENGTH  8
            RESIDUES IN THIS DATA SET.
   #  RESIDUE AA STRUCTURE BP1 BP2 ACC ..ETC..FOR EACH RESIDUE I:
   #  RESIDUE - TWO COLUMNS OF RESIDUE NUMBERS. FIRST COLUMN IS DSSP'S
            SEQUENTIAL RESIDUE NUMBER, STARTING AT THE FIRST
            RESIDUE ACTUALLY IN THE DATA SET AND INCLUDING CHAIN BREAKS;
            THIS NUMBER IS USED TO REFER TO RESIDUES THROUGHOUT. SECOND
            COLUMN GIVES CRYSTALLOGRAPHERS' 'RESIDUE SEQUENCE
            NUMBER','INSERTION CODE' AND 'CHAIN IDENTIFIER' (SEE PROTEIN
            DATA BANK FILE RECORD FORMAT MANUAL), GIVEN FOR REFERENCE
            ONLY AND NOT USED FURTHER..
   AA -     ONE LETTER AMINO ACID CODE, LOWER CASE FOR SS-BRIDGE CYS.
   STRUCTURE - SEE BIOPOLYMERS
   BP1 BP2  - RESIDUE NUMBER OF FIRST AND SECOND BRIDGE PARTNER
            FOLLOWED BY ONE LETTER SHEET LABEL
   ACC -    NUMBER OF WATER MOLECULES IN CONTACT WITH THIS RESIDUE *10.
            OR RESIDUE WATER EXPOSED SURFACE IN ANGSTROM**2.
   N-H-->O ETC. -  HYDROGEN BONDS. E.G. -3,-1.4 MEANS: IF THIS RESIDUE
            IS RESIDUE I THEN N-H OF I IS H-BONDED TO C=O OF I-3
            WITH AN ELECTROSTATIC H-BOND ENERGY OF -1.4 KCAL/MOL.
   TCO -    COSINE OF ANGLE BETWEEN C=O OF RESIDUE I AND C=O OF
            RESIDUE I-1. FOR ALPHA-HELICES, TCO IS NEAR +1, FOR
            BETA-SHEETS TCO IS NEAR -1. NOT USED FOR STRUCTURE
            DEFINITION.
   KAPPA -  VIRTUAL BOND ANGLE (BEND ANGLE) DEFINED BY THE THREE
            C-ALPHA ATOMS OF RESIDUES I-2,I,I+2. USED TO DEFINE
            BEND (STRUCTURE CODE 'S').
   ALPHA -  VIRTUAL TORSION ANGLE (DIHEDRAL ANGLE) DEFINED BY THE FOUR
            C-ALPHA ATOMS OF RESIDUES I-1,I,I+1,I+2. USED TO DEFINE
            CHIRALITY (STRUCTURE CODE '+' OR '-').
   PHI PSI - IUPAC PEPTIDE BACKBONE TORSION ANGLES
   X-CA Y-CA Z-CA -  ECHO OF C-ALPHA ATOM COORDINATES              */
/*--------------------------------------------------------------------*/
/* WORDS OF CAUTION */
/* THE VALUES FOR SOLVENT EXPOSURE MAY NOT MEAN WHAT YOU THINK!
    (A) EFFECTS LEADING TO LARGER THAN EXPECTED VALUES:
     SOLVENT EXPOSURE CALCULATION IGNORES UNUSUAL RESIDUES, LIKE ACE,
     OR RESIDUES WITH INCOMPLETE BACKBONE, LIKE ALA 1 OF DATA SET 1CPA.
     IT ALSO IGNORES HETATOMS, LIKE A HEME OR METAL LIGANDS.
     ALSO, SIDE CHAINS MAY BE INCOMPLETE (AN ERROR MESSAGE IS WRITTEN).
    (B) EFFECTS LEADING TO SMALLER THAN EXPECTED VALUES:
     IF YOU APPLY THIS PROGRAM TO PROTEIN DATA BANK DATA SETS
     CONTAINING OLIGOMERS, SOLVENT EXPOSURE IS FOR THE ENTIRE ASSEMBLY,
     NOT FOR THE MONOMER. ALSO, ATOM OXT OF C-TERMINAL RESIDUES IS
     TREATED LIKE A SIDE CHAIN ATOM IF IT IS LISTED AS PART OF THE LAST
     RESIDUE. ALSO, PEPTIDE SUBSTRATES, WHEN LISTED AS ATOMS RATHER THAN
     HETATOMS, ARE TREATED AS PART OF THE PROTEIN, E.G. RESIDUES 499 S
     AND 500 S IN 1CPA.                                               */
/* UNKNOWN OR UNUSUAL RESIDUES ARE NAMED X ON OUTPUT AND THEY ARE
   NOT CHECKED FOR STANDARD NUMBER OF SIDECHAIN ATOMS.                */
/* ALL EXPLICIT WATER MOLECULES, LIKE OTHER HETATOMS, ARE IGNORED.    */
/* END OF INTRODUCTORY COMMENTS */
/**********************************************************************/


/************************************************************
*                                                           *
*          p2c.h                                            *
*                                                           *
************************************************************/


#ifndef P2C_H
#define P2C_H


/* Header file for code generated by "p2c", the Pascal-to-C translator */

/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version 1.19.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */


#include <stdio.h>



/* If the following heuristic fails, compile -DBSD=0 for non-BSD systems,
   or -DBSD=1 for BSD systems. */

#ifdef M_XENIX
# define BSD 0
#endif

#ifdef vms
# define BSD 0
# ifndef __STDC__
#  define __STDC__ 1
# endif
#endif

#ifdef __TURBOC__
# define MSDOS 1
#endif

#ifdef MSDOS
# define BSD 0
#endif

#ifdef FILE       /* a #define in BSD, a typedef in SYSV (hp-ux, at least) */
# ifndef BSD	  /*  (a convenient, but horrible kludge!) */
#  define BSD 1
# endif
#endif

#ifdef BSD
# if !BSD
#  undef BSD
# endif
#endif


#ifdef __STDC__
# include <stddef.h>
# include <stdlib.h>
# define HAS_STDLIB
# ifdef vms
#  define __ID__(a)a
# endif
#else
# ifndef BSD
#  ifndef __TURBOC__
#   include <memory.h>
#  endif
# endif
# ifdef hpux
#  ifdef _INCLUDE__STDC__
#   include <stddef.h>
#   include <stdlib.h>
#  endif
# endif
# include <sys/types.h>
# if !defined(MSDOS) || defined(__TURBOC__)
#  define __ID__(a)a
# endif
#endif

#ifdef __ID__
# define __CAT__(a,b)__ID__(a)b
#else
# define __CAT__(a,b)a##b
#endif


#ifdef BSD
# include <strings.h>
# define memcpy(a,b,n) (bcopy(b,a,n),a)
# define memcmp(a,b,n) bcmp(a,b,n)
# define strchr(s,c) index(s,c)
# define strrchr(s,c) rindex(s,c)
#else
# include <string.h>
#endif

#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <assert.h>


#ifdef vms

#define LACK_LABS
#define LACK_MEMMOVE
#define LACK_MEMCPY

#else

#define LACK_LABS       /* Undefine these if your library has these */
#define LACK_MEMMOVE

#endif


typedef struct __p2c_jmp_buf {
    struct __p2c_jmp_buf *next;
    jmp_buf jbuf;
} __p2c_jmp_buf;


/* Warning: The following will not work if setjmp is used simultaneously.
   This also violates the ANSI restriction about using vars after longjmp,
   but a typical implementation of longjmp will get it right anyway. */

#ifndef FAKE_TRY
# define TRY(x)         do { __p2c_jmp_buf __try_jb;  \
			     __try_jb.next = __top_jb;  \
			     if (!setjmp((__top_jb = &__try_jb)->jbuf)) {
# define RECOVER(x)	__top_jb = __try_jb.next; } else {
# define RECOVER2(x,L)  __top_jb = __try_jb.next; } else {  \
			     if (0) { L: __top_jb = __try_jb.next; }
# define ENDTRY(x)      } } while (0) 
#else
# define TRY(x)         if (1) {
# define RECOVER(x)     } else do {
# define RECOVER2(x,L)  } else do { L: ;
# define ENDTRY(x)      } while (0)
#endif



#ifdef M_XENIX  /* avoid compiler bug */
# define SHORT_MAX  (32767)
# define SHORT_MIN  (-32768)
#endif


/* The following definitions work only on twos-complement machines */
#ifndef SHORT_MAX
# define SHORT_MAX  ((short)(((unsigned short) -1) >> 1))
# define SHORT_MIN  (~SHORT_MAX)
#endif

#ifndef INT_MAX
# define INT_MAX    ((int)(((unsigned int) -1) >> 1))
# define INT_MIN    (~INT_MAX)
#endif

#ifndef LONG_MAX
# define LONG_MAX   ((long)(((unsigned long) -1) >> 1))
# define LONG_MIN   (~LONG_MAX)
#endif

#ifndef SEEK_SET
# define SEEK_SET   0
# define SEEK_CUR   1
# define SEEK_END   2
#endif

#ifndef EXIT_SUCCESS
# ifdef vms
#  define EXIT_SUCCESS  1
#  define EXIT_FAILURE  (02000000000L)
# else
#  define EXIT_SUCCESS  0
#  define EXIT_FAILURE  1
# endif
#endif


#define SETBITS  32


#ifdef __STDC__
# ifndef vms
#  define Signed    signed
# else
#  define Signed
# endif
# define Void       void      /* Void f() = procedure */
# ifndef Const
#  define Const     const
# endif
# ifndef Volatile
# define Volatile  volatile
# endif
# define PP(x)      x         /* function prototype */
# define PV()       (void)    /* null function prototype */
typedef void *Anyptr;
#else
# define Signed
# define Void       void
# ifndef Const
#  define Const
# endif
# ifndef Volatile
#  define Volatile
# endif
# define PP(x)      ()
# define PV()       ()
typedef char *Anyptr;
#endif

#ifdef __GNUC__
# define Inline     inline
#else
# define Inline
#endif

#define Register    register  /* Register variables */
#define Char        char      /* Characters (not bytes) */

#ifndef Static
# define Static     static    /* Private global funcs and vars */
#endif

#ifndef Local
# define Local      static    /* Nested functions */
#endif

typedef Signed   char schar;
typedef unsigned char uchar;
typedef unsigned char boolean;

#ifndef true
# define true    1
# define false   0
#endif


typedef struct {
    Anyptr proc, link;
} _PROCEDURE;

#ifndef _FNSIZE
# define _FNSIZE  120
#endif


extern Void    PASCAL_MAIN  PP( (int, Char **) );
extern Char    **P_argv;
extern int     P_argc;
extern short   P_escapecode;
extern int     P_ioresult;
extern __p2c_jmp_buf *__top_jb;


#ifdef P2C_H_PROTO   /* if you have Ansi C but non-prototyped header files */
extern Char    *strcat      PP( (Char *, Const Char *) );
extern Char    *strchr      PP( (Const Char *, int) );
extern int      strcmp      PP( (Const Char *, Const Char *) );
extern Char    *strcpy      PP( (Char *, Const Char *) );
extern size_t   strlen      PP( (Const Char *) );
extern Char    *strncat     PP( (Char *, Const Char *, size_t) );
extern int      strncmp     PP( (Const Char *, Const Char *, size_t) );
extern Char    *strncpy     PP( (Char *, Const Char *, size_t) );
extern Char    *strrchr     PP( (Const Char *, int) );

extern Anyptr   memchr      PP( (Const Anyptr, int, size_t) );
extern Anyptr   memmove     PP( (Anyptr, Const Anyptr, size_t) );
extern Anyptr   memset      PP( (Anyptr, int, size_t) );
#ifndef memcpy
extern Anyptr   memcpy      PP( (Anyptr, Const Anyptr, size_t) );
extern int      memcmp      PP( (Const Anyptr, Const Anyptr, size_t) );
#endif

extern int      atoi        PP( (Const Char *) );
extern double   atof        PP( (Const Char *) );
extern long     atol        PP( (Const Char *) );
extern double   strtod      PP( (Const Char *, Char **) );
extern long     strtol      PP( (Const Char *, Char **, int) );
#endif /*P2C_H_PROTO*/

#ifndef HAS_STDLIB
extern Anyptr   malloc      PP( (size_t) );
extern Void     free        PP( (Anyptr) );
#endif

extern int      _OutMem     PV();
extern int      _CaseCheck  PV();
extern int      _NilCheck   PV();
extern int	_Escape     PP( (int) );
extern int	_EscIO      PP( (int) );

extern long     ipow        PP( (long, long) );
extern Char    *strsub      PP( (Char *, Char *, int, int) );
extern Char    *strltrim    PP( (Char *) );
extern Char    *strrtrim    PP( (Char *) );
extern Char    *strrpt      PP( (Char *, Char *, int) );
extern Char    *strpad      PP( (Char *, Char *, int, int) );
extern int      strpos2     PP( (Char *, Char *, int) );
extern long     memavail    PV();
extern int      P_peek      PP( (FILE *) );
extern int      P_eof       PP( (FILE *) );
extern int      P_eoln      PP( (FILE *) );
extern Void     P_readpaoc  PP( (FILE *, Char *, int) );
extern Void     P_readlnpaoc PP( (FILE *, Char *, int) );
extern long     P_maxpos    PP( (FILE *) );
extern Char    *P_trimname  PP( (Char *, int) );
extern long    *P_setunion  PP( (long *, long *, long *) );
extern long    *P_setint    PP( (long *, long *, long *) );
extern long    *P_setdiff   PP( (long *, long *, long *) );
extern long    *P_setxor    PP( (long *, long *, long *) );
extern int      P_inset     PP( (unsigned, long *) );
extern int      P_setequal  PP( (long *, long *) );
extern int      P_subset    PP( (long *, long *) );
extern long    *P_addset    PP( (long *, unsigned) );
extern long    *P_addsetr   PP( (long *, unsigned, unsigned) );
extern long    *P_remset    PP( (long *, unsigned) );
extern long    *P_setcpy    PP( (long *, long *) );
extern long    *P_expset    PP( (long *, long) );
extern long     P_packset   PP( (long *) );
extern int      P_getcmdline PP( (int l, int h, Char *line) );
extern Void     TimeStamp   PP( (int *Day, int *Month, int *Year,
				 int *Hour, int *Min, int *Sec) );
extern Void	P_sun_argv  PP( (char *, int, int) );


/* I/O error handling */
#define _CHKIO(cond,ior,val,def)  ((cond) ? P_ioresult=0,(val)  \
					  : P_ioresult=(ior),(def))
#define _SETIO(cond,ior)          (P_ioresult = (cond) ? 0 : (ior))

/* Following defines are suitable for the HP Pascal operating system */
#define FileNotFound     10
#define FileNotOpen      13
#define FileWriteError   38
#define BadInputFormat   14
#define EndOfFile        30

/* Creating temporary files */
#if (defined(BSD) || defined(NO_TMPFILE)) && !defined(HAVE_TMPFILE)
# define tmpfile()  (fopen(tmpnam(NULL), "w+"))
#endif

/* File buffers */
#define FILEBUF(f,sc,type) sc int __CAT__(f,_BFLAGS);   \
			   sc type __CAT__(f,_BUFFER)
#define FILEBUFNC(f,type)  int __CAT__(f,_BFLAGS);   \
			   type __CAT__(f,_BUFFER)

#define RESETBUF(f,type)   (__CAT__(f,_BFLAGS) = 1)
#define SETUPBUF(f,type)   (__CAT__(f,_BFLAGS) = 0)

#define GETFBUF(f,type)    (*((__CAT__(f,_BFLAGS) == 1 &&   \
			       ((__CAT__(f,_BFLAGS) = 2),   \
				fread(&__CAT__(f,_BUFFER),  \
				      sizeof(type),1,(f)))),\
			      &__CAT__(f,_BUFFER)))
#define AGETFBUF(f,type)   ((__CAT__(f,_BFLAGS) == 1 &&   \
			     ((__CAT__(f,_BFLAGS) = 2),   \
			      fread(__CAT__(f,_BUFFER),  \
				    sizeof(type),1,(f)))),\
			    __CAT__(f,_BUFFER))

#define PUTFBUF(f,type,v)  (GETFBUF(f,type) = (v))
#define CPUTFBUF(f,v)      (PUTFBUF(f,char,v))
#define APUTFBUF(f,type,v) (memcpy(AGETFBUF(f,type), (v),  \
				   sizeof(__CAT__(f,_BUFFER))))

#define GET(f,type)        (__CAT__(f,_BFLAGS) == 1 ?   \
			    fread(&__CAT__(f,_BUFFER),sizeof(type),1,(f)) :  \
			    (__CAT__(f,_BFLAGS) = 1))

#define PUT(f,type)        (fwrite(&__CAT__(f,_BUFFER),sizeof(type),1,(f)),  \
			    (__CAT__(f,_BFLAGS) = 0))
#define CPUT(f)            (PUT(f,char))

#define BUFEOF(f)	   (__CAT__(f,_BFLAGS) != 2 && P_eof(f))
#define BUFFPOS(f)	   (ftell(f) - (__CAT__(f,_BFLAGS) == 2))

typedef struct {
    FILE *f;
    FILEBUFNC(f,Char);
    Char name[_FNSIZE];
} _TEXT;

/* Memory allocation */
#ifdef __GCC__
# define Malloc(n)  (malloc(n) ?: (Anyptr)_OutMem())
#else
extern Anyptr __MallocTemp__;
# define Malloc(n)  ((__MallocTemp__ = malloc(n)) ? __MallocTemp__ : (Anyptr)_OutMem())
#endif
#define FreeR(p)    (free((Anyptr)(p)))    /* used if arg is an rvalue */
#define Free(p)     (free((Anyptr)(p)), (p)=NULL)

/* sign extension */
#define SEXT(x,n)   ((x) | -(((x) & (1L<<((n)-1))) << 1))

/* packed arrays */   /* BEWARE: these are untested! */
#define P_getbits_UB(a,i,n,L)   ((int)((a)[(i)>>(L)-(n)] >>   \
				       (((~(i))&((1<<(L)-(n))-1)) << (n)) &  \
				       (1<<(1<<(n)))-1))

#define P_getbits_SB(a,i,n,L)   ((int)((a)[(i)>>(L)-(n)] <<   \
				       (16 - ((((~(i))&((1<<(L)-(n))-1))+1) <<\
					      (n)) >> (16-(1<<(n))))))

#define P_putbits_UB(a,i,x,n,L) ((a)[(i)>>(L)-(n)] |=   \
				 (x) << (((~(i))&((1<<(L)-(n))-1)) << (n)))

#define P_putbits_SB(a,i,x,n,L) ((a)[(i)>>(L)-(n)] |=   \
				 ((x) & (1<<(1<<(n)))-1) <<   \
				 (((~(i))&((1<<(L)-(n))-1)) << (n)))

#define P_clrbits_B(a,i,n,L)    ((a)[(i)>>(L)-(n)] &=   \
				 ~( ((1<<(1<<(n)))-1) <<   \
				   (((~(i))&((1<<(L)-(n))-1)) << (n))) )

/* small packed arrays */
#define P_getbits_US(v,i,n)     ((int)((v) >> ((i)<<(n)) & (1<<(1<<(n)))-1))
#define P_getbits_SS(v,i,n)     ((int)((long)(v) << (SETBITS - (((i)+1) << (n))) >> (SETBITS-(1<<(n)))))
#define P_putbits_US(v,i,x,n)   ((v) |= (x) << ((i) << (n)))
#define P_putbits_SS(v,i,x,n)   ((v) |= ((x) & (1<<(1<<(n)))-1) << ((i)<<(n)))
#define P_clrbits_S(v,i,n)      ((v) &= ~( ((1<<(1<<(n)))-1) << ((i)<<(n)) ))

#define P_max(a,b)   ((a) > (b) ? (a) : (b))
#define P_min(a,b)   ((a) < (b) ? (a) : (b))


/* Fix ANSI-isms */

#ifdef LACK_LABS
# ifndef labs
#  define labs  my_labs
   extern long my_labs PP( (long) );
# endif
#endif

#ifdef LACK_MEMMOVE
# ifndef memmove
#  define memmove  my_memmove
   extern Anyptr my_memmove PP( (Anyptr, Const Anyptr, size_t) );
# endif
#endif

#ifdef LACK_MEMCPY
# ifndef memcpy
#  define memcpy  my_memcpy
   extern Anyptr my_memcpy PP( (Anyptr, Const Anyptr, size_t) );
# endif
# ifndef memcmp
#  define memcmp  my_memcmp
   extern int my_memcmp PP( (Const Anyptr, Const Anyptr, size_t) );
# endif
# ifndef memset
#  define memset  my_memset
   extern Anyptr my_memset PP( (Anyptr, int, size_t) );
# endif
#endif

/* Fix toupper/tolower on Suns and other stupid BSD systems */
#ifdef toupper
# undef toupper
# undef tolower
# define toupper(c)   my_toupper(c)
# define tolower(c)   my_tolower(c)
#endif

#ifndef _toupper
# if 'A' == 65 && 'a' == 97
#  define _toupper(c)  ((c)-'a'+'A')
#  define _tolower(c)  ((c)-'A'+'a')
# else
#  ifdef toupper
#   undef toupper   /* hope these are shadowing real functions, */
#   undef tolower   /* because my_toupper calls _toupper! */
#  endif
#  define _toupper(c)  toupper(c)
#  define _tolower(c)  tolower(c)
# endif
#endif


#endif    /* P2C_H */



/* End. */




/************************************************************
*                                                           *
*          p2clib.c                                           *
*                                                           *
************************************************************/

/* Run-time library for use with "p2c", the Pascal to C translator */

/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */



 /* #include "p2c.h" */


#ifndef NO_TIME
# include <time.h>
#endif


#define Isspace(c)  isspace(c)      /* or "((c) == ' ')" if preferred */




int P_argc;
char **P_argv;

short P_escapecode;
int P_ioresult;

long EXCP_LINE;    /* Used by Pascal workstation system */

Anyptr __MallocTemp__;

__p2c_jmp_buf *__top_jb;




void PASCAL_MAIN(argc, argv)
int argc;
char **argv;
{
    P_argc = argc;
    P_argv = argv;
    __top_jb = NULL;

#ifdef LOCAL_INIT
    LOCAL_INIT();
#endif
}





/* In case your system lacks these... */

long my_labs(x)
long x;
{
    return((x > 0) ? x : -x);
}


#ifdef __STDC__
Anyptr my_memmove(Anyptr d, Const Anyptr s, size_t n)
#else
Anyptr my_memmove(d, s, n)
Anyptr d, s;
register long n;
#endif
{
    register char *dd = (char *)d, *ss = (char *)s;
    if (dd < ss || dd - ss >= n) {
	memcpy(dd, ss, n);
    } else if (n > 0) {
	dd += n;
	ss += n;
	while (--n >= 0)
	    *--dd = *--ss;
    }
    return d;
}


#ifdef __STDC__
Anyptr my_memcpy(Anyptr d, Const Anyptr s, size_t n)
#else
Anyptr my_memcpy(d, s, n)
Anyptr d, s;
register long n;
#endif
{
    register char *ss = (char *)s, *dd = (char *)d;
    while (--n >= 0)
	*dd++ = *ss++;
    return d;
}

#ifdef __STDC__
int my_memcmp(Const Anyptr s1, Const Anyptr s2, size_t n)
#else
int my_memcmp(s1, s2, n)
Anyptr s1, s2;
register long n;
#endif
{
    register char *a = (char *)s1, *b = (char *)s2;
    register int i;
    while (--n >= 0)
	if ((i = (*a++) - (*b++)) != 0)
	    return i;
    return 0;
}

#ifdef __STDC__
Anyptr my_memset(Anyptr d, int c, size_t n)
#else
Anyptr my_memset(d, c, n)
Anyptr d;
register int c;
register long n;
#endif
{
    register char *dd = (char *)d;
    while (--n >= 0)
	*dd++ = c;
    return d;
}


int my_toupper(c)
int c;
{
    if (islower(c))
	return _toupper(c);
    else
	return c;
}


int my_tolower(c)
int c;
{
    if (isupper(c))
	return _tolower(c);
    else
	return c;
}




long ipow(a, b)
long a, b;
{
    long v;

    if (a == 0 || a == 1)
	return a;
    if (a == -1)
	return (b & 1) ? -1 : 1;
    if (b < 0)
	return 0;
    if (a == 2)
	return 1 << b;
    v = (b & 1) ? a : 1;
    while ((b >>= 1) > 0) {
	a *= a;
	if (b & 1)
	    v *= a;
    }
    return v;
}




/* Common string functions: */

/* Store in "ret" the substring of length "len" starting from "pos" (1-based).
   Store a shorter or null string if out-of-range.  Return "ret". */

char *strsub(ret, s, pos, len)
register char *ret, *s;
register int pos, len;
{
    register char *s2;

    if (--pos < 0 || len <= 0) {
        *ret = 0;
        return ret;
    }
    while (pos > 0) {
        if (!*s++) {
            *ret = 0;
            return ret;
        }
        pos--;
    }
    s2 = ret;
    while (--len >= 0) {
        if (!(*s2++ = *s++))
            return ret;
    }
    *s2 = 0;
    return ret;
}


/* Return the index of the first occurrence of "pat" as a substring of "s",
   starting at index "pos" (1-based).  Result is 1-based, 0 if not found. */

int strpos2(s, pat, pos)
char *s;
register char *pat;
register int pos;
{
    register char *cp, ch;
    register int slen;

    if (--pos < 0)
        return 0;
    slen = strlen(s) - pos;
    cp = s + pos;
    if (!(ch = *pat++))
        return 0;
    pos = strlen(pat);
    slen -= pos;
    while (--slen >= 0) {
        if (*cp++ == ch && !strncmp(cp, pat, pos))
            return cp - s;
    }
    return 0;
}


/* Case-insensitive version of strcmp. */

int strcicmp(s1, s2)
register char *s1, *s2;
{
    register unsigned char c1, c2;

    while (*s1) {
	if (*s1++ != *s2++) {
	    if (!s2[-1])
		return 1;
	    c1 = toupper(s1[-1]);
	    c2 = toupper(s2[-1]);
	    if (c1 != c2)
		return c1 - c2;
	}
    }
    if (*s2)
	return -1;
    return 0;
}




/* HP and Turbo Pascal string functions: */

/* Trim blanks at left end of string. */

char *strltrim(s)
register char *s;
{
    while (Isspace(*s++)) ;
    return s - 1;
}


/* Trim blanks at right end of string. */

char *strrtrim(s)
register char *s;
{
    register char *s2 = s;

    if (!*s)
	return s;
    while (*++s2) ;
    while (s2 > s && Isspace(*--s2))
        *s2 = 0;
    return s;
}


/* Store in "ret" "num" copies of string "s".  Return "ret". */

char *strrpt(ret, s, num)
char *ret;
register char *s;
register int num;
{
    register char *s2 = ret;
    register char *s1;

    while (--num >= 0) {
        s1 = s;
        while ((*s2++ = *s1++)) ;
        s2--;
    }
    return ret;
}


/* Store in "ret" string "s" with enough pad chars added to reach "size". */

char *strpad(ret, s, padchar, num)
char *ret;
register char *s;
register int padchar, num;
{
    register char *d = ret;

    if (s == d) {
	while (*d++) ;
    } else {
	while ((*d++ = *s++)) ;
    }
    num -= (--d - ret);
    while (--num >= 0)
	*d++ = padchar;
    *d = 0;
    return ret;
}


/* Copy the substring of length "len" from index "spos" of "s" (1-based)
   to index "dpos" of "d", lengthening "d" if necessary.  Length and
   indices must be in-range. */

void strmove(len, s, spos, d, dpos)
register char *s, *d;
register int len, spos, dpos;
{
    s += spos - 1;
    d += dpos - 1;
    while (*d && --len >= 0)
	*d++ = *s++;
    if (len > 0) {
	while (--len >= 0)
	    *d++ = *s++;
	*d = 0;
    }
}


/* Delete the substring of length "len" at index "pos" from "s".
   Delete less if out-of-range. */

void strdelete(s, pos, len)
register char *s;
register int pos, len;
{
    register int slen;

    if (--pos < 0)
        return;
    slen = strlen(s) - pos;
    if (slen <= 0)
        return;
    s += pos;
    if (slen <= len) {
        *s = 0;
        return;
    }
    while ((*s = s[len])) s++;
}


/* Insert string "src" at index "pos" of "dst". */

void strinsert(src, dst, pos)
register char *src, *dst;
register int pos;
{
    register int slen, dlen;

    if (--pos < 0)
        return;
    dlen = strlen(dst);
    dst += dlen;
    dlen -= pos;
    if (dlen <= 0) {
        strcpy(dst, src);
        return;
    }
    slen = strlen(src);
    do {
        dst[slen] = *dst;
        --dst;
    } while (--dlen >= 0);
    dst++;
    while (--slen >= 0)
        *dst++ = *src++;
}




/* File functions */

/* Peek at next character of input stream; return EOF at end-of-file. */

int P_peek(f)
FILE *f;
{
    int ch;

    ch = getc(f);
    if (ch == EOF)
	return EOF;
    ungetc(ch, f);
    return (ch == '\n') ? ' ' : ch;
}


/* Check if at end of file, using Pascal "eof" semantics.  End-of-file for
   stdin is broken; remove the special case for it to be broken in a
   different way. */

int P_eof(f)
FILE *f;
{
    register int ch;

    if (feof(f))
	return 1;
    if (f == stdin)
	return 0;    /* not safe to look-ahead on the keyboard! */
    ch = getc(f);
    if (ch == EOF)
	return 1;
    ungetc(ch, f);
    return 0;
}


/* Check if at end of line (or end of entire file). */

int P_eoln(f)
FILE *f;
{
    register int ch;

    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return (ch == '\n');
}


/* Read a packed array of characters from a file. */

Void P_readpaoc(f, s, len)
FILE *f;
char *s;
int len;
{
    int ch;

    for (;;) {
	if (len <= 0)
	    return;
	ch = getc(f);
	if (ch == EOF || ch == '\n')
	    break;
	*s++ = ch;
	--len;
    }
    while (--len >= 0)
	*s++ = ' ';
    if (ch != EOF)
	ungetc(ch, f);
}

Void P_readlnpaoc(f, s, len)
FILE *f;
char *s;
int len;
{
    int ch;

    for (;;) {
	ch = getc(f);
	if (ch == EOF || ch == '\n')
	    break;
	if (len > 0) {
	    *s++ = ch;
	    --len;
	}
    }
    while (--len >= 0)
	*s++ = ' ';
}


/* Compute maximum legal "seek" index in file (0-based). */

long P_maxpos(f)
FILE *f;
{
    long savepos = ftell(f);
    long val;

    if (fseek(f, 0L, SEEK_END))
        return -1;
    val = ftell(f);
    if (fseek(f, savepos, SEEK_SET))
        return -1;
    return val;
}


/* Use packed array of char for a file name. */

Char *P_trimname(fn, len)
register Char *fn;
register int len;
{
    static Char fnbuf[256];
    register Char *cp = fnbuf;
    
    while (--len >= 0 && *fn && !isspace(*fn))
	*cp++ = *fn++;
    *cp = 0;
    return fnbuf;
}




/* Pascal's "memavail" doesn't make much sense in Unix with virtual memory.
   We fix memory size as 10Meg as a reasonable compromise. */

long memavail()
{
    return 10000000;            /* worry about this later! */
}

long maxavail()
{
    return memavail();
}




/* Sets are stored as an array of longs.  S[0] is the size of the set;
   S[N] is the N'th 32-bit chunk of the set.  S[0] equals the maximum
   I such that S[I] is nonzero.  S[0] is zero for an empty set.  Within
   each long, bits are packed from lsb to msb.  The first bit of the
   set is the element with ordinal value 0.  (Thus, for a "set of 5..99",
   the lowest five bits of the first long are unused and always zero.) */

/* (Sets with 32 or fewer elements are normally stored as plain longs.) */

long *P_setunion(d, s1, s2)         /* d := s1 + s2 */
register long *d, *s1, *s2;
{
    long *dbase = d++;
    register int sz1 = *s1++, sz2 = *s2++;
    while (sz1 > 0 && sz2 > 0) {
        *d++ = *s1++ | *s2++;
	sz1--, sz2--;
    }
    while (--sz1 >= 0)
	*d++ = *s1++;
    while (--sz2 >= 0)
	*d++ = *s2++;
    *dbase = d - dbase - 1;
    return dbase;
}


long *P_setint(d, s1, s2)           /* d := s1 * s2 */
register long *d, *s1, *s2;
{
    long *dbase = d++;
    register int sz1 = *s1++, sz2 = *s2++;
    while (--sz1 >= 0 && --sz2 >= 0)
        *d++ = *s1++ & *s2++;
    while (--d > dbase && !*d) ;
    *dbase = d - dbase;
    return dbase;
}


long *P_setdiff(d, s1, s2)          /* d := s1 - s2 */
register long *d, *s1, *s2;
{
    long *dbase = d++;
    register int sz1 = *s1++, sz2 = *s2++;
    while (--sz1 >= 0 && --sz2 >= 0)
        *d++ = *s1++ & ~*s2++;
    if (sz1 >= 0) {
        while (sz1-- >= 0)
            *d++ = *s1++;
    }
    while (--d > dbase && !*d) ;
    *dbase = d - dbase;
    return dbase;
}


long *P_setxor(d, s1, s2)         /* d := s1 / s2 */
register long *d, *s1, *s2;
{
    long *dbase = d++;
    register int sz1 = *s1++, sz2 = *s2++;
    while (sz1 > 0 && sz2 > 0) {
        *d++ = *s1++ ^ *s2++;
	sz1--, sz2--;
    }
    while (--sz1 >= 0)
	*d++ = *s1++;
    while (--sz2 >= 0)
	*d++ = *s2++;
    while (--d > dbase && !*d) ;
    *dbase = d - dbase;
    return dbase;
}


int P_inset(val, s)                 /* val IN s */
register unsigned val;
register long *s;
{
    register int bit;
    bit = val % SETBITS;
    val /= SETBITS;
    if (val < *s++ && ((1<<bit) & s[val]))
	return 1;
    return 0;
}


long *P_addset(s, val)              /* s := s + [val] */
register long *s;
register unsigned val;
{
    register long *sbase = s;
    register int bit, size;
    bit = val % SETBITS;
    val /= SETBITS;
    size = *s;
    if (++val > size) {
        s += size;
        while (val > size)
            *++s = 0, size++;
        *sbase = size;
    } else
        s += val;
    *s |= 1<<bit;
    return sbase;
}


long *P_addsetr(s, v1, v2)              /* s := s + [v1..v2] */
register long *s;
register unsigned v1, v2;
{
    register long *sbase = s;
    register int b1, b2, size;
    if ((int)v1 > (int)v2)
	return sbase;
    b1 = v1 % SETBITS;
    v1 /= SETBITS;
    b2 = v2 % SETBITS;
    v2 /= SETBITS;
    size = *s;
    v1++;
    if (++v2 > size) {
        while (v2 > size)
            s[++size] = 0;
        s[v2] = 0;
        *s = v2;
    }
    s += v1;
    if (v1 == v2) {
        *s |= (~((-2)<<(b2-b1))) << b1;
    } else {
        *s++ |= (-1) << b1;
        while (++v1 < v2)
            *s++ = -1;
        *s |= ~((-2) << b2);
    }
    return sbase;
}


long *P_remset(s, val)              /* s := s - [val] */
register long *s;
register unsigned val;
{
    register int bit;
    bit = val % SETBITS;
    val /= SETBITS;
    if (++val <= *s) {
	if (!(s[val] &= ~(1<<bit)))
	    while (*s && !s[*s])
		(*s)--;
    }
    return s;
}


int P_setequal(s1, s2)              /* s1 = s2 */
register long *s1, *s2;
{
    register int size = *s1++;
    if (*s2++ != size)
        return 0;
    while (--size >= 0) {
        if (*s1++ != *s2++)
            return 0;
    }
    return 1;
}


int P_subset(s1, s2)                /* s1 <= s2 */
register long *s1, *s2;
{
    register int sz1 = *s1++, sz2 = *s2++;
    if (sz1 > sz2)
        return 0;
    while (--sz1 >= 0) {
        if (*s1++ & ~*s2++)
            return 0;
    }
    return 1;
}


long *P_setcpy(d, s)                /* d := s */
register long *d, *s;
{
    register long *save_d = d;

#ifdef SETCPY_MEMCPY
    memcpy(d, s, (*s + 1) * sizeof(long));
#else
    register int i = *s + 1;
    while (--i >= 0)
        *d++ = *s++;
#endif
    return save_d;
}


/* s is a "smallset", i.e., a 32-bit or less set stored
   directly in a long. */

long *P_expset(d, s)                /* d := s */
register long *d;
register long s;
{
    if (s) {
	d[1] = s;
	*d = 1;
    } else
        *d = 0;
    return d;
}


long P_packset(s)                   /* convert s to a small-set */
register long *s;
{
    if (*s++)
        return *s;
    else
        return 0;
}





/* Oregon Software Pascal extensions, courtesy of William Bader */

int P_getcmdline(l, h, line)
int l, h;
Char *line;
{
    int i, len;
    char *s;
    
    h = h - l + 1;
    len = 0;
    for(i = 1; i < P_argc; i++) {
	s = P_argv[i];
	while (*s) {
	    if (len >= h) return len;
	    line[len++] = *s++;
	}
	if (len >= h) return len;
	line[len++] = ' ';
    }
    return len;
}

Void TimeStamp(Day, Month, Year, Hour, Min, Sec)
int *Day, *Month, *Year, *Hour, *Min, *Sec;
{
#ifndef NO_TIME
    struct tm *tm;
    long clock;

    time(&clock);
    tm = localtime(&clock);
    *Day = tm->tm_mday;
    *Month = tm->tm_mon + 1;		/* Jan = 0 */
    *Year = tm->tm_year;
    if (*Year < 1900)
	*Year += 1900;     /* year since 1900 */
    *Hour = tm->tm_hour;
    *Min = tm->tm_min;
    *Sec = tm->tm_sec;
#endif
}




/* SUN Berkeley Pascal extensions */

Void P_sun_argv(s, len, n)
register char *s;
register int len, n;
{
    register char *cp;

    if ((unsigned)n < P_argc)
	cp = P_argv[n];
    else
	cp = "";
    while (*cp && --len >= 0)
	*s++ = *cp++;
    while (--len >= 0)
	*s++ = ' ';
}




int _OutMem()
{
    return _Escape(-2);
}

int _CaseCheck()
{
    return _Escape(-9);
}

int _NilCheck()
{
    return _Escape(-3);
}





/* The following is suitable for the HP Pascal operating system.
   It might want to be revised when emulating another system. */

char *_ShowEscape(buf, code, ior, prefix)
char *buf, *prefix;
int code, ior;
{
    char *bufp;

    if (prefix && *prefix) {
        strcpy(buf, prefix);
	strcat(buf, ": ");
        bufp = buf + strlen(buf);
    } else {
        bufp = buf;
    }
    if (code == -10) {
        sprintf(bufp, "Pascal system I/O error %d", ior);
        switch (ior) {
            case 3:
                strcat(buf, " (illegal I/O request)");
                break;
            case 7:
                strcat(buf, " (bad file name)");
                break;
            case FileNotFound:   /*10*/
                strcat(buf, " (file not found)");
                break;
            case FileNotOpen:    /*13*/
                strcat(buf, " (file not open)");
                break;
            case BadInputFormat: /*14*/
                strcat(buf, " (bad input format)");
                break;
            case 24:
                strcat(buf, " (not open for reading)");
                break;
            case 25:
                strcat(buf, " (not open for writing)");
                break;
            case 26:
                strcat(buf, " (not open for direct access)");
                break;
            case 28:
                strcat(buf, " (string subscript out of range)");
                break;
            case EndOfFile:      /*30*/
                strcat(buf, " (end-of-file)");
                break;
            case FileWriteError: /*38*/
		strcat(buf, " (file write error)");
		break;
        }
    } else {
        sprintf(bufp, "Pascal system error %d", code);
        switch (code) {
            case -2:
                strcat(buf, " (out of memory)");
                break;
            case -3:
                strcat(buf, " (reference to NIL pointer)");
                break;
            case -4:
                strcat(buf, " (integer overflow)");
                break;
            case -5:
                strcat(buf, " (divide by zero)");
                break;
            case -6:
                strcat(buf, " (real math overflow)");
                break;
            case -8:
                strcat(buf, " (value range error)");
                break;
            case -9:
                strcat(buf, " (CASE value range error)");
                break;
            case -12:
                strcat(buf, " (bus error)");
                break;
            case -20:
                strcat(buf, " (stopped by user)");
                break;
        }
    }
    return buf;
}


int _Escape(code)
int code;
{
    char buf[100];

    P_escapecode = code;
    if (__top_jb) {
	__p2c_jmp_buf *jb = __top_jb;
	__top_jb = jb->next;
	longjmp(jb->jbuf, 1);
    }
    if (code == 0)
        exit(EXIT_SUCCESS);
    if (code == -1)
        exit(EXIT_FAILURE);
    fprintf(stderr, "%s\n", _ShowEscape(buf, P_escapecode, P_ioresult, ""));
    exit(EXIT_FAILURE);
}

int _EscIO(code)
int code;
{
    P_ioresult = code;
    return _Escape(-10);
}




/* End. */





/************************************************************
*                                                           *
*          date.c                                           *
*                                                           *
************************************************************/



/*#include <time.h>
#include <stdio.h>*/
static  char *mon[]={"JAN","FEB","MAR","APR","MAY","JUN",
                "JUL","AUG","SEP","OCT","NOV","DEC"};
/* PROCEDURE DATE(VAR DATESTRING:PACKED ARRAY[1..11] OF CHAR);EXTERN; */
/* activate DATE by removing comment brackets if necessary */
/***/

void Date(string)
char* string;
{
  time_t tt;
  struct tm *t;
  time(&tt);
  t=localtime(&tt);
  sprintf(string,"%d-%s-19%d\n",t->tm_mday,mon[t->tm_mon],t->tm_year);
}


/* p2c: dssp.p, line 295: 
 * Note: Unexpected name "tapein" in program header [262] */
/* p2c: dssp.p, line 295: 
 * Note: Unexpected name "tapeout" in program header [262] */


/*--------------------------------------------------------------------*/
/* PROGRAM FATAL ERROR EXIT LABEL */
/*******************  MATHEMATICAL CONSTANTS  **************************
 YVERTEX, - ARE Y,Z-COMPONENTS OF THE FIRST ICOSAHEDRON VERTEX. THE
 ZVERTEX    X-COMPONENT IS 0.
 EPS      - NUMERICAL TOLERANCE
  --------------------------------------------------------------------*/

#define PIHALF          1.570796
#define PI              3.141593
#define TWOPI           6.283185
#define FOURPI          12.56637
#define RADIAN          57.29578
#define YVERTEX         0.8506508
#define ZVERTEX         0.5257311
#define EPS             0.00001
/***/
/***************  ARRAY DIMENSIONING CONSTANTS  ***********************
 NMAX     - MAXIMUM NUMBER OF AMINOACID RESIDUES IN ARRAY CHAIN
 MAXATOM  - MAXIMUM NUMBER OF SIDECHAIN ATOMS IN ARRAY SIDECHAIN
 MAXBRIDGE- MAXIMUM NUMBER OF BRIDGES IN ARRAY BRIDGETABLE
 NFACE,   - NUMBER OF FACES OF POLYHEDRON. THE COORDINATES OF THE CENTRE
 ORDER      OF EACH TRIANGULAR FACE ARE STORED IN ARRAY P, THE AREA
             IS STORED IN ARRAY WP IN PROCEDURE FLAGACCESS. NFACE MUST
             BE OF THE FORM NFACE=20*(4**ORDER), ORDER=0,1,2,...
             THE ACCURACY OF THE SOLVENT ACCESSIBLE SURFACE OF EACH
             AMINOACID RESIDUE IS ONE ANGSTROM**2 FOR ORDER=2,NFACE=320.
 MAXPACK  - MAXIMUM NUMBER OF PROTEIN ATOMS WHICH CAN INTRUDE INTO
             SOLVENT AROUND ANY GIVEN TEST ATOM. THE COORDINATES OF
             THESE ATOMS ARE STORED IN ARRAY X, THEIR RADII IN ARRAY RX
             IN PROCEDURE SURFACE.
 MAXHIST  - NUMBER OF SLOTS IN ARRAYS HELIXHIST AND BETAHIST USED FOR
             LENGTH STATISTICS OF SECONDARY STRUCTURE.
 MAXSS    - MAXIMUM NUMBER OF SSBOND RECORDS ON INPUT FILE. THE
             DISULFIDE BOND ARE SAVED IN ARRAY SSBONDS.
  --------------------------------------------------------------------*/

#define NMAX            6000

#define MAXATOM         40000L

#define MAXBRIDGE       300
#define NFACE           320
#define ORDER           2
#define MAXPACK         200
#define MAXHIST         30
#define MAXSS           100
/***/
/*********************  PHYSICAL CONSTANTS   **************************
 RN       - RADIUS OF PEPTIDE NITROGEN ATOM
 RCA      - RADIUS OF PEPTIDE ALPHA-CARBON ATOM
 RC       - RADIUS OF PEPTIDE C'-CARBON ATOM
 RO       - RADIUS OF PEPTIDE OXYGEN ATOM
 RSIDEATOM- RADIUS OF SIDECHAIN ATOM
 RWATER   - RADIUS OF WATER MOLECULE
 SSDIST   - MAXIMUM ALLOWED DISTANCE OF DISULFIDE BRIDGE
 BREAKDIST- MAXIMUM ALLOWED PEPTIDE BOND LENGTH. IF DISTANCE IS
             GREATER A POLYPEPTIDE CHAIN INTERRUPTION IS ASSUMED.
 RESRAD   - MAXIMUM RADIUS OF A SPHERE AROUND C-ALPHA CONTAINING
             ALL ATOMS OF A RESIDUE
 CADIST   - MINIMUM DISTANCE BETWEEN ALPHA-CARBON ATOMS SUCH THAT NO
             BACKBONE HYDROGEN BONDS CAN BE FORMED
 DIST     - SMALLEST ALLOWED DISTANCE BETWEEN ANY ATOMS
 MAXDIST  - LARGEST ALLOWED DISTANCE BETWEEN SIDECHAIN ATOM AND C-ALPHA
             WITHIN A RESIDUE
 Q        - COUPLING CONSTANT FOR ELECTROSTATIC ENERGY
                    Q=-332*0.42*0.2*1000.0
 HBLOW    - LOWEST ALLOWED  ENERGY OF A HYDROGEN BOND IN CAL/MOL
 HBHIGH   - HIGHEST ALLOWED ENERGY OF A HYDROGEN BOND IN CAL/MOL
  --------------------------------------------------------------------*/

#define RN              1.65
#define RCA             1.87
#define RC              1.76
#define RO              1.4
#define RSIDEATOM       1.8
#define RWATER          1.4
#define SSDIST          3.0
#define BREAKDIST       2.5
#define RESRAD          10.0
#define CADIST          9.0
#define DIST            0.5
#define MAXDIST         10.0
#define Q               (-27888.0)

#define HBLOW           (-9900)
#define HBHIGH          (-500)


/***/
/***************** GLOBAL DATA TYPE DEFINITIONS ***********************/

typedef double vector[3];
typedef Char char4[4];
typedef Char char6[6];
typedef enum {
  parallel, antiparallel, nobridge
} bridgetyp;
typedef enum {
  symbol, turn3, turn4, turn5, bend, chirality, beta1, beta2
} structure;

typedef struct hydrogenbond {
  long residue, energy;
} hydrogenbond;

typedef hydrogenbond bonds[2];

typedef struct backbone {
  char6 aaident;
  Char sheetlabel, aa;
  char4 threelettercode;
  Char ss[(long)beta2 - (long)symbol + 1];
  long partner[(long)beta2 - (long)beta1 + 1];
  long access;
  double alpha, kappa;
  bonds acceptor, donor;
  vector h, n, ca, c, o;
  long atompointer, nsideatoms;
} backbone;

typedef struct bridge {
  Char sheetname, laddername;
  bridgetyp btyp;
  long linkset[MAXBRIDGE / 32 + 2];
  long ib, ie, jb, je, from, towards;
} bridge;


Static int noaccFlag,silentFlag;
Static long nss, nssintra, nssinter, lchain, nbridge;
Static char6 ssbonds[MAXSS][2];
Static backbone chain[NMAX + 1];
Static FILE *tapein, *tapeout;
Static vector sidechain[MAXATOM];
Static bridge bridgetable[MAXBRIDGE];

Static Void VecCopy(dest,source)
double* dest;
double* source;
{
  dest[0]=source[0];
  dest[1]=source[1];
  dest[2]=source[2];
}

Static Void StrCopy(dest,source,n)
char* dest;
char* source;
int n;
{
  int i;
  for(i=0;i<n;i++)
    dest[i]=source[i];
}

Static double Atan2(y, x)
double y, x;
{
  double z;

  if (x != 0.0)
    z = atan(y / x);
  else if (y > 0.0)
    z = PIHALF;
  else if (y < 0.0)
    z = -PIHALF;
  else
    z = TWOPI;
  if (x >= 0.0)
    return z;
  if (y > 0.0)
    z += PI;
  else
    z -= PI;
  return z;
}  /* Atan2 */


/***/

Static Void Diff(x, y, z)
double *x, *y, *z;
{
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  z[2] = x[2] - y[2];
}  /* Diff */


/***/

Static double Dot(x, y)
double *x, *y;
{
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}  /* Dot */


/***/

Static Void Cross(x, y, z)
double *x, *y, *z;
{
  z[0] = x[1] * y[2] - y[1] * x[2];
  z[1] = x[2] * y[0] - y[2] * x[0];
  z[2] = x[0] * y[1] - y[0] * x[1];
}  /* Cross */


/***/

Static Void Norm(x, xnorm)
double *x;
double *xnorm;
{
  /* RETURNS INPUT VECTOR X NORMALIZED TO UNIT LENGTH.
     XNORM IS THE ORIGINAL LENGTH OF X.                         */
  double TEMP, TEMP1, TEMP2;

  TEMP = x[0];
  TEMP1 = x[1];
  TEMP2 = x[2];
  *xnorm = TEMP * TEMP + TEMP1 * TEMP1 + TEMP2 * TEMP2;
  if (*xnorm <= 0.0)
    return;
  *xnorm = sqrt(*xnorm);
  x[0] /= *xnorm;
  x[1] /= *xnorm;
  x[2] /= *xnorm;
}  /* Norm */


/***/

Static double Dihedralangle(v1, v2, v3, v4)
double *v1, *v2, *v3, *v4;
{
  /*CALCULATES TORSION ANGLE OF A SET OF 4 ATOMS V1-V2-V3-V4.
    DIHEDRALANGLE IS THE ANGLE BETWEEN THE PROJECTION OF
    V1-V2 AND THE PROJECTION OF V4-V3 ONTO A PLANE NORMAL TO
    BOND V2-V3.*/
  /***/
  double Result, u, v;
  vector v12, v43, x, y, z, p;

  Diff(v1, v2, v12);
  Diff(v4, v3, v43);
  Diff(v2, v3, z);
  Cross(z, v12, p);
  Cross(z, v43, x);
  Cross(z, x, y);
  u = Dot(x, x);
  v = Dot(y, y);
  Result = 360.0;
  if (u <= 0.0 || v <= 0.0)
    return Result;
  u = Dot(p, x) / sqrt(u);
  v = Dot(p, y) / sqrt(v);
  if (u != 0.0 || v != 0.0)
    return (Atan2(v, u) * RADIAN);
  return Result;
}  /* Dihedralangle */


/***/

Static double Cosangle(v1, v2, v3, v4)
double *v1, *v2, *v3, *v4;
{
  vector u, v;
  double x;

  Diff(v1, v2, u);
  Diff(v3, v4, v);
  x = Dot(u, u) * Dot(v, v);
  if (x > 0.0)
    return (Dot(u, v) / sqrt(x));
  else
    return 0.0;
}  /* Cosangle */


/***/

Static double Distance(u, v)
double *u, *v;
{
  double TEMP, TEMP1, TEMP2;

  TEMP = u[0] - v[0];
  TEMP1 = u[1] - v[1];
  TEMP2 = u[2] - v[2];
  return sqrt(TEMP * TEMP + TEMP1 * TEMP1 + TEMP2 * TEMP2);
}  /* Distance */


/***/

Static double Distsq(u, v)
double *u, *v;
{
  double TEMP, TEMP1, TEMP2;

  TEMP = u[0] - v[0];
  TEMP1 = u[1] - v[1];
  TEMP2 = u[2] - v[2];
  return (TEMP * TEMP + TEMP1 * TEMP1 + TEMP2 * TEMP2);
}  /* Distsq */


/*--------------------------------------------------------------------*/

Static boolean Nochainbreak(i, j)
long i, j;
{
  long k;
  boolean test;

  test = (i >= 1 && j <= NMAX && i <= j);
  k = i;
  while (test && k <= j) {
    if (chain[k].aa == '!')
      test = false;
    else
      k++;
  }
  return test;
}  /* Nochainbreak */


/***/
/*--------------------------------------------------------------------*/

Static Void Writeresidue(res)
backbone res;
{
  long i;

  for (i = 0; i <= 3; i++)
    putchar(res.threelettercode[i]);
  for (i = 0; i <= 5; i++)
    putchar(res.aaident[i]);
}  /* Writeresidue */


#define MAXSIDEATOMS    20


typedef enum {
  headercard, compndcard, sourcecard, authorcard, ssbondcard, atomcard,
  tercard, endcard, othercard
} cardtype;
/***/

typedef struct cardcontents {
  cardtype art;
  union {
    Char z[128];
    char6 r[2];
    struct {
      char4 atomname, aaname;
      Char altloc, residuename;
      char6 reseqnum;
      vector coordinates;
    } U5;
    Char ch;
  } UU;
} cardcontents;   /* CARDCONTENTS TYPE DEFINITION */

/***/


Static jmp_buf _JL99;

/* Local variables for Inputcoordinates: */
struct LOC_Inputcoordinates {
  long *lchain, latom, hatoms;
  boolean nmissing, camissing, cmissing, omissing, corelimit;
  vector sidecoordinates[MAXSIDEATOMS];
  double dco;
  char4 sideatomnames[MAXSIDEATOMS];
  backbone reszero, resinfo;
} ;

/***/

Local Char Onelettercode(aaa, LINK)
Char *aaa;
struct LOC_Inputcoordinates *LINK;
{
  Char aasymbol[50];
  Char aminoacid[150];
  Char string[5][30];
  long i, l, k;
  Char a;

  StrCopy(aasymbol, "ARNDCEQGHILKMFPSTWYVBZXXXXXXXXXXXXXXXX--CCCCIPPPW-", 50L);
  StrCopy(string[0], "ALAARGASNASPCYSGLUGLNGLYHISILE", 30L);
  StrCopy(string[1], "LEULYSMETPHEPROSERTHRTRPTYRVAL", 30L);
  StrCopy(string[2], "ASXGLXACDALBALIABUAROBASBETHSE", 30L);
  StrCopy(string[3], "HYPHYLORNPCASARTAUTHYUNKACEFOR", 30L);
  StrCopy(string[4], "CYHCSHCSSCYXILUPRZPR0CPRTRYHOH", 30L);
  l = 0;
  for (k = 0; k <= 4; k++) {
    for (i = 0; i <= 29; i++) {
      l++;
      aminoacid[l - 1] = string[k][i];
    }
  }
  a = '-';
  i = 1;
  k = 1;
  while (k < 51 && a == '-') {
    if (aminoacid[i - 1] == aaa[0]) {
      if (aminoacid[i] == aaa[1]) {
	if (aminoacid[i + 1] == aaa[2])
	  a = aasymbol[k - 1];
      }
    }
    i += 3;
    k++;
  }
  return a;
}  /* Onelettercode */

/* Local variables for Checksideatoms: */
struct LOC_Checksideatoms {
  struct LOC_Inputcoordinates *LINK;
} ;

/***/

Local Void Checkdist(resinfo, LINK)
backbone *resinfo;
struct LOC_Checksideatoms *LINK;
{
  long i, j, FORLIM;

  i = 1;
  while (i <= resinfo->nsideatoms) {
    if (Distance(resinfo->ca, LINK->LINK->sidecoordinates[i - 1]) <= MAXDIST) {
      i++;
      continue;
    }
    printf(" !!! RESIDUE ");
    Writeresidue(*resinfo);
    printf(" HAS ILLEGAL SIDECHAIN ATOM NAMED ");
    for (j = 0; j <= 3; j++)
      putchar(LINK->LINK->sideatomnames[i - 1][j]);
    printf(".\n");
    printf("     THIS ATOM WILL BE IGNORED !!!\n\n");
    FORLIM = resinfo->nsideatoms;
    for (j = i + 1; j <= FORLIM; j++) {
      StrCopy(LINK->LINK->sideatomnames[j - 2],
	     LINK->LINK->sideatomnames[j - 1], sizeof(char4));
      VecCopy(LINK->LINK->sidecoordinates[j - 2],
	     LINK->LINK->sidecoordinates[j - 1]);
    }
    resinfo->nsideatoms--;
  }
}  /* Checkdist */

/***/

Local Void Checksideatoms(resinfo, LINK)
backbone *resinfo;
struct LOC_Inputcoordinates *LINK;
{
  struct LOC_Checksideatoms V;
  long i, j;
  Char c;

  /***/

  V.LINK = LINK;
  Checkdist(resinfo, &V);
  i = -1;
  c = resinfo->aa;
  if (c == 'G')
    i = 0;
  if (c == 'A')
    i = 1;
  if (c == 'S' || c == 'C')
    i = 2;
  if (c == 'P' || c == 'T' || c == 'V')
    i = 3;
  if (c == 'B' || c == 'M' || c == 'L' || c == 'I' || c == 'D' || c == 'N')
    i = 4;
  if (c == 'Z' || c == 'K' || c == 'Q' || c == 'E')
    i = 5;
  if (c == 'H')
    i = 6;
  if (c == 'F' || c == 'R')
    i = 7;
  if (c == 'Y')
    i = 8;
  if (c == 'W')
    i = 10;
  if (resinfo->nsideatoms < i) {
    printf(" !!! RESIDUE ");
    Writeresidue(*resinfo);
    printf(" HAS%3ld INSTEAD OF EXPECTED ", resinfo->nsideatoms);
    printf("%3ld SIDECHAIN ATOMS.\n", i);
    printf(
      "     CALCULATED SOLVENT ACCESSIBILITY REFERS TO INCOMPLETE SIDECHAIN !!!\n\n");
  }
  if (i == -1 || resinfo->nsideatoms <= i)
    return;
  printf(" !!! RESIDUE ");
  Writeresidue(*resinfo);
  printf(" HAS%3ld INSTEAD OF EXPECTED ", resinfo->nsideatoms);
  printf("%3ld SIDECHAIN ATOMS.\n", i);
  printf("     LAST SIDECHAIN ATOM NAME IS ");
  for (j = 0; j <= 3; j++)
    putchar(LINK->sideatomnames[resinfo->nsideatoms - 1][j]);
  printf("\n     CALCULATED SOLVENT ACCESSIBILITY INCLUDES EXTRA ATOMS !!!\n\n");
}  /* Checksideatoms */

/***/

Local Void Putresidue(LINK)
struct LOC_Inputcoordinates *LINK;
{
  /* insert residue into protein chain */
  long i;
  boolean complete;
  long FORLIM;

  complete = !(LINK->nmissing || LINK->camissing || LINK->cmissing ||
	       LINK->omissing);
  if (!complete &&
      strncmp(LINK->reszero.aaident, LINK->resinfo.aaident, sizeof(char6))) {
    printf(" !!! BACKBONE INCOMPLETE FOR RESIDUE ");
    Writeresidue(LINK->resinfo);
    printf("\n     RESIDUE WILL BE IGNORED !!!\n\n");
  }
  LINK->corelimit = (LINK->latom + LINK->resinfo.nsideatoms > MAXATOM ||
		     *LINK->lchain > NMAX - 2);
  if (complete && !LINK->corelimit) {
    Checksideatoms(&LINK->resinfo, LINK);
    VecCopy(LINK->resinfo.h, LINK->resinfo.n);
    if (Nochainbreak(*LINK->lchain, *LINK->lchain)) {
      if (Distance(chain[*LINK->lchain].c, LINK->resinfo.n) > BREAKDIST)
	  /* keep ! at LCHAIN */
	  {  /* CS Oct 1987 */
	printf(" !!! EXCESSIVE C TO N DISTANCE ");
	printf("% .5E>% .5E\n",
	       Distance(chain[*LINK->lchain].c, LINK->resinfo.n), BREAKDIST);
	printf("     BEFORE RESIDUE ");
	Writeresidue(LINK->resinfo);
	printf(". CHAIN BREAK RESIDUE INSERTED !!!\n\n");
	(*LINK->lchain)++;
      }
    }
    if (Nochainbreak(*LINK->lchain, *LINK->lchain) && LINK->resinfo.aa != 'P') {
      LINK->dco = Distance(chain[*LINK->lchain].c, chain[*LINK->lchain].o);
      for (i = 0; i <= 2; i++)
	LINK->resinfo.h[i] = LINK->resinfo.n[i] +
	    (chain[*LINK->lchain].c[i] - chain[*LINK->lchain].o[i]) / LINK->dco;
    }
    (*LINK->lchain)++;
    chain[*LINK->lchain] = LINK->resinfo;
    FORLIM = LINK->resinfo.nsideatoms;
    for (i = 0; i < FORLIM; i++)
      VecCopy(sidechain[LINK->latom + i], LINK->sidecoordinates[i]);
    LINK->latom += LINK->resinfo.nsideatoms;
  }
  if (Nochainbreak(*LINK->lchain, *LINK->lchain) && !complete)
    (*LINK->lchain)++;
  LINK->resinfo = LINK->reszero;
  LINK->nmissing = true;
  LINK->camissing = true;
  LINK->cmissing = true;
  LINK->omissing = true;
}  /* Putresidue */

/***/

Local Void Getresidue(atomname, coordinates, LINK)
Char *atomname;
double *coordinates;
struct LOC_Inputcoordinates *LINK;
{
  boolean hydrogenatom;

  hydrogenatom = ((atomname[0] == '9' || atomname[0] == '8' ||
		   atomname[0] == '7' || atomname[0] == '6' ||
		   atomname[0] == '5' || atomname[0] == '4' ||
		   atomname[0] == '3' || atomname[0] == '2' ||
		   atomname[0] == '1' || atomname[0] == '0' ||
		   atomname[0] == ' ') &&
		  (atomname[1] == 'D' || atomname[1] == 'H'));
  if (hydrogenatom) {
    LINK->hatoms++;
    return;
  }
  if (!strncmp(atomname, " N  ", sizeof(char4))) {
    LINK->nmissing = false;
    VecCopy(LINK->resinfo.n, coordinates);
    return;
  }
  if (!strncmp(atomname, " CA ", sizeof(char4))) {
    LINK->camissing = false;
    VecCopy(LINK->resinfo.ca, coordinates);
    return;
  }
  if (!strncmp(atomname, " C  ", sizeof(char4))) {
    LINK->cmissing = false;
    VecCopy(LINK->resinfo.c, coordinates);
    return;
  }
  if (!strncmp(atomname, " O  ", sizeof(char4))) {
    LINK->omissing = false;
    VecCopy(LINK->resinfo.o, coordinates);
    return;
  }
  if (LINK->resinfo.nsideatoms >= MAXSIDEATOMS)
    return;
  LINK->resinfo.nsideatoms++;
  VecCopy(LINK->sidecoordinates[LINK->resinfo.nsideatoms - 1], coordinates);
  StrCopy(LINK->sideatomnames[LINK->resinfo.nsideatoms - 1], atomname,
	 sizeof(char4));
}  /* Getresidue */

/***/

Local Void Readcard(cardinfo, LINK)
cardcontents *cardinfo;
struct LOC_Inputcoordinates *LINK;
{
  Char c;
  long k, l, m;
  char6 key;

  cardinfo->art = othercard;
  do {
    if (!P_eof(tapein)) {
      *key = getc(tapein);
      if (key[0] == '\n')
	key[0] = ' ';
    }
  } while (!(isupper(key[0]) | P_eof(tapein)));
  if (P_eof(tapein)) {
    cardinfo->art = endcard;
    return;
  }
  for (l = 1; l <= 5; l++) {
    if (!P_eoln(tapein)) {
      key[l] = getc(tapein);
      if (key[l] == '\n')
	key[l] = ' ';
    }
  }
  if (!strncmp(key, "HEADER", sizeof(char6)))
    cardinfo->art = headercard;
  if (!strncmp(key, "COMPND", sizeof(char6)))
    cardinfo->art = compndcard;
  if (!strncmp(key, "SOURCE", sizeof(char6)))
    cardinfo->art = sourcecard;
  if (!strncmp(key, "AUTHOR", sizeof(char6)))
    cardinfo->art = authorcard;
  if (!strncmp(key, "SSBOND", sizeof(char6)))
    cardinfo->art = ssbondcard;
  if (!strncmp(key, "ATOM  ", sizeof(char6)))
    cardinfo->art = atomcard;
  if (!strncmp(key, "TER   ", sizeof(char6)))
    cardinfo->art = tercard;
  if (!strncmp(key, "END   ", sizeof(char6)))
    cardinfo->art = endcard;
  switch (cardinfo->art) {

  case headercard:
  case compndcard:
  case sourcecard:
  case authorcard:
    for (l = 0; l <= 5; l++)
      cardinfo->UU.z[l] = key[l];
    for (l = 6; l <= 126; l++)
      cardinfo->UU.z[l] = ' ';
    cardinfo->UU.z[127] = '.';
    if (cardinfo->art == headercard)
      m = 66;
    else
      m = 70;
    for (l = 6; l < m; l++) {
      if (!P_eoln(tapein)) {
	cardinfo->UU.z[l] = getc(tapein);
	if (cardinfo->UU.z[l] == '\n')
	  cardinfo->UU.z[l] = ' ';
      }
    }
    break;

  case ssbondcard:
    for (l = 7; l <= 8; l++) {
      c = getc(tapein);
      if (c == '\n')
	c = ' ';
    }
    for (k = 0; k <= 1; k++) {
      for (l = 1; l <= 7; l++) {
	c = getc(tapein);
	if (c == '\n')
	  c = ' ';
      }
      cardinfo->UU.r[k][5] = getc(tapein);
      c = getc(tapein);
      if (cardinfo->UU.r[k][5] == '\n')
	cardinfo->UU.r[k][5] = ' ';
      if (c == '\n')
	c = ' ';
      /* minor modification suggested by Steven Sheriff */
      for (l = 0; l <= 3; l++) {
	cardinfo->UU.r[k][l] = getc(tapein);
	if (cardinfo->UU.r[k][l] == '\n')
	  cardinfo->UU.r[k][l] = ' ';
      }
      if (P_eoln(tapein))
	cardinfo->UU.r[k][4] = ' ';
      else {
	cardinfo->UU.r[k][4] = getc(tapein);
	if (cardinfo->UU.r[k][4] == '\n')
	  cardinfo->UU.r[k][4] = ' ';
      }
    }
    /* end minor modification suggested by Steven Sheriff */
    break;

  case atomcard:
    for (l = 7; l <= 12; l++) {
      c = getc(tapein);
      if (c == '\n')
	c = ' ';
    }
    for (l = 0; l <= 3; l++) {
      cardinfo->UU.U5.atomname[l] = getc(tapein);
      if (cardinfo->UU.U5.atomname[l] == '\n')
	cardinfo->UU.U5.atomname[l] = ' ';
    }
    cardinfo->UU.U5.altloc = getc(tapein);
    if (cardinfo->UU.U5.altloc == '\n')
      cardinfo->UU.U5.altloc = ' ';
    for (l = 0; l <= 2; l++) {
      cardinfo->UU.U5.aaname[l] = getc(tapein);
      if (cardinfo->UU.U5.aaname[l] == '\n')
	cardinfo->UU.U5.aaname[l] = ' ';
    }
    cardinfo->UU.U5.aaname[3] = ' ';
    cardinfo->UU.U5.residuename = Onelettercode(cardinfo->UU.U5.aaname, LINK);
    c = getc(tapein);
    cardinfo->UU.U5.reseqnum[5] = getc(tapein);
    if (c == '\n')
      c = ' ';
    if (cardinfo->UU.U5.reseqnum[5] == '\n')
      cardinfo->UU.U5.reseqnum[5] = ' ';
    for (l = 0; l <= 4; l++) {
      cardinfo->UU.U5.reseqnum[l] = getc(tapein);
      if (cardinfo->UU.U5.reseqnum[l] == '\n')
	cardinfo->UU.U5.reseqnum[l] = ' ';
    }
    for (l = 0; l <= 2; l++)
      fscanf(tapein, "%lf", &cardinfo->UU.U5.coordinates[l]);
    break;

  case tercard:
  case endcard:
  case othercard:
    /* blank case */
    break;
  }
  fscanf(tapein, "%*[^\n]");
  getc(tapein);
}  /* Readcard */


/***/
/*--------------------------------------------------------------------*/
/* SEE BROOKHAVEN PROTEIN DATA BANK ATOMIC COORDINATE ENTRY FORMAT
                                    OF DEC. 1981.
   -------------------------------------------------------------------*/

Static Void Inputcoordinates(lchain_)
long *lchain_;
{
  struct LOC_Inputcoordinates V;
  Char datestring[30];
  long i, j;
  boolean finish;
  structure s;
  cardtype ctype;
  cardcontents cardinfo;
  long cardhist[(long)othercard - (long)headercard + 1];

  /***/

  V.lchain = lchain_;
  nss = 0;
  V.latom = 0;
  V.hatoms = 0;
  for (j = 0; j <= 5; j++)
    V.reszero.aaident[j] = ' ';
  V.reszero.aa = '!';
  V.reszero.access = 0;
  StrCopy(V.reszero.threelettercode, "    ", sizeof(char4));
  for (s = symbol; (long)s <= (long)beta2; s = (structure)((long)s + 1))
    V.reszero.ss[(long)s - (long)symbol] = ' ';
  V.reszero.sheetlabel = ' ';
  V.reszero.partner[0] = 0;
  V.reszero.partner[(long)beta2 - (long)beta1] = 0;
  V.reszero.alpha = 360.0;
  V.reszero.kappa = 360.0;
  for (j = 0; j <= 1; j++) {
    V.reszero.acceptor[j].residue = 0;
    V.reszero.acceptor[j].energy = 0;
    V.reszero.donor[j].residue = 0;
    V.reszero.donor[j].energy = 0;
  }
  V.reszero.atompointer = 0;
  V.reszero.nsideatoms = 0;
  for (j = 0; j <= 2; j++) {
    V.reszero.h[j] = 0.0;
    V.reszero.n[j] = 0.0;
    V.reszero.ca[j] = 0.0;
    V.reszero.c[j] = 0.0;
    V.reszero.o[j] = 0.0;
  }
  for (i = 0; i <= NMAX; i++)
    chain[i] = V.reszero;
  Date(datestring);   /* DATE(DAY-MONTH-YEAR); */
  /* comment out this line if necessary */
  fprintf(tapeout, "**** SECONDARY STRUCTURE DEFINITION ");
  fprintf(tapeout, "BY THE PROGRAM DSSP, VERSION OCT. 1988 ****");
  fprintf(tapeout, " DATE=%.11s", datestring);
  for (i = 106; i <= 127; i++)
    putc(' ', tapeout);
  fprintf(tapeout, ".\n");
  fprintf(tapeout, "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS ");
  fprintf(tapeout, "22 (1983) 2577-2637");
  for (i = 66; i <= 127; i++)
    putc(' ', tapeout);
  fprintf(tapeout, ".\n");
  for (ctype = headercard;
       (long)ctype <= (long)othercard;
       ctype = (cardtype)((long)ctype + 1))
    cardhist[(long)ctype - (long)headercard] = 0;
  V.corelimit = false;
  finish = false;
  V.resinfo = V.reszero;
  V.nmissing = true;
  V.camissing = true;
  V.cmissing = true;
  V.omissing = true;
  do {
    Readcard(&cardinfo, &V);
    cardhist[(long)cardinfo.art - (long)headercard]++;
    switch (cardinfo.art) {

    case headercard:
    case compndcard:
    case sourcecard:
    case authorcard:
      if (cardhist[(long)cardinfo.art - (long)headercard] == 1) {
	for (i = 0; i <= 127; i++)
	  putc(cardinfo.UU.z[i], tapeout);
	putc('\n', tapeout);
      }
      break;

    case ssbondcard:
      nss++;
      for (i = 0; i <= 1; i++)
	StrCopy(ssbonds[nss - 1][i], cardinfo.UU.r[i], sizeof(char6));
      break;

    case atomcard:
      if (cardinfo.UU.U5.residuename != '-' &&
	  (cardinfo.UU.U5.altloc == 'A' || cardinfo.UU.U5.altloc == ' ')) {
	if (strncmp(V.resinfo.aaident, cardinfo.UU.U5.reseqnum, sizeof(char6))) {
	  Putresidue(&V);
	  V.resinfo.atompointer = V.latom;
	  StrCopy(V.resinfo.aaident, cardinfo.UU.U5.reseqnum, sizeof(char6));
	  V.resinfo.aa = cardinfo.UU.U5.residuename;
	  StrCopy(V.resinfo.threelettercode, cardinfo.UU.U5.aaname,
		 sizeof(char4));
	}
	Getresidue(cardinfo.UU.U5.atomname, cardinfo.UU.U5.coordinates, &V);
      }
      if (cardinfo.UU.U5.residuename == '-') {
	printf(" !!! RESIDUE ");
	for (i = 0; i <= 3; i++)
	  putchar(cardinfo.UU.U5.aaname[i]);
	for (i = 0; i <= 5; i++)
	  putchar(cardinfo.UU.U5.reseqnum[i]);
	printf(" HAS NONSTANDARD NAME.\n");
	printf("     RESIDUE WILL BE ");
	printf("IGNORED !!!\n");
      }
      if (cardinfo.UU.U5.altloc != 'A' && cardinfo.UU.U5.altloc != ' ') {
	printf(" !!! IN RESIDUE");
	for (i = 0; i <= 3; i++)
	  printf(" %c", cardinfo.UU.U5.aaname[i]);
	for (i = 0; i <= 5; i++)
	  putchar(cardinfo.UU.U5.reseqnum[i]);
	printf(" ALTERNATE LOCATION INDICATOR ");
	printf("IS %c AND\n", cardinfo.UU.U5.altloc);
	printf("     NOT BLANK OR A. ATOM ");
	printf("NAMED ");
	for (i = 0; i <= 3; i++)
	  putchar(cardinfo.UU.U5.atomname[i]);
	printf(" WILL BE IGNORED !!!\n\n");
      }
      break;

    case tercard:
      Putresidue(&V);
      break;

    case endcard:
      finish = true;
      Putresidue(&V);
      break;

    case othercard:
      /* blank case */
      break;
    }
  } while (!(V.corelimit || finish));
  if (V.corelimit) {
    printf(" !!! NUMBER OF ATOMS OR RESIDUES EXCEEDS ");
    printf("STORAGE CAPACITY !!!\n");
  }
  if (!Nochainbreak(*V.lchain, *V.lchain))
    (*V.lchain)--;
  if (V.hatoms > 0) {
    printf(" !!! %12ld HYDROGEN OR DEUTERIUM ATOMS WERE IGNORED\n", V.hatoms);
    printf("     IN THE CALCULATION OF SIDE CHAIN SOLVENT \n");
    printf("     ACCESSIBILITY !!!\n");
  }
  if (cardhist[0] < 1)
    printf(" !!! HEADER-CARD MISSING !!!\n");
  if (cardhist[(long)compndcard - (long)headercard] < 1)
    printf(" !!! COMPOUND-CARD MISSING !!!\n");
  if (cardhist[(long)sourcecard - (long)headercard] < 1)
    printf(" !!! SOURCE-CARD MISSING !!!\n");
  if (cardhist[(long)authorcard - (long)headercard] < 1)
    printf(" !!! AUTHOR CARD MISSING !!!\n");
  if (*V.lchain < 1) {
    printf(" !!! NO RESIDUE WITH COMPLETE BACKBONE !!!\n");
    longjmp(_JL99, 1);
  }
  if (V.latom == 0)
    printf(" !!! ALL SIDECHAIN COORDINATES MISSING !!!\n");
}  /* Inputcoordinates */

#undef MAXSIDEATOMS


/***/
/*--------------------------------------------------------------------*/

Static boolean Testbond(i, j)
long i, j;
{
  /* TESTBOND IS TRUE IF I IS DONOR[=NH] TO J, OTHERWISE FALSE */
  backbone *WITH;

  WITH = &chain[i];
  return (WITH->acceptor[0].residue == j && WITH->acceptor[0].energy < HBHIGH ||
	  WITH->acceptor[1].residue == j && WITH->acceptor[1].energy < HBHIGH);
}  /* Testbond */


/***/

Local boolean Testssbond(i, j)
long i, j;
{
  boolean ssbond;
  long k;

  ssbond = false;
  k = 1;
  if (!(Nochainbreak(i, i) & Nochainbreak(j, j)))
    return ssbond;
  while (!(ssbond || k > nss)) {
    ssbond = ((!strncmp(chain[i].aaident, ssbonds[k - 1][0], sizeof(char6)) &&
	       !strncmp(chain[j].aaident, ssbonds[k - 1][1], sizeof(char6))) ||
	      (!strncmp(chain[i].aaident, ssbonds[k - 1][1], sizeof(char6)) &&
	       !strncmp(chain[j].aaident, ssbonds[k - 1][0], sizeof(char6))));
    k++;
  }
  return ssbond;
}  /* Testssbond */


/***/
/*--------------------------------------------------------------------*/

Static Void Flagssbonds()
{
  boolean ssbond;
  Char cc;
  long i, j, ii, jj;
  double d;
  long FORLIM;
  backbone *WITH;
  long FORLIM1;

  /***/

  nssintra = 0;
  nssinter = 0;
  cc = '`';
  FORLIM = lchain - 2;
  for (i = 1; i <= FORLIM; i++) {
    if (chain[i].aa == 'C' && chain[i].nsideatoms > 1) {
      ii = chain[i].atompointer + 2;
      j = i + 1;
      do {
	j++;
	ssbond = false;
	if (chain[j].nsideatoms > 1 && chain[j].aa == 'C')
	  jj = chain[j].atompointer + 2;
	else
	  jj = 0;
	if (jj > 0)
	  ssbond = (Distance(sidechain[ii - 1], sidechain[jj - 1]) < SSDIST);
      } while (!(ssbond || j == lchain));
      if (ssbond & (!Testssbond(i, j))) {
	printf(" !!! ADDITIONAL SSBOND FOUND BETWEEN ");
	printf("RESIDUES ");
	Writeresidue(chain[i]);
	printf(" AND ");
	Writeresidue(chain[j]);
	printf(" !!!\n\n");
      }
    }
  }
  if (nss > 0) {
    FORLIM = lchain - 2;
    for (i = 1; i <= FORLIM; i++) {
      WITH = &chain[i];
      if (WITH->aa == 'C') {
	FORLIM1 = lchain;
	for (j = i + 2; j <= FORLIM1; j++) {
	  if (chain[j].aa == 'C') {
	    if (Testssbond(i, j)) {
	      if (cc == 'z') {
		printf(" !!! SS-BRIDGE LABEL RESTART AT a !!!\n");
		cc = '`';
	      }
	      cc++;
	      WITH->aa = cc;
	      chain[j].aa = cc;
	      if (Nochainbreak(i, j))
		nssintra++;
	      else
		nssinter++;
	      if (WITH->nsideatoms > 1) {
		if (chain[j].nsideatoms > 1) {
		  jj = chain[j].atompointer + 2;
		  ii = WITH->atompointer + 2;
		  d = Distance(sidechain[ii - 1], sidechain[jj - 1]);
		  if (d > SSDIST) {
		    printf(" !!! SSBOND DISTANCE IS%5.1f BETWEEN RESIDUES", d);
		    Writeresidue(chain[i]);
		    printf(" AND ");
		    Writeresidue(chain[j]);
		    printf(" !!!\n\n");
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if (nss != nssintra + nssinter)
    printf(" !!! ERROR IN SSBOND DATA RECORDS !!!\n");
}  /* Flagssbonds */


/***/
/*--------------------------------------------------------------------*/

Static Void Flagchirality()
{
  long i;
  double ckap, skap;
  long FORLIM;
  backbone *WITH;

  FORLIM = lchain - 2;
  for (i = 2; i <= FORLIM; i++) {
    WITH = &chain[i];
    if (Nochainbreak(i - 1, i + 2)) {
      WITH->alpha = Dihedralangle(chain[i - 1].ca, WITH->ca, chain[i + 1].ca,
				  chain[i + 2].ca);
      if (WITH->alpha < 0.0)
	WITH->ss[(long)chirality - (long)symbol] = '-';
      else
	WITH->ss[(long)chirality - (long)symbol] = '+';
    }
  }
  FORLIM = lchain - 2;
  /***/
  for (i = 3; i <= FORLIM; i++) {
    WITH = &chain[i];
    if (Nochainbreak(i - 2, i + 2)) {
      ckap = Cosangle(chain[i].ca, chain[i - 2].ca, chain[i + 2].ca,
		      chain[i].ca);
      skap = sqrt(1.0 - ckap * ckap);
      WITH->kappa = RADIAN * Atan2(skap, ckap);
    }
  }
}  /* Flagchirality */


/***/

Local long Bondenergy(i, j)
long i, j;
{
  /*RESIDUE I IS DONOR[=NH],J IS ACCEPTOR[=CO] OF THE PROTON IN THE
     HYDROGEN BOND. THE BONDENERGY IS IN CAL/MOL */
  double dho, dhc, dnc, dno;
  long hbe;
  backbone *WITH;

  hbe = 0;
  WITH = &chain[i];
  if (WITH->aa == 'P')
    return hbe;
  dho = Distance(WITH->h, chain[j].o);
  dhc = Distance(WITH->h, chain[j].c);
  dnc = Distance(WITH->n, chain[j].c);
  dno = Distance(WITH->n, chain[j].o);
  if (dho < DIST || dhc < DIST || dnc < DIST || dno < DIST)
    hbe = HBLOW;
  else
    hbe = (long)floor(Q / dho - Q / dhc + Q / dnc - Q / dno + 0.5);
  if (hbe > HBLOW)
    return hbe;
  printf(" !!! CONTACT BETWEEN RESIDUES ");
  Writeresidue(chain[i]);
  printf(" AND ");
  Writeresidue(chain[j]);
  printf("  TOO CLOSE !!!\n");
  hbe = HBLOW;
  return hbe;
}  /* Bondenergy */

/***/

Local Void Updatebonds(b, hb)
hydrogenbond *b;
hydrogenbond hb;
{
  if (hb.energy < b[0].energy) {
    b[1] = b[0];
    b[0] = hb;
  } else if (hb.energy < b[1].energy)
    b[1] = hb;
}  /* Updatebonds */

/***/

Local Void Setbonds(i, j)
long i, j;
{
  /*I IS NH, J IS CO*/
  hydrogenbond hb;

  hb.energy = Bondenergy(i, j);
  hb.residue = j;
  /* CO(J) IS ACCEPTOR OF NH(I) */
  Updatebonds(chain[i].acceptor, hb);
  hb.residue = i;
  Updatebonds(chain[j].donor, hb);
}  /* Setbond */


/***/
/*--------------------------------------------------------------------*/

Static Void Flaghydrogenbonds()
{
  long i, j, FORLIM;
  backbone *WITH;
  long FORLIM1;

  /***/

  FORLIM = lchain;
  for (i = 1; i <= FORLIM; i++) {
    if (Nochainbreak(i, i)) {
      WITH = &chain[i];
      FORLIM1 = lchain;
      for (j = i + 1; j <= FORLIM1; j++) {
	if (Nochainbreak(j, j)) {
	  if (Distance(WITH->ca, chain[j].ca) < CADIST) {
	    Setbonds(i, j);
	    if (j != i + 1)
	      Setbonds(j, i);
	  }
	}
      }
    }
  }
}  /* Flaghydrogenbonds */


/***/

Local Void Ladder(i, j, b)
long i, j;
bridgetyp b;
{
  long k;
  boolean found;
  bridge *WITH;

  found = false;
  k = 1;
  if (b == nobridge || i >= j)
    return;
  do {
    WITH = &bridgetable[k - 1];
    if (WITH->ib == 0) {
      WITH->ib = i;
      WITH->ie = i;
      WITH->jb = j;
      WITH->je = j;
      WITH->from = 0;
      WITH->towards = 0;
      WITH->btyp = b;
      nbridge++;
      found = true;
    } else {
      found = (WITH->btyp == b && i == WITH->ie + 1) & Nochainbreak(WITH->ie,
		i) & (((j == WITH->je + 1 && b == parallel) &
		       Nochainbreak(WITH->je, j)) | ((j == WITH->jb - 1 &&
			  b == antiparallel) & Nochainbreak(j, WITH->jb)));
/* p2c: dssp.p, line 1609: Note:
 * Line breaker spent 1.1+0.26 seconds, 3126 tries on line 1540 [251] */
      if (found) {
	WITH->ie++;
	if (b == parallel)
	  WITH->je++;
	else
	  WITH->jb--;
      } else {
	k++;
	if (k > MAXBRIDGE) {
	  printf(" !!! BRIDGETABLE OVERFLOW !!!\n");
	  longjmp(_JL99, 1);
	}
      }
    }
  } while (!found);   /* Ladder */
}

/***/

Local Void Testbridge(i)
long i;
{
  long j1, j2, j;
  bridgetyp b;

  /***/

  j1 = 0;
  j2 = 0;
  j = i + 3;
  if (!Nochainbreak(i - 1, i + 1))
    return;
  while (j2 == 0 && j < lchain) {
    if (Nochainbreak(j - 1, j + 1)) {
      if ((Testbond(i + 1, j) & Testbond(j, i - 1)) |
	  (Testbond(j + 1, i) & Testbond(i, j - 1)))
	b = parallel;
      else if ((Testbond(i + 1, j - 1) & Testbond(j + 1, i - 1)) |
	       (Testbond(j, i) & Testbond(i, j)))
	b = antiparallel;
      else
	b = nobridge;
      if (b != nobridge) {
	if (j1 == 0) {
	  j1 = j;
	  Ladder(i, j, b);
	} else if (j != j1) {
	  j2 = j;
	  Ladder(i, j, b);
	}
      }
    }
    j++;
  }
}  /* Testbridge */

/***/

Local Void Extendladder()
{
  long i, j, ib1, jb1, je1;
  boolean bulge;
  long FORLIM;
  bridge *WITH;
  long SET[11];

  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++) {
    WITH = &bridgetable[i - 1];
    j = i + 1;
    while (j <= nbridge && WITH->towards == 0) {
      ib1 = bridgetable[j - 1].ib;
      jb1 = bridgetable[j - 1].jb;
      je1 = bridgetable[j - 1].je;
      bulge = (Nochainbreak(WITH->ie, ib1) && ib1 - WITH->ie < 6 &&
	       bridgetable[j - 1].btyp == WITH->btyp &&
	       bridgetable[j - 1].from == 0);
      if (bulge) {
	switch (WITH->btyp) {

	case parallel:
	  bulge = (jb1 - WITH->je < 6 && ib1 - WITH->ie < 3 ||
		   jb1 - WITH->je < 3) & Nochainbreak(WITH->je, jb1);
	  break;

	case antiparallel:
	  bulge = (WITH->jb - je1 < 6 && ib1 - WITH->ie < 3 ||
		   WITH->jb - je1 < 3) & Nochainbreak(je1, WITH->jb);
	  break;
	}
      }
      if (bulge) {
	WITH->towards = j;
	bridgetable[j - 1].from = i;
      }
      j++;
    }
  }
  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++) {
    WITH = &bridgetable[i - 1];
    if (WITH->from == 0) {
      P_expset(WITH->linkset, 0L);
      j = i;
      do {
	P_addset(WITH->linkset, (int)j);
	j = bridgetable[j - 1].towards;
      } while (j != 0);
      j = WITH->towards;
      while (j != 0) {
	P_setcpy(bridgetable[j - 1].linkset, WITH->linkset);
	j = bridgetable[j - 1].towards;
      }
    }
  }
}  /* Extendladder */

/* Local variables for Sheet: */
struct LOC_Sheet {
  long ladderset[MAXBRIDGE / 32 + 2], sheetset[MAXBRIDGE / 32 + 2];
} ;

/***/

Local boolean Link(l1, l2)
long l1, l2;
{
  /* LINK IS TRUE IF THERE IS A COMMON RESIDUE IN LADDERS L1 AND L2 */
  long ib1, ie1, jb1, je1, ib2, ie2, jb2, je2;

  ib1 = bridgetable[l1 - 1].ib;
  ie1 = bridgetable[l1 - 1].ie;
  jb1 = bridgetable[l1 - 1].jb;
  je1 = bridgetable[l1 - 1].je;
  ib2 = bridgetable[l2 - 1].ib;
  ie2 = bridgetable[l2 - 1].ie;
  jb2 = bridgetable[l2 - 1].jb;
  je2 = bridgetable[l2 - 1].je;
  return (ie1 >= ib2 && ib1 <= ie2 || ie1 >= jb2 && ib1 <= je2 ||
	  je1 >= ib2 && jb1 <= ie2 || je1 >= jb2 && jb1 <= je2);
}  /* Link */

/***/

Local Void Findsheet(LINK)
struct LOC_Sheet *LINK;
{
  long l1, l2;
  boolean finish;
  long FORLIM, FORLIM1;

  /***/

  P_expset(LINK->sheetset, 0L);
  l1 = 0;
  if (*LINK->ladderset != 0L) {
    do {
      l1++;
    } while (!P_inset((int)l1, LINK->ladderset));
  }
  if (l1 > 0)
    P_setcpy(LINK->sheetset, bridgetable[l1 - 1].linkset);
  if (l1 <= 0)
    return;
  do {
    finish = true;
    FORLIM = nbridge;
    for (l1 = 1; l1 <= FORLIM; l1++) {
      if (P_inset((int)l1, LINK->sheetset)) {
	FORLIM1 = nbridge;
	for (l2 = 1; l2 <= FORLIM1; l2++) {
	  if (P_inset((int)l2, LINK->ladderset)) {
	    if (Link(l1, l2)) {
	      P_setunion(LINK->sheetset, LINK->sheetset,
			 bridgetable[l2 - 1].linkset);
	      P_setdiff(LINK->ladderset, LINK->ladderset,
			bridgetable[l2 - 1].linkset);
	      finish = false;
	    }
	  }
	}
      }
    }
  } while (!finish);   /* Findsheet */
}

/***/

Local Void Sheet()
{
  struct LOC_Sheet V;
  long asci, i, j;
  Char ccs;
  long SET[11];
  long FORLIM;
  bridge *WITH;

  /***/

  P_expset(V.ladderset, 0L);
  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++)
    P_addset(V.ladderset, (int)i);
  ccs = '@';
  asci = 64;
  while (*V.ladderset != 0L) {
    ccs++;
    if (ccs > 'z') {
      printf(" !!! SHEET LABEL RESTART AT A !!!\n");
      ccs = 'A';
    }
    Findsheet(&V);
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
      WITH = &bridgetable[i - 1];
      if (P_inset((int)i, V.sheetset) && WITH->from == 0) {
	if (asci == 90) {
	  printf(" !!! STRAND LABEL RESTART AT A !!!\n");
	  asci = 64;
	}
	asci++;
	if (WITH->btyp == parallel)
	  WITH->laddername = (Char)(asci + 32);
	else
	  WITH->laddername = (Char)asci;
	WITH->sheetname = ccs;
	P_setcpy(WITH->linkset, V.sheetset);
	j = WITH->towards;
	while (j != 0) {
	  bridgetable[j - 1].laddername = WITH->laddername;
	  bridgetable[j - 1].sheetname = WITH->sheetname;
	  P_setcpy(bridgetable[j - 1].linkset, V.sheetset);
	  j = bridgetable[j - 1].towards;
	}
      }
    }
  }
}  /* Sheet */

/***/

Local Void Markstrands()
{
  long i, j, l, ib0, ie0, jb0, je0;
  structure beta, betai, betaj;
  long iset[(long)beta2 - (long)beta1 + 1][9],
       jset[(long)beta2 - (long)beta1 + 1][9];
  Char cc;
  long FORLIM, FORLIM1;
  long SET[9];
  bridge *WITH;
  backbone *WITH1;
  long SET1[9];
  long SET2[3];
  long SET3[3];

  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++) {
    if (bridgetable[i - 1].from == 0) {
      j = i;
      for (beta = beta1;
	   (long)beta <= (long)beta2;
	   beta = (structure)((long)beta + 1)) {
	P_setcpy(iset[(long)beta - (long)beta1], P_expset(SET, 0L));
	P_setcpy(jset[(long)beta - (long)beta1], P_expset(SET, 0L));
      }
      ib0 = lchain;
      ie0 = 0;
      jb0 = lchain;
      je0 = 0;
      do {
	WITH = &bridgetable[j - 1];
	FORLIM1 = WITH->ie;
	for (l = WITH->ib; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  for (beta = beta1;
	       (long)beta <= (long)beta2;
	       beta = (structure)((long)beta + 1))
	    P_setcpy(iset[(long)beta - (long)beta1], P_setunion(SET1,
		       iset[(long)beta - (long)beta1],
		       P_addset(P_expset(SET, 0L),
				WITH1->ss[(long)beta - (long)symbol])));
	}
	FORLIM1 = WITH->je;
	for (l = WITH->jb; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  for (beta = beta1;
	       (long)beta <= (long)beta2;
	       beta = (structure)((long)beta + 1))
	    P_setcpy(jset[(long)beta - (long)beta1], P_setunion(SET1,
		       jset[(long)beta - (long)beta1],
		       P_addset(P_expset(SET, 0L),
				WITH1->ss[(long)beta - (long)symbol])));
	}
	if (WITH->ib < ib0)
	  ib0 = WITH->ib;
	if (WITH->ie > ie0)
	  ie0 = WITH->ie;
	if (WITH->jb < jb0)
	  jb0 = WITH->jb;
	if (WITH->je > je0)
	  je0 = WITH->je;
	j = WITH->towards;
      } while (j != 0);
      j = i;
      if (P_setequal(iset[0], P_addset(P_expset(SET2, 0L), ' ')))
	betai = beta1;
      else
	betai = beta2;
      if (P_setequal(jset[0], P_addset(P_expset(SET2, 0L), ' ')))
	betaj = beta1;
      else
	betaj = beta2;
      if ((!P_setequal(iset[(long)betai - (long)beta1],
		       P_addset(P_expset(SET2, 0L), ' '))) |
	  (!P_setequal(jset[(long)betaj - (long)beta1],
		       P_addset(P_expset(SET3, 0L), ' '))))
	printf(" !!! STRAND COLUMN OVERWRITTEN !!!\n");
      do {
	WITH = &bridgetable[j - 1];
	FORLIM1 = WITH->ie;
	for (l = WITH->ib; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  WITH1->ss[(long)betai - (long)symbol] = WITH->laddername;
	  if (WITH->btyp == parallel)
	    WITH1->partner[(long)betai - (long)beta1] = WITH->jb + l - WITH->ib;
	  else
	    WITH1->partner[(long)betai - (long)beta1] = WITH->je - l + WITH->ib;
	}
	FORLIM1 = WITH->je;
	for (l = WITH->jb; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  WITH1->ss[(long)betaj - (long)symbol] = WITH->laddername;
	  if (WITH->btyp == parallel)
	    WITH1->partner[(long)betaj - (long)beta1] = WITH->ib + l - WITH->jb;
	  else
	    WITH1->partner[(long)betaj - (long)beta1] = WITH->ie - l + WITH->jb;
	}
	j = WITH->towards;
      } while (j != 0);
      if (ib0 == ie0)
	cc = 'B';
      else
	cc = 'E';
      for (j = ib0; j <= ie0; j++) {
	WITH1 = &chain[j];
	if (WITH1->ss[0] != 'E')
	  WITH1->ss[0] = cc;
      }
      for (j = jb0; j <= je0; j++) {
	WITH1 = &chain[j];
	if (WITH1->ss[0] != 'E')
	  WITH1->ss[0] = cc;
      }
    }
  }
  FORLIM = nbridge;
  for (j = 0; j < FORLIM; j++) {
    WITH = &bridgetable[j];
    FORLIM1 = WITH->ie;
    for (l = WITH->ib; l <= FORLIM1; l++)
      chain[l].sheetlabel = WITH->sheetname;
    FORLIM1 = WITH->je;
    for (l = WITH->jb; l <= FORLIM1; l++)
      chain[l].sheetlabel = WITH->sheetname;
  }
}  /* Markstrands */


/***/
/*--------------------------------------------------------------------*/

Static Void Flagbridge()
{
  long i, FORLIM;
  bridge *WITH;

  /***/

  for (i = 0; i < MAXBRIDGE; i++) {
    WITH = &bridgetable[i];
    WITH->ib = 0;
    WITH->ie = 0;
    WITH->jb = 0;
    WITH->je = 0;
    WITH->btyp = nobridge;
  }
  nbridge = 0;
  FORLIM = lchain;
  for (i = 2; i < FORLIM; i++)
    Testbridge(i);
  if (nbridge <= 0)
    return;
  Extendladder();
  Sheet();
  Markstrands();
}  /* Flagbridge */


/***/

Local Void Flagsymbol()
{
  /* FLAGS ALPHA HELICES AND TURNS IN SYMBOL COLUMN */
  long i, j, k;
  Char cc;
  long nhset[9];
  structure turn;
  boolean empty;
  long FORLIM;
  backbone *WITH;

  P_addset(P_expset(nhset, 0L), '>');
  P_addset(nhset, 'X');
  FORLIM = lchain - 4;
  for (i = 2; i <= FORLIM; i++) {
    if (P_inset(chain[i - 1].ss[(long)turn4 - (long)symbol], nhset) &
	P_inset(chain[i].ss[(long)turn4 - (long)symbol], nhset)) {
      for (j = i; j <= i + 3; j++)
	chain[j].ss[0] = 'H';
    }
  }
  FORLIM = lchain - 3;
  for (i = 2; i <= FORLIM; i++) {
    if (P_inset(chain[i - 1].ss[(long)turn3 - (long)symbol], nhset) &
	P_inset(chain[i].ss[(long)turn3 - (long)symbol], nhset)) {
      empty = true;
      for (j = i; j <= i + 2; j++) {
	WITH = &chain[j];
	if (WITH->ss[0] != 'G' && WITH->ss[0] != ' ')
	  empty = false;
      }
      if (empty) {
	for (j = i; j <= i + 2; j++)
	  chain[j].ss[0] = 'G';
      }
    }
  }
  FORLIM = lchain - 5;
  for (i = 2; i <= FORLIM; i++) {
    if (P_inset(chain[i - 1].ss[(long)turn5 - (long)symbol], nhset) &
	P_inset(chain[i].ss[(long)turn5 - (long)symbol], nhset)) {
      empty = true;
      for (j = i; j <= i + 4; j++) {
	WITH = &chain[j];
	if (WITH->ss[0] != 'I' && WITH->ss[0] != ' ')
	  empty = false;
      }
      if (empty) {
	for (j = i; j <= i + 4; j++)
	  chain[j].ss[0] = 'I';
      }
    }
  }
  FORLIM = lchain;
  for (i = 2; i < FORLIM; i++) {
    WITH = &chain[i];
    if (WITH->ss[0] == ' ') {
      cc = ' ';
      j = 1;
      for (turn = turn3;
	   (long)turn <= (long)turn5;
	   turn = (structure)((long)turn + 1)) {
	j++;
	for (k = 1; k <= j; k++) {
	  if (i > k) {
	    if (P_inset(chain[i - k].ss[(long)turn - (long)symbol], nhset))
	      cc = 'T';
	  }
	}
      }
      if (cc == ' ')
	cc = WITH->ss[(long)bend - (long)symbol];
      WITH->ss[0] = cc;
    }
  }
}  /* Flagsymbol */


/***/
/*--------------------------------------------------------------------*/

Static Void Flagturn()
{
  long i, j, k;
  structure turn;
  Char cc;
  long FORLIM1;
  backbone *WITH;

  /***/

  k = 2;
  cc = '2';
  for (turn = turn3; (long)turn <= (long)turn5; turn = (structure)((long)turn + 1)) {
    k++;
    cc++;
    FORLIM1 = lchain - k;
    for (i = 1; i <= FORLIM1; i++) {
      if (Nochainbreak(i, i + k)) {
	if (Testbond(i + k, i)) {
	  chain[i + k].ss[(long)turn - (long)symbol] = '<';
	  for (j = 1; j < k; j++) {
	    WITH = &chain[i + j];
	    if (WITH->ss[(long)turn - (long)symbol] == ' ')
	      WITH->ss[(long)turn - (long)symbol] = cc;
	  }
	  WITH = &chain[i];
	  if (WITH->ss[(long)turn - (long)symbol] == '<')
	    WITH->ss[(long)turn - (long)symbol] = 'X';
	  else
	    WITH->ss[(long)turn - (long)symbol] = '>';
	}
      }
    }
  }
  FORLIM1 = lchain;
  for (i = 1; i <= FORLIM1; i++) {
    WITH = &chain[i];
    if (WITH->kappa != 360.0 && WITH->kappa > 70.0)
      WITH->ss[(long)bend - (long)symbol] = 'S';
  }
  Flagsymbol();
}  /* Flagturn */


/* Local variables for Flagaccess: */
struct LOC_Flagaccess {
  long np;
  vector p[NFACE];
  double wp[NFACE];
} ;

/* Local variables for Polyeder: */
struct LOC_Polyeder {
  struct LOC_Flagaccess *LINK;
} ;

/***/

Local Void Triangle(x1, x2, x3, level, LINK)
double *x1, *x2, *x3;
long level;
struct LOC_Polyeder *LINK;
{
  long k, level1;
  double xnorm;
  vector x4, x5, x6;

  if (level > 0) {
    level1 = level - 1;
    for (k = 0; k <= 2; k++) {
      x4[k] = x1[k] + x2[k];
      x5[k] = x2[k] + x3[k];
      x6[k] = x1[k] + x3[k];
    }
    Norm(x4, &xnorm);
    Norm(x5, &xnorm);
    Norm(x6, &xnorm);
    Triangle(x1, x4, x6, level1, LINK);
    Triangle(x4, x2, x5, level1, LINK);
    Triangle(x4, x5, x6, level1, LINK);
    Triangle(x5, x3, x6, level1, LINK);
    return;
  }
  for (k = 0; k <= 2; k++)
    x6[k] = x1[k] + x2[k] + x3[k];
  Norm(x6, &xnorm);
  LINK->LINK->np++;
  VecCopy(LINK->LINK->p[LINK->LINK->np - 1], x6);
  Diff(x3, x1, x5);
  Diff(x2, x1, x4);
  Cross(x5, x4, x6);
  Norm(x6, &xnorm);
  LINK->LINK->wp[LINK->LINK->np - 1] = xnorm / 2.0;
}  /* Triangle */

/***/

Local Void Polyeder(LINK)
struct LOC_Flagaccess *LINK;
{  /* GENERATES ALL 12 VERTICES OF ICOSAHEDRON */
  struct LOC_Polyeder V;
  vector v[12];
  double a, b;
  long i, j, k, level, FORLIM;

  /***/

  V.LINK = LINK;
  k = 0;
  a = YVERTEX;
  b = ZVERTEX;
  for (i = 1; i <= 2; i++) {
    a = -a;
    for (j = 1; j <= 2; j++) {
      b = -b;
      k++;
      v[k - 1][0] = 0.0;
      v[k - 1][1] = a;
      v[k - 1][2] = b;
      k++;
      v[k - 1][0] = b;
      v[k - 1][1] = 0.0;
      v[k - 1][2] = a;
      k++;
      v[k - 1][0] = a;
      v[k - 1][1] = b;
      v[k - 1][2] = 0.0;
    }
  }
  LINK->np = 0;
  level = ORDER;
  /* GET ALL 20 FACES OF ICOSAHEDRON */
  for (i = 0; i <= 9; i++) {   /* FIND INTEGRATION POINTS */
    for (j = i + 1; j <= 10; j++) {
      if (Distance(v[i], v[j]) < 1.1) {
	for (k = j + 1; k <= 11; k++) {
	  if ((Distance(v[i], v[k]) < 1.1) & (Distance(v[j], v[k]) < 1.1))
	    Triangle(v[i], v[j], v[k], level, &V);
	}
      }
    }
  }
  a = 0.0;
  FORLIM = LINK->np;
  for (i = 0; i < FORLIM; i++)
    a += LINK->wp[i];
  a = FOURPI / a;
  FORLIM = LINK->np;
  for (i = 0; i < FORLIM; i++)
    LINK->wp[i] *= a;
}  /* Polyeder (eugiM oremodla) */

/* Local variables for Surface: */
struct LOC_Surface {
  struct LOC_Flagaccess *LINK;
  long nx;
  vector x[MAXPACK];
  double rx[MAXPACK];
} ;

/***/

Local boolean Step(xx, LINK)
double *xx;
struct LOC_Surface *LINK;
{
  long k;
  boolean one;
  double TEMP;

  one = true;
  k = 1;
  while (k <= LINK->nx && one) {
    TEMP = LINK->rx[k - 1] + RWATER;
    if (Distsq(xx, LINK->x[k - 1]) < TEMP * TEMP)
      one = false;
    else
      k++;
  }
  return one;
}  /* Step */

/* Local variables for Liste: */
struct LOC_Liste {
  struct LOC_Surface *LINK;
} ;

/***/

Local Void Listentry(xx, yy, d, r, LINK)
double *xx, *yy;
double d, r;
struct LOC_Liste *LINK;
{
  vector zz;
  double delta;

  delta = Distance(xx, yy);
  if (delta >= d + r)
    return;
  if (delta <= EPS)
    return;
  LINK->LINK->nx++;
  if (LINK->LINK->nx > MAXPACK) {
    printf(" !!! TABLE OVERFLOW IN FLAGACCESS !!!\n");
    longjmp(_JL99, 1);
    return;
  }
  Diff(yy, xx, zz);
  VecCopy(LINK->LINK->x[LINK->LINK->nx - 1], zz);
  LINK->LINK->rx[LINK->LINK->nx - 1] = r;
}  /* Listentry */

/***/

Local Void Liste(xx, rxx, LINK)
double *xx;
double rxx;
struct LOC_Surface *LINK;
{
  struct LOC_Liste V;
  long i, k;
  double d;
  long FORLIM;
  backbone *WITH;
  long FORLIM1;

  /***/

  V.LINK = LINK;
  LINK->nx = 0;
  d = rxx + RWATER + RWATER;
  FORLIM = lchain;
  for (i = 1; i <= FORLIM; i++) {
    if (Nochainbreak(i, i)) {
      WITH = &chain[i];
      if (Distance(xx, WITH->ca) < d + RESRAD) {
	Listentry(xx, WITH->n, d, RN, &V);
	Listentry(xx, WITH->ca, d, RCA, &V);
	Listentry(xx, WITH->c, d, RC, &V);
	Listentry(xx, WITH->o, d, RO, &V);
	if (WITH->nsideatoms > 0) {
	  FORLIM1 = WITH->nsideatoms;
	  for (k = 0; k < FORLIM1; k++)
	    Listentry(xx, sidechain[WITH->atompointer + k], d, RSIDEATOM, &V);
	}
      }
    }
  }
}  /* Liste */

/***/

Local double Surface(xatom, ratom, LINK)
double *xatom;
double ratom;
struct LOC_Flagaccess *LINK;
{
  struct LOC_Surface V;
  long i, j;
  double f, radius;
  vector xx;
  long FORLIM;

  /***/

  V.LINK = LINK;
  Liste(xatom, ratom, &V);
  radius = ratom + RWATER;
  f = 0.0;
  FORLIM = LINK->np;
  for (i = 0; i < FORLIM; i++) {
    for (j = 0; j <= 2; j++)
      xx[j] = LINK->p[i][j] * radius;
    if (Step(xx, &V))
      f += LINK->wp[i];
  }
  return (radius * radius * f);
}  /* Surface */


/***/
/*--------------------------------------------------------------------*/

Static Void Flagaccess()
{
  struct LOC_Flagaccess V;
  long i, k;
  double f;
  long FORLIM;
  backbone *WITH;
  long FORLIM1;

  /***/

  Polyeder(&V);
  FORLIM = lchain;
  for (i = 1; i <= FORLIM; i++) {
    if (Nochainbreak(i, i)) {
      WITH = &chain[i];
      f = Surface(WITH->n, RN, &V) + Surface(WITH->ca, RCA, &V) +
	  Surface(WITH->c, RC, &V) + Surface(WITH->o, RO, &V);
      if (WITH->nsideatoms > 0) {
	FORLIM1 = WITH->nsideatoms;
	for (k = 0; k < FORLIM1; k++)
	  f += Surface(sidechain[WITH->atompointer + k], RSIDEATOM, &V);
      }
      WITH->access = (long)floor(f + 0.5);
    }
  }
}  /* Flagaccess */


/***/

Local Void Statistics()
{
  long i, j, k, nchain, nres, nhbond, lhelix;
  bridgetyp b;
  Char cc;
  double Surface;
  long nhbturn[11];
  long ladderset[MAXBRIDGE / 32 + 2];
  long hbridge[(long)antiparallel - (long)parallel + 1];
  long helixhist[MAXHIST], sheethist[MAXHIST];
  long betahist[(long)antiparallel - (long)parallel + 1][MAXHIST];
  long FORLIM, FORLIM1;
  backbone *WITH;
  bridge *WITH1;
  long SET[11];
  long SET1[257];

  lhelix = 0;
  nhbond = 0;
  nchain = 0;
  nres = 0;
  for (i = 0; i < MAXHIST; i++) {
    for (b = parallel; (long)b <= (long)antiparallel; b = (bridgetyp)((long)b + 1))
      betahist[(long)b - (long)parallel][i] = 0;
    helixhist[i] = 0;
    sheethist[i] = 0;
  }
  Surface = 0.0;
  for (k = 0; k <= 10; k++)
    nhbturn[k] = 0;
  for (b = parallel; (long)b <= (long)antiparallel; b = (bridgetyp)((long)b + 1))
    hbridge[(long)b - (long)parallel] = 0;
  FORLIM = lchain;
  for (i = 0; i <= FORLIM; i++) {
    WITH = &chain[i];
    if (Nochainbreak(i, i)) {
      nres++;
      Surface += WITH->access;
      for (j = 0; j <= 1; j++) {
	if (WITH->donor[j].energy < HBHIGH) {
	  nhbond++;
	  k = WITH->donor[j].residue - i;
	  if (labs(k) < 6)
	    nhbturn[k + 5]++;
	}
      }
    } else
      nchain++;
    if (WITH->ss[0] == 'H')
      lhelix++;
    else if (lhelix > 0) {
      if (lhelix > MAXHIST)
	lhelix = MAXHIST;
      helixhist[lhelix - 1]++;
      lhelix = 0;
    }
  }
  if (nbridge > 0) {
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
      WITH1 = &bridgetable[i - 1];
      hbridge[(long)WITH1->btyp - (long)parallel] += WITH1->ie - WITH1->ib + 2;
      if (WITH1->from == 0) {
	j = i;
	k = 0;
	do {
	  k += bridgetable[j - 1].ie - bridgetable[j - 1].ib + 1;
	  j = bridgetable[j - 1].towards;
	} while (j != 0);
	if (k > MAXHIST)
	  k = MAXHIST;
	betahist[(long)WITH1->btyp - (long)parallel][k - 1]++;
      }
    }
  }
  if (nbridge > 0) {
    P_expset(ladderset, 0L);
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++)
      P_addset(ladderset, (int)i);
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
      WITH1 = &bridgetable[i - 1];
      if ((WITH1->from == 0) & P_inset((int)i, ladderset)) {
	if (!P_setequal(P_addset(P_expset(SET1, 0L), (int)i), WITH1->linkset) ||
	    WITH1->ie > WITH1->ib) {
	  k = 0;
	  FORLIM1 = nbridge;
	  for (j = 1; j <= FORLIM1; j++) {
	    if ((bridgetable[j - 1].from == 0) & P_inset((int)j, WITH1->linkset))
	      k++;
	  }
	  sheethist[k - 1]++;
	}
	P_setdiff(ladderset, ladderset, WITH1->linkset);
      }
    }
  }
  fprintf(tapeout,
    "%5ld%3ld%3ld%3ld%3ld TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .\n",
    nres, nchain, nssinter + nssintra, nssintra, nssinter);
  fprintf(tapeout,
    "%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .\n",
    Surface);
  fprintf(tapeout,
    "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .\n",
    nhbond, 100.0 * nhbond / nres);
  i = hbridge[0];
  j = hbridge[(long)antiparallel - (long)parallel];
  fprintf(tapeout,
    "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n",
    i, 100.0 * i / nres);
  fprintf(tapeout,
    "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n",
    j, 100.0 * j / nres);
  for (i = -5; i <= 5; i++) {
    if (i < 0)
      cc = '-';
    else
      cc = '+';
    k = labs(i);
    fprintf(tapeout,
      "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%ld), SAME NUMBER PER 100 RESIDUES                              .\n",
      nhbturn[i + 5], 100.0 * nhbturn[i + 5] / nres, cc, k);
  }
  for (i = 1; i <= MAXHIST; i++)
    fprintf(tapeout, "%3ld", i);
  fprintf(tapeout, "     *** HISTOGRAMS OF ***           .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", helixhist[i]);
  fprintf(tapeout, "    RESIDUES PER ALPHA HELIX         .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", betahist[0][i]);
  fprintf(tapeout, "    PARALLEL BRIDGES PER LADDER      .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", betahist[(long)antiparallel - (long)parallel][i]);
  fprintf(tapeout, "    ANTIPARALLEL BRIDGES PER LADDER  .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", sheethist[i]);
  fprintf(tapeout, "    LADDERS PER SHEET                .\n");
}  /* Statistics */

/***/

Local Void Writehb(i, hb)
long i;
hydrogenbond hb;
{
  double e;

  if (hb.residue != 0)
    hb.residue -= i;
  e = hb.energy / 1000.0;
  fprintf(tapeout, "%4ld,%4.1f", hb.residue, e);
}  /* Writehb */


/***/
/*--------------------------------------------------------------------*/

Static Void Printout()
{
  long i, j;
  structure s;
  double phi, psi, tco;
  long FORLIM;
  backbone *WITH;

  /***/

  Statistics();
  fprintf(tapeout,
    "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA \n");
  FORLIM = lchain;
  for (i = 1; i <= FORLIM; i++) {
    WITH = &chain[i];
    fprintf(tapeout, "%5ld ", i);
    for (j = 0; j <= 5; j++)
      putc(WITH->aaident[j], tapeout);
    fprintf(tapeout, " %c  %c ", WITH->aa, WITH->ss[0]);
    for (s = turn3; (long)s <= (long)beta2; s = (structure)((long)s + 1))
      putc(WITH->ss[(long)s - (long)symbol], tapeout);
    for (s = beta1; (long)s <= (long)beta2; s = (structure)((long)s + 1))
      fprintf(tapeout, "%4ld", WITH->partner[(long)s - (long)beta1]);
    fprintf(tapeout, "%c%4ld ", WITH->sheetlabel, WITH->access);
    for (j = 0; j <= 1; j++) {
      Writehb(i, WITH->acceptor[j]);
      Writehb(i, WITH->donor[j]);
    }
    phi = 360.0;
    psi = 360.0;
    tco = 0.0;
    if (Nochainbreak(i - 1, i)) {
      phi = Dihedralangle(chain[i - 1].c, WITH->n, WITH->ca, WITH->c);
      tco = Cosangle(WITH->c, WITH->o, chain[i - 1].c, chain[i - 1].o);
    }
    if (Nochainbreak(i, i + 1))
      psi = Dihedralangle(WITH->n, WITH->ca, WITH->c, chain[i + 1].n);
    fprintf(tapeout, "%8.3f%6.1f%6.1f%6.1f%6.1f%7.1f%7.1f%7.1f\n",
	    tco, WITH->kappa, WITH->alpha, phi, psi, WITH->ca[0], WITH->ca[1],
	    WITH->ca[2]);
  }
}  /* Printout */

Static Void Usage()
{
  fprintf(stderr,"Usage: dssp [-na] pdb_file dssp_file\n");
  fprintf(stderr,"the -na flag disables the calculation of accessible surface\n");
}

/***/
/*--------------------------------------------------------------------*/

main(argc, argv)
int argc;
Char *argv[];
{
  int tt; /*TIMELOCK*/
  PASCAL_MAIN(argc, argv);
  if (setjmp(_JL99))
    goto _L99;
  if(argc<3)
    {
      Usage();
      exit(-1);
    }
  while(*argv[1]=='-')
    {
      if(strcmp("-na",argv[1])==0)
	noaccFlag=1;
/*
      else if(strcmp("-s",argv[1])==0)
        silentFlag=1;
*/
      else
	{
	  Usage();
	  exit(-1);
	}
      argv++;
    }
  tapeout=fopen(argv[2],"w");
  tapein=fopen(argv[1],"r");
  rewind(tapein);
  tt=time(0); /*TIMELOCK*/
  printf(" \n");
  printf("                           DSSP\n");
  printf("            by Wolfgang Kabsch and Chris Sander\n");
  printf("\n");
  printf("Defines secondary structure and solvent exposure of proteins from\n");
  printf("atomic coordinates as given in Protein Data Bank format.   \n");
  printf("\n");
  printf("_________________________________________________________________________\n");
  printf("This version licensed to Ethan Benatan at Univ_Pittsburgh               \n");
  printf("for academic purposes.                                                 \n");
  printf("Do not redistribute. \n");
  printf("\n");
  printf("Commercial licenses available on request.                                \n");
  printf("\n");
  printf("Copyright by Wolfgang Kabsch and Chris Sander, 1983, 1985, 1988. \n");
  printf("Fax: +49-6221-387 306\n");
  printf("\n");
  printf("Algorithm version October 1988. Refer to Biopolymers 22(1983) 2577-2637.\n");
  printf("Do report errors if you find any.\n");
  printf("\n");
  printf("Email: Sander@embl-heidelberg.de \n");
  printf("       Kabsch@embl-heidelberg.de \n");
  printf("\n");
  printf("_________________________________________________________________________\n");
  printf("Related databases and datasets available from the Protein\n");
  printf("Design Group at EMBL via anonymous ftp from embl-heidelberg.de:\n");
  printf("\n");
  printf("pdb_select   Representative set of protein structures.\n");
  printf("             By U. Hobohm, C. Sander, M. Scharf and R. Schneider.\n");
  printf("             See Protein Science 1, 409-417.\n");
  printf("DSSP         Dictionary of secondary structures of proteins. \n");
  printf("HSSP         Database of sequence-homology derived protein families.\n");
  printf("             By C. Sander and R.Schneider.\n");
  printf("             See Proteins 9, 56-68 (1991).\n");
  printf("FSSP         Database of protein structure families with \n");
  printf("             common folding motifs. \n");
  printf("             L.Holm, C. Ouzounis, C. Sander, G.Tuparev, G. Vriend\n");
  printf("             See Protein Science, in the press (1992).\n");
  printf("In the XSSP databases, there is one dataset for each unique or\n");
  printf("             representative PDB protein, e.g., 1PPT.HSSP etc.\n");
  printf("\n");
  printf("Restrictions:Commercial users must apply for a license. \n");
  printf("             Not to be used for classified research.\n");
  
  lchain = 0;
  Inputcoordinates(&lchain);
  if (!Nochainbreak(1L, lchain))
    printf(" !!! POLYPEPTIDE CHAIN INTERRUPTED !!!\n");
  printf("INPUTCOORDINATES DONE%12ld\n", lchain);
  Flagssbonds();
  printf("FLAGSSBONDS DONE\n");
  Flagchirality();
  printf("FLAGCHIRALITY DONE\n");
  Flaghydrogenbonds();
  printf("FLAGHYDROGENBONDS DONE\n");
  Flagbridge();
  printf("FLAGBRIDGE DONE\n");
  Flagturn();
  printf("FLAGTURN DONE\n");
  if(noaccFlag==1)
    {
      printf("*ACCESSIBLE SURFACE *NOT* CALCULATED*\n");
    }
  else
    {
      Flagaccess();
      printf("FLAGACCESS DONE\n");
    }
  Printout();
  printf("PRINTOUT DONE\n");
_L99:
  if (tapein != NULL)
    fclose(tapein);
  if (tapeout != NULL)
    fclose(tapeout);
  exit(0);
}  /* END OF PROGRAM DSSP */

