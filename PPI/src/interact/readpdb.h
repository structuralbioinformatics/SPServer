#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#define MAXS 200		/* Maximum number of characters in a line */
#define MAXATOM 100000	/* Maximum number of atoms in a protein */
#define MAXPROT 2

typedef struct input
  {
    char pdb_file[200];
  } INPUT;

typedef struct atom
  {
    char name[5];
    char res[5];
    char res_code[1];
    char chain[5];
    int  res_number;
    char a_res_number[6];
    float x;
    float y;
    float z;
    float bfact;
  } ATOM;


typedef struct prot
  {
    int number_of_atoms;
    int number_of_res;
    ATOM atoms_of_prot[MAXATOM];
  } PROT;

