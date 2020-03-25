#include     "interact.h"

PROT  readpdb(pdb_file)
char pdb_file[];
{

/*---------- Variables definitions ----------*/

  FILE *pdb;      	        /* Pointer to DataBase pdb file */
  char buffer[MAXS];
  char label[MAXS];
  char codigo[MAXS];
  PROT protein;
  int i;
  int cont;
  int j;
  int number;
  int ic;
  int ini;
  int fin;
  char a_name[5];		
  char a_res[5];
  char a_chain[5];
  char a_res_number[6];
  char seq[1];
  float a_x;
  float a_y;
  float a_z;
  float a_bfact;
  
  int pr_number_of_atoms;
  int pr_number_of_res;


/*---------- Initialize structures -----------*/

      for (j=0;j<MAXATOM;j++) {
	  memset(protein.atoms_of_prot[j].name,'\0',5);
          memset(protein.atoms_of_prot[j].res,'\0',5);
          memset(protein.atoms_of_prot[j].res_code,'\0',2);
          memset(protein.atoms_of_prot[j].a_res_number,'\0',6);
	  memset(protein.atoms_of_prot[j].chain,'\0',5);
          protein.atoms_of_prot[j].res_number = 0;
          protein.atoms_of_prot[j].x = 0;
	  protein.atoms_of_prot[j].y = 0;
	  protein.atoms_of_prot[j].z = 0;
	  protein.atoms_of_prot[j].bfact = 0;
	  }

/*---------- Open files ----------*/
/*---------- Running functions ----------*/

  cont=0;
   memset(buffer,'\0',MAXS);
   memset(label,'\0',MAXS);
   memset(a_name,'\0',5);
   memset(seq,'\0',1);
   memset(a_res,'\0',5);
   memset(a_res_number,'\0',6);
   memset(a_chain,'\0',5);
    if ((pdb=fopen(pdb_file, "r")) == NULL)			/* Open of pdb file */
    {
      printf("Can't open the pdb file %s \n",pdb_file);
      cont=1;
      /*
      printf("Enter a new name of pdb file: ");
      scanf("%s",&pdb_file);
      */
    }
  if (cont==0){ 
    printf("Open the pdb file %s \n",pdb_file);
    j=0;
    pr_number_of_atoms = 0;
    pr_number_of_res= 1;
    while (fgets(buffer,MAXS,pdb)!=NULL)
    {
      memset(label,'\0',MAXS);
      sscanf(buffer+0,"%s",label);
      if (strcmp(label,"ATOM") == 0)
      {
	memset(codigo,'\0',MAXS);
	sscanf(buffer+72,"%s",codigo);
        sscanf(buffer+13,"%c%c%c%c",&a_name[0],&a_name[1],&a_name[2],&a_name[3]);
	 /* printf("%s\n",buffer); */
        strcpy(protein.atoms_of_prot[j].name,a_name);
        sscanf(buffer+17,"%c%c%c",&a_res[0],&a_res[1],&a_res[2]);
        strcpy(protein.atoms_of_prot[j].res,a_res);
       if (strncmp(a_res, "ACE",3)){
        strncpy(seq,"X",1);
	if(!strncmp(a_res, "ALA",3))strncpy(seq,"A",1);
	if(!strncmp(a_res, "ARG",3))strncpy(seq,"R",1);
	if(!strncmp(a_res, "ASN",3))strncpy(seq,"N",1);
	if(!strncmp(a_res, "ASP",3))strncpy(seq,"D",1);
	if(!strncmp(a_res, "CYS",3))strncpy(seq,"C",1);
	if(!strncmp(a_res, "GLN",3))strncpy(seq,"Q",1);
	if(!strncmp(a_res, "GLU",3))strncpy(seq,"E",1);
	if(!strncmp(a_res, "GLY",3))strncpy(seq,"G",1);
	if(!strncmp(a_res, "HIS",3))strncpy(seq,"H",1);
	if(!strncmp(a_res, "ILE",3))strncpy(seq,"I",1);
	if(!strncmp(a_res, "LEU",3))strncpy(seq,"L",1);
	if(!strncmp(a_res, "LYS",3))strncpy(seq,"K",1);
	if(!strncmp(a_res, "MET",3))strncpy(seq,"M",1);
	if(!strncmp(a_res, "PHE",3))strncpy(seq,"F",1);
	if(!strncmp(a_res, "PRO",3))strncpy(seq,"P",1);
	if(!strncmp(a_res, "SER",3))strncpy(seq,"S",1);
	if(!strncmp(a_res, "THR",3))strncpy(seq,"T",1);
	if(!strncmp(a_res, "TRP",3))strncpy(seq,"W",1);
	if(!strncmp(a_res, "TYR",3))strncpy(seq,"Y",1);
	if(!strncmp(a_res, "VAL",3))strncpy(seq,"V",1);
        memset(protein.atoms_of_prot[j].res_code,'\0',2);
        strncpy(protein.atoms_of_prot[j].res_code,seq,1);
        sscanf(buffer+20,"%c%c",&a_chain[0],&a_chain[1]);
	if (a_chain[0]==' ' || a_chain[1]==' '){
	    if (a_chain[0]==' '&& a_chain[1]!=' ') {a_chain[0]=a_chain[1];a_chain[1]='\0';}
	    if (a_chain[0]==' '&& a_chain[1]==' ') {a_chain[0]=a_chain[1];a_chain[1]='\0';}
        }
        strcpy(protein.atoms_of_prot[j].chain,a_chain);
        sscanf(buffer+22,"%c%c%c%c%c",&a_res_number[0],&a_res_number[1],&a_res_number[2],&a_res_number[3],&a_res_number[4]);
        strcpy(protein.atoms_of_prot[j].a_res_number, a_res_number);
	if (j>0){
	if (strncmp(protein.atoms_of_prot[j-1].a_res_number,protein.atoms_of_prot[j].a_res_number,5) ) {
	 pr_number_of_res++;
	}}
	protein.atoms_of_prot[j].res_number=pr_number_of_res;
        sscanf(buffer+30,"%f",&a_x);
        protein.atoms_of_prot[j].x = a_x;
        sscanf(buffer+38,"%f",&a_y);
        protein.atoms_of_prot[j].y = a_y;
        sscanf(buffer+46,"%f",&a_z);
        protein.atoms_of_prot[j].z = a_z;
        sscanf(buffer+61,"%f",&a_bfact);
        protein.atoms_of_prot[j].bfact = a_bfact;
        pr_number_of_atoms++;
        j++;
       }
      }
      memset(buffer,'\0',MAXS);
    }
    fclose(pdb);					/* Close the PDB files */
    protein.number_of_atoms=pr_number_of_atoms;
    protein.number_of_res=pr_number_of_res;
    printf("Close the pdb file %s Total of Atoms=%d Total of Residues=%d\n",pdb_file,pr_number_of_atoms,pr_number_of_res); 
   }
 
    return protein;

}

