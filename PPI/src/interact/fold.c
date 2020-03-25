#include     "interact.h"

main (int argc, char *argv[])
{
 FILE *out;
 int  i,j,ic,gap,n,m,cont,help,skp;
 int  ia[MAXATOM],ja[MAXATOM],nc[MAXCHAIN];
 char potencial[3],pdb_file[1000],output[1000];
 PROT prot;
 PROT readpdb();
 float cutoff,x,y,z,dx,dy,dz;
 float **d;
 float **Fmatrix();
 void  free_Fmatrix();


  memset(potencial,'\0',3);
  memset(pdb_file,'\0',1000);
  memset(output,'\0',1000);
  
  help=0;
  cont=0;
  for (i=0;i<argc;i++){
     if (strcmp(argv[i], "-i") == 0) { /* Protein A*/
         strcpy(pdb_file,argv[i+1]);
         cont++;
     }	
     if (strcmp(argv[i], "-c") == 0) { /* Cutoff distance of interactions*/
         sscanf(argv[i+1], "%f", &cutoff);
         cont++;
     }	
     if (strcmp(argv[i], "-g") == 0) { /* Cutoff distance of interactions*/
         sscanf(argv[i+1], "%d", &gap);
         cont++;
     }	
     if (strcmp(argv[i], "-p") == 0) { /* Potential type : ca, cb, min */
         strcpy(potencial,argv[i+1]);
         cont++;
     }
     if (strcmp(argv[i], "-o") == 0) { /* Output file*/
         strcpy(output,argv[i+1]);
         cont++;
     }	
     if (strcmp(argv[i], "-h") == 0) {help=1;}		/* HELP */
  }


  if (help==1 || cont<5){
    printf("\t -i Protein \n\t -g Gap of residues to take the interaction\n\t -c Cut-Off distance\n\t -p Potential type (ca,cb,min)\n\t -o Output file\n\t -h Help\n");
    exit(0);
  }
 
 
  prot=readpdb(pdb_file);

  d=Fmatrix(0,prot.number_of_res,0,prot.number_of_res);

  x=0.0;
  y=0.0;
  z=0.0;
  dx=0.0;
  dy=0.0;
  dz=0.0;
  ic=0;

  for (i=0;i<prot.number_of_atoms;i++){
      if (!strcmp(prot.atoms_of_prot[i].name,"CA  ") && (!strcmp(potencial,"ca")||!strcmp(potencial,"CA")) ){
         x=prot.atoms_of_prot[i].x;
         y=prot.atoms_of_prot[i].y;
         z=prot.atoms_of_prot[i].z;
      }
      if (   (!strcmp(prot.atoms_of_prot[i].name,"CB  ") || (!strcmp(prot.atoms_of_prot[i].name,"CA  ") && !strcmp(prot.atoms_of_prot[i].res,"GLY")) )
          && (!strcmp(potencial,"cb")||!strcmp(potencial,"CB")) ){
         x=prot.atoms_of_prot[i].x;
         y=prot.atoms_of_prot[i].y;
         z=prot.atoms_of_prot[i].z;
      }
      if (   ( (strcmp(prot.atoms_of_prot[i].name,"C   ") && strcmp(prot.atoms_of_prot[i].name,"N   ") && strcmp(prot.atoms_of_prot[i].name,"O   ") && strcmp(prot.atoms_of_prot[i].name,"CA  "))
             ||(!strcmp(prot.atoms_of_prot[i].res,"GLY")) )
          && (!strcmp(potencial,"min")||!strcmp(potencial,"min")) ){
         x=prot.atoms_of_prot[i].x;
         y=prot.atoms_of_prot[i].y;
         z=prot.atoms_of_prot[i].z;
      }
      if (strcmp(prot.atoms_of_prot[i].chain,prot.atoms_of_prot[i+1].chain)){
                 nc[ic]=prot.atoms_of_prot[i].res_number;
                 ic++;
                 if (ic==MAXCHAIN){printf("Too many chains\n");}
      }
      for (j=i+1;j<prot.number_of_atoms;j++){
        if ( abs(prot.atoms_of_prot[i].res_number - prot.atoms_of_prot[j].res_number) > gap || strcmp(prot.atoms_of_prot[i].chain,prot.atoms_of_prot[j].chain)){
          if (!strcmp(prot.atoms_of_prot[j].name,"CA  ") && (!strcmp(potencial,"ca")||!strcmp(potencial,"CA")) &&
              !strcmp(prot.atoms_of_prot[i].name,"CA  ")){
             dx=prot.atoms_of_prot[j].x - x ;
             dy=prot.atoms_of_prot[j].y - y;
             dz=prot.atoms_of_prot[j].z - z;
             d[prot.atoms_of_prot[i].res_number][prot.atoms_of_prot[j].res_number]=sqrt(dx*dx + dy*dy + dz*dz);
             ia[prot.atoms_of_prot[i].res_number]=i;
             ja[prot.atoms_of_prot[j].res_number]=j;
          }
          if (   (!strcmp(prot.atoms_of_prot[j].name,"CB  ") || (!strcmp(prot.atoms_of_prot[j].name,"CA  ") && !strcmp(prot.atoms_of_prot[j].res,"GLY")) ) &&
                 (!strcmp(prot.atoms_of_prot[i].name,"CB  ") || (!strcmp(prot.atoms_of_prot[i].name,"CA  ") && !strcmp(prot.atoms_of_prot[i].res,"GLY")) )
              && (!strcmp(potencial,"cb")||!strcmp(potencial,"CB")) ){
             dx=prot.atoms_of_prot[j].x - x ;
             dy=prot.atoms_of_prot[j].y - y;
             dz=prot.atoms_of_prot[j].z - z;
             d[prot.atoms_of_prot[i].res_number][prot.atoms_of_prot[j].res_number]=sqrt(dx*dx + dy*dy + dz*dz);
             ia[prot.atoms_of_prot[i].res_number]=i;
             ja[prot.atoms_of_prot[j].res_number]=j;
          }
          if (   ( ( strcmp(prot.atoms_of_prot[j].name,"C   ") &&  strcmp(prot.atoms_of_prot[j].name,"N   ") &&  strcmp(prot.atoms_of_prot[j].name,"O   ") &&  strcmp(prot.atoms_of_prot[j].name,"CA  "))
                 ||(!strcmp(prot.atoms_of_prot[j].res,"GLY")) ) &&
                 ( (strcmp(prot.atoms_of_prot[i].name,"C   ") && strcmp(prot.atoms_of_prot[i].name,"N   ") && strcmp(prot.atoms_of_prot[i].name,"O   ") && strcmp(prot.atoms_of_prot[i].name,"CA  "))
                 ||(!strcmp(prot.atoms_of_prot[i].res,"GLY")) )
              && (!strcmp(potencial,"min")||!strcmp(potencial,"min")) ){
             dx=prot.atoms_of_prot[j].x - x ;
             dy=prot.atoms_of_prot[j].y - y;
             dz=prot.atoms_of_prot[j].z - z;
             if ( d[prot.atoms_of_prot[i].res_number][prot.atoms_of_prot[j].res_number] == 0.0) {
                  d[prot.atoms_of_prot[i].res_number][prot.atoms_of_prot[j].res_number] = sqrt(dx*dx + dy*dy + dz*dz);
                  ia[prot.atoms_of_prot[i].res_number]=i;
                  ja[prot.atoms_of_prot[j].res_number]=j;
             }
             if ( sqrt(dx*dx + dy*dy +  dz*dz) < d[prot.atoms_of_prot[i].res_number][prot.atoms_of_prot[j].res_number]){
                  d[prot.atoms_of_prot[i].res_number][prot.atoms_of_prot[j].res_number]=sqrt(dx*dx + dy*dy + dz*dz);
                  ia[prot.atoms_of_prot[i].res_number]=i;
                  ja[prot.atoms_of_prot[j].res_number]=j;
             }
          }
        }
      }
  }

  out=fopen(output,"w");
  for (n=1;n<=prot.number_of_res;n++){
  skp=0;
  for (i=0;i<ic;i++){if (n==nc[i]){skp=1;break;}}
  for (m=n+1;m<=prot.number_of_res;m++){
   if (m>=n+gap || skp==1){
    if (d[n][m]<=cutoff && d[n][m]>0.0){ fprintf(out,"%5d\t%5d\t%10.5e\t%c\t%c\t%s\t%s\n",n,m,d[n][m],prot.atoms_of_prot[ia[n]].res_code[0],prot.atoms_of_prot[ja[m]].res_code[0],prot.atoms_of_prot[ia[n]].name,prot.atoms_of_prot[ja[m]].name);}
   }
  }}
  fclose(out);
 
}



