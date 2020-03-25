#include     "interact.h"

main (int argc, char *argv[])
{
 FILE *out;
 int  i,j,k,n,m,cont,help;
 int  ia[MAXATOM],ja[MAXATOM];
 char potencial[3],pdb_file_a[1000],pdb_file_b[1000],output[1000];
 PROT prota,protb;
 PROT readpdb();
 float cutoff,x,y,z,dx,dy,dz;
 float **d;
 float **Fmatrix();
 void  free_Fmatrix();


  memset(potencial,'\0',3);
  memset(pdb_file_a,'\0',1000);
  memset(pdb_file_b,'\0',1000);
  memset(output,'\0',1000);
  
  help=0;
  cont=0;
  for (i=0;i<argc;i++){
     if (strcmp(argv[i], "-a") == 0) { /* Protein A*/
         strcpy(pdb_file_a,argv[i+1]);
         cont++;
     }	
     if (strcmp(argv[i], "-b") == 0) { /* Protein B*/
         strcpy(pdb_file_b,argv[i+1]);
         cont++;
     }	
     if (strcmp(argv[i], "-c") == 0) { /* Cutoff distance of interactions*/
         sscanf(argv[i+1], "%f", &cutoff);
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
    printf("\t -a Protein A\n\t -b Protein B\n\t -c Cut-Off distance\n\t -p Potential type (ca,cb,min)\n\t -o Output file\n\t -h Help\n");
    exit(0);
  }
 

  prota=readpdb(pdb_file_a);
  protb=readpdb(pdb_file_b);

  d=Fmatrix(0,prota.number_of_res,0,protb.number_of_res);

  x=0.0;
  y=0.0;
  z=0.0;
  dx=0.0;
  dy=0.0;
  dz=0.0;

  for (i=0;i<prota.number_of_atoms;i++){
      if (!strcmp(prota.atoms_of_prot[i].name,"CA  ") && (!strcmp(potencial,"ca")||!strcmp(potencial,"CA")) ){
         x=prota.atoms_of_prot[i].x;
         y=prota.atoms_of_prot[i].y;
         z=prota.atoms_of_prot[i].z;
         ia[prota.atoms_of_prot[i].res_number]=i;
      }
      if (   (!strcmp(prota.atoms_of_prot[i].name,"CB  ") || (!strcmp(prota.atoms_of_prot[i].name,"CA  ") && !strcmp(prota.atoms_of_prot[i].res,"GLY")) )
          && (!strcmp(potencial,"cb")||!strcmp(potencial,"CB")) ){
         x=prota.atoms_of_prot[i].x;
         y=prota.atoms_of_prot[i].y;
         z=prota.atoms_of_prot[i].z;
         ia[prota.atoms_of_prot[i].res_number]=i;
      }
     //if (   ( (strcmp(prota.atoms_of_prot[i].name,"C   ") && strcmp(prota.atoms_of_prot[i].name,"N   ") && strcmp(prota.atoms_of_prot[i].name,"O   ") && strcmp(prota.atoms_of_prot[i].name,"CA  "))
     //       ||(!strcmp(prota.atoms_of_prot[i].res,"GLY")) )
     //    && (!strcmp(potencial,"min")||!strcmp(potencial,"MIN")) ){
      if (!strcmp(potencial,"min")||!strcmp(potencial,"MIN")) {
         x=prota.atoms_of_prot[i].x;
         y=prota.atoms_of_prot[i].y;
         z=prota.atoms_of_prot[i].z;
         ia[prota.atoms_of_prot[i].res_number]=i;
      }
      for (j=0;j<protb.number_of_atoms;j++){
          if (   !strcmp(protb.atoms_of_prot[j].name,"CA  ") && (!strcmp(potencial,"ca")||!strcmp(potencial,"CA"))  
              && !strcmp(prota.atoms_of_prot[i].name,"CA  ")){
             dx=protb.atoms_of_prot[j].x - x ;
             dy=protb.atoms_of_prot[j].y - y;
             dz=protb.atoms_of_prot[j].z - z;
             d[prota.atoms_of_prot[i].res_number][protb.atoms_of_prot[j].res_number]=sqrt(dx*dx + dy*dy + dz*dz);
             ja[protb.atoms_of_prot[j].res_number]=j;
          }
          if (   (!strcmp(protb.atoms_of_prot[j].name,"CB  ") || (!strcmp(protb.atoms_of_prot[j].name,"CA  ") && !strcmp(protb.atoms_of_prot[j].res,"GLY")) ) &&
                 (!strcmp(prota.atoms_of_prot[i].name,"CB  ") || (!strcmp(prota.atoms_of_prot[i].name,"CA  ") && !strcmp(prota.atoms_of_prot[i].res,"GLY")) )
              && (!strcmp(potencial,"cb")||!strcmp(potencial,"CB"))){
             dx=protb.atoms_of_prot[j].x - x ;
             dy=protb.atoms_of_prot[j].y - y;
             dz=protb.atoms_of_prot[j].z - z;
             d[prota.atoms_of_prot[i].res_number][protb.atoms_of_prot[j].res_number]=sqrt(dx*dx + dy*dy + dz*dz);
             ja[protb.atoms_of_prot[j].res_number]=j;
          }
         //if (   ( ( strcmp(protb.atoms_of_prot[j].name,"C   ") &&  strcmp(protb.atoms_of_prot[j].name,"N   ") &&  strcmp(protb.atoms_of_prot[j].name,"O   ") &&  strcmp(protb.atoms_of_prot[j].name,"CA  "))
         //       ||(!strcmp(protb.atoms_of_prot[j].res,"GLY")) ) &&
         //       ( (strcmp(prota.atoms_of_prot[i].name,"C   ") && strcmp(prota.atoms_of_prot[i].name,"N   ") && strcmp(prota.atoms_of_prot[i].name,"O   ") && strcmp(prota.atoms_of_prot[i].name,"CA  "))
         //       ||(!strcmp(prota.atoms_of_prot[i].res,"GLY")) )
         //    && (!strcmp(potencial,"min")||!strcmp(potencial,"MIN")) ){
          if (!strcmp(potencial,"min")||!strcmp(potencial,"MIN")) {
             dx=protb.atoms_of_prot[j].x - x ;
             dy=protb.atoms_of_prot[j].y - y;
             dz=protb.atoms_of_prot[j].z - z;
             if ( d[prota.atoms_of_prot[i].res_number][protb.atoms_of_prot[j].res_number] == 0.0) {
                  d[prota.atoms_of_prot[i].res_number][protb.atoms_of_prot[j].res_number] = sqrt(dx*dx + dy*dy + dz*dz);
                  ja[protb.atoms_of_prot[j].res_number]=j;
             }
             if ( sqrt(dx*dx + dy*dy +  dz*dz) < d[prota.atoms_of_prot[i].res_number][protb.atoms_of_prot[j].res_number]){
                  d[prota.atoms_of_prot[i].res_number][protb.atoms_of_prot[j].res_number]=sqrt(dx*dx + dy*dy + dz*dz);
                  ja[protb.atoms_of_prot[j].res_number]=j;
             }
          }
      }
  }

  out=fopen(output,"w");
  for (n=1;n<=prota.number_of_res;n++){
  for (m=1;m<=protb.number_of_res;m++){
   if (d[n][m]<=cutoff && d[n][m]>0.0){ fprintf(out,"%5d\t%5d\t%10.5e\t%c\t%c\t%s\t%s\n",n,m,d[n][m],prota.atoms_of_prot[ia[n]].res_code[0],protb.atoms_of_prot[ja[m]].res_code[0],prota.atoms_of_prot[ia[n]].name,protb.atoms_of_prot[ja[m]].name);}
  }}
  fclose(out);

}


 
