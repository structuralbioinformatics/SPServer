#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

void rotateAtom(float, float, float, float *, float *, float *, float, float, float, int rev = 0);
void creatPDB(float, float, float, float, float, float, float, float, float, float, 
              float, float, float, float, float, int, int, int, int, float, int);

/******************************************************/
/* Function rotateAtom: rotates around 3 euler angles */
/******************************************************/

void rotateAtom (float oldX, float oldY, float oldZ, float *newX, float *newY, float *newZ, float psi, float theta, float phi, int rev) {
    float r11, r21, r31, r12, r22, r32, r13, r23, r33;
    if (rev == 0) {
        r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
        r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
        r31 = sin(theta)*sin(phi);
      
        r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
        r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
        r32 = sin(theta)*cos(phi);

        r13 = sin(psi)*sin(theta);
        r23 = -cos(psi)*sin(theta);
        r33 = cos(theta);
    } else { // if we are performing a reverse then need the transpose of the matrix
        r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
        r12 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
        r13 = sin(theta)*sin(phi);
      
        r21 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
        r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
        r23 = sin(theta)*cos(phi);

        r31 = sin(psi)*sin(theta);
        r32 = -cos(psi)*sin(theta);
        r33 = cos(theta);
    }

    *newX = r11 * oldX + r12 * oldY + r13 * oldZ;
    *newY = r21 * oldX + r22 * oldY + r23 * oldZ;
    *newZ = r31 * oldX + r32 * oldY + r33 * oldZ;
}

/*********************************************/
/* Function createPDB: create the ligand PDB */
/*********************************************/

void createPDB(float rrand1, float rrand2, float rrand3,
               float rand1, float rand2, float rand3,
               float r1, float r2, float r3, float l1, float l2, float l3, 
               float a1, float a2, float a3, int t1, int t2, int t3, int N,
               float spacing, int rot_rec) {
    string tmpbuf;
    float tx1, ty1, tz1, tx2, ty2, tz2, tx3, ty3, tz3;
    string x_string, y_string, z_string;
  
    while (getline(cin, tmpbuf)) {
        if ((tmpbuf.substr(0, 4) != "ATOM") && (tmpbuf.substr(0, 6) != "HETATM")) continue;
        float xcoord = atof(tmpbuf.substr(30, 8).c_str()) - l1;
        float ycoord = atof(tmpbuf.substr(38, 8).c_str()) - l2;
        float zcoord = atof(tmpbuf.substr(46, 8).c_str()) - l3;
      
        // rotate for initial randomization and the pose rotation
        rotateAtom(xcoord, ycoord, zcoord, &tx1, &ty1, &tz1, rand1, rand2, rand3);
        rotateAtom(tx1, ty1, tz1, &tx2, &ty2, &tz2, a1, a2, a3);
        // adjust so coordinates are in the box
        if (t1 >= N/2) t1 -= N;
        if (t2 >= N/2) t2 -= N;  
        if (t3 >= N/2) t3 -= N;
        // rot_rec == 1, need to rotate to the receptor's frame before translating back to original
        if (rot_rec == 0) {
            tx3 = tx2 - t1*spacing + r1;
            ty3 = ty2 - t2*spacing + r2;
            tz3 = tz2 - t3*spacing + r3;
        } else { 
            rotateAtom(tx2 - t1*spacing + r1, ty2 - t2*spacing + r2, tz2 - t3*spacing + r3, &tx3, &ty3, &tz3, rrand1, rrand2, rrand3, 1);
        }

        printf("%s%8.3f%8.3f%8.3f%s%s", tmpbuf.substr(0, 30).c_str(), tx3, ty3, tz3, tmpbuf.substr(54).c_str(), "\n");
    }
}

/**********************************************************************************************/
/* Function createPDBrev: create the ligand PDB based on the rotated and translated receptor! */
/**********************************************************************************************/

void createPDBrev(float rrand1, float rrand2, float rrand3,
                  float rand1, float rand2, float rand3,
                  float r1, float r2, float r3, float l1, float l2, float l3, 
                  float a1, float a2, float a3, int t1, int t2, int t3, int N,
                  float spacing) {
    string tmpbuf;
    float tx1, ty1, tz1, tx2, ty2, tz2, tx3, ty3, tz3;
    string x_string, y_string, z_string;
  
    while (getline(cin, tmpbuf)) {
        if ((tmpbuf.substr(0, 4) != "ATOM") && (tmpbuf.substr(0, 6) != "HETATM")) continue;
        float xcoord = atof(tmpbuf.substr(30, 8).c_str());
        float ycoord = atof(tmpbuf.substr(38, 8).c_str());
        float zcoord = atof(tmpbuf.substr(46, 8).c_str());
      
        // rotate the ligand
        rotateAtom(xcoord, ycoord, zcoord, &tx1, &ty1, &tz1, rrand1, rrand2, rrand3);
        // translate the ligand to the origin coords
        tx1 -= r1;
        ty1 -= r2;
        tz1 -= r3;
        // adjust so coordinates are in the box
        if (t1 >= N/2) t1 -= N;
        if (t2 >= N/2) t2 -= N;  
        if (t3 >= N/2) t3 -= N;  
        // translate to the pose position
        tx1 += t1*spacing;
        ty1 += t2*spacing;
        tz1 += t3*spacing;
        // rotate for initial pose rotation and initial rotation
        rotateAtom(tx1, ty1, tz1, &tx2, &ty2, &tz2, a1, a2, a3, 1);
        rotateAtom(tx2, ty2, tz2, &tx3, &ty3, &tz3, rand1, rand2, rand3, 1);
        tx3 += l1;
        ty3 += l2;
        tz3 += l3;
      
        printf("%s%8.3f%8.3f%8.3f%s%s", tmpbuf.substr(0, 30).c_str(), tx3, ty3, tz3, tmpbuf.substr(54).c_str(), "\n");
    }
}

/*****************/
/* MAIN FUNCTION */
/*****************/

int main(int argc, char **argv ) { 

    int t1, t2, t3;
    float r1, r2, r3, l1, l2, l3, a1, a2, a3, rand1, rand2, rand3, spacing;
    float recrand1 = 0.0, recrand2 = 0.0, recrand3 = 0.0;
    int N, reverse = -1, rot_rec = 0;
  
    if (argc == 18) {
        rand1 = atof(argv[1]);
        rand2 = atof(argv[2]);
        rand3 = atof(argv[3]);
        r1 = atof(argv[4]);
        r2 = atof(argv[5]);
        r3 = atof(argv[6]);
        l1 = atof(argv[7]);
        l2 = atof(argv[8]);
        l3 = atof(argv[9]);
        a1 = atof(argv[10]);
        a2 = atof(argv[11]);
        a3 = atof(argv[12]);
        t1 = atoi(argv[13]);
        t2 = atoi(argv[14]);
        t3 = atoi(argv[15]);
        N = atoi(argv[16]);
        spacing = atof(argv[17]);
    } else if (argc == 22) {
        reverse = atoi(argv[1]);
        recrand1 = atof(argv[2]);
        recrand2 = atof(argv[3]);
        recrand3 = atof(argv[4]);
        rand1 = atof(argv[5]);
        rand2 = atof(argv[6]);
        rand3 = atof(argv[7]);
        r1 = atof(argv[8]);
        r2 = atof(argv[9]);
        r3 = atof(argv[10]);
        l1 = atof(argv[11]);
        l2 = atof(argv[12]);
        l3 = atof(argv[13]);
        a1 = atof(argv[14]);
        a2 = atof(argv[15]);
        a3 = atof(argv[16]);
        t1 = atoi(argv[17]);
        t2 = atoi(argv[18]);
        t3 = atoi(argv[19]);
        N = atoi(argv[20]);
        spacing = atof(argv[21]);
        rot_rec = 1;
    } else {
        exit(1);
    }
  
    if ((reverse == -1) || (reverse == 0)) {
        createPDB(recrand1, recrand2, recrand3, rand1, rand2, rand3, r1, r2, r3, l1, l2, l3, a1, a2, a3, t1, t2, t3, N, spacing, rot_rec);
    } else { 
        createPDBrev(recrand1, recrand2, recrand3, rand1, rand2, rand3, r1, r2, r3, l1, l2, l3, a1, a2, a3, t1, t2, t3, N, spacing); // receptor is rotated in this case
    }
  return 0;
} 
