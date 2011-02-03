/*
 
   Program information
   
      ThING Is Not Gaussian is a rudimentary program to calculate
      matrix elements of gaussian contractions using recursion
      relations.
      
   Copyright notice

      Copyright 2011, 2011 Sebastian Wouters
      <sebastianwouters@gmail.com>
      
   Copyright permission
   
      This file is part of ThING Is Not Gaussian.

      ThING Is Not Gaussian is free software: you can redistribute
      it and/or modify it under the terms of the GNU General 
      Public License as published by the Free Software Foundation,
      either version 3 of the License, or (at your option) any
      later version.

      ThING Is Not Gaussian is distributed in the hope that it will
      be useful, but WITHOUT ANY WARRANTY; without even the implied
      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
      See the GNU General Public License for more details.

      You should have received a copy of the GNU General Public
      License along with ThING Is Not Gaussian. If not, see 
      <http://www.gnu.org/licenses/>.
      
   File information MxElem.h and MxElem.cpp
   
      Class to fetch & store the matrix elements of a chemical
      problem. It contains:
         - the overlap matrix elements
         - the one body matrix elements
         - the kinetic energy matrix elements (so that nuclear charge
           rescaling can be used in HF if desired)
         - the teo body matrix elements

      Basic operations concerning the overlap matrix (Lodwin tfo,
      canonical tfo, inverse overlap matrix) are available
      "to a certain extent".
    
*/

#include "preamble.h"
#include "input.h"
#include "R.h"
#include "MxElem.h"
#include "MxElemFiller.h"


/**
 * Constructor for the MxElem class
 * @param readin the problem to be solved
 */
MxElem::MxElem(input & readin){

   this->NOrbTot = CalcTotalNumberOfOrbitals(readin);
   allocate();

}


/**
 * Copy constructor for the MxElem class
 * @param tocopy the MxElem object to be copied
 */
MxElem::MxElem(MxElem & tocopy){

   NOrbTot = tocopy.gNOrbTot();
   allocate();
   for (int i=0; i<NOrbTot; i++){
      for (int j=i; j<NOrbTot; j++){
         setTelem(i,j,tocopy.gTelem(i,j));
         setKEseparate(i,j,tocopy.gKEseparate(i,j));
         setSoverlap(i,j,tocopy.gSoverlap(i,j));
         for (int k=i; k<NOrbTot; k++)
            for (int l=j; l<NOrbTot; l++)
               setVelem(i,j,k,l,tocopy.gVelem(i,j,k,l));
      }
   }

}


/**
 * Standard destructor
 */
MxElem::~MxElem(){

   delete [] Telem;
   delete [] KEseparate;
   delete [] Soverlap;
   for (int i=0; i<NOrbTot; i++){
      for (int x=0; x<NOrbTot-i+1; x++){
         for (int y=0; y<NOrbTot-i+1; y++){
            delete [] Velem[i][x][y];
         }
         delete [] Velem[i][x];
      }
      delete [] Velem[i];
   }
   delete [] Velem;

}


/**
 * Function to allocate all the necassary memory
 */
void MxElem::allocate(){

   Telem = new double[(NOrbTot*(NOrbTot+1))/2];
   KEseparate = new double[(NOrbTot*(NOrbTot+1))/2];
   Soverlap = new double[(NOrbTot*(NOrbTot+1))/2];
   Velem = new double***[NOrbTot];
   for (int i=0; i<NOrbTot; i++){
      Velem[i] = new double**[NOrbTot-i+1];
      for (int x=0; x<NOrbTot-i+1; x++){
         Velem[i][x] = new double*[NOrbTot-i+1];
         for (int y=0; y<NOrbTot-i+1; y++){
            Velem[i][x][y] = new double[NOrbTot-i-x+1];
         }
      }
   }

}


/**
 * Function to find the orbital angular momentum corresponding to a subshell
 * @param type the type of the subshell: s, p, d, f, g, h
 * @return L the angular momentum
 */
int MxElem::GetLofType(char type){

   int L = -1;

   switch (type){
      case 's':
      case 'S':
         L = 0;
         break;
      case 'p':
      case 'P':
         L = 1;
         break;
      case 'd':
      case 'D':
         L = 2;
         break;
      case 'f':
      case 'F':
         L = 3;
         break;
      case 'g':
      case 'G':
         L = 4;
         break;
      case 'h':
      case 'H':
         L = 5;
         break;
      case 'i':
      case 'I':
         L = 6;
         break;
      case 'j':
      case 'J':
         L = 7;
         break;
      default:
         cout << "Subshell orbitals don't go higher than J in this program. Given : " << type << "." << endl;
         break;
   }

   assert(L>-1);
   return L;

}


/**
 * Function to find the total number of orbitals corresponding to a problem
 * @param readin the problem to be solved
 * @return the number of different orbitals
 */
int MxElem::CalcTotalNumberOfOrbitals(input & readin){

   int counter = 0;
   int Ncores = readin.gNcores();

   for (int cnt=0; cnt<Ncores; cnt++){

      Gauss * atom = readin.gGaussInfo(cnt);
      int Ntypes = atom->gNtypes();
      
      for (int cnt2=0; cnt2<Ntypes; cnt2++){

         char type = atom->gtype(cnt2);
         int L = GetLofType(type);
         counter += ((L+1)*(L+2))/2;

      }

   }

   return counter;

}


/**
 * Getter of the total number of orbitals
 * @return the total number of orbitals
 */
int MxElem::gNOrbTot(){

   return NOrbTot;

}


/**
 * Getter of the one body matrix element (i|O|j)
 * @param i the first orbital
 * @param j the second orbital
 * @return (i|O|j)
 */
double MxElem::gTelem(int i, int j){

   if (i>j)
      return Telem[j + (i*(i+1))/2];

   return Telem[i + (j*(j+1))/2];

}


/**
 * Setter of the one body matrix element (i|O|j)
 * @param i the first orbital
 * @param j the second orbital
 * @param v the new value of the mx element
 */
void MxElem::setTelem(int i, int j, double v){

   if (i>j)
      Telem[j + (i*(i+1))/2] = v;
   else
      Telem[i + (j*(j+1))/2] = v;

}


/**
 * Getter of the kinetic energy matrix element (i|T|j)
 * @param i the first orbital
 * @param j the second orbital
 * @return (i|T|j)
 */
double MxElem::gKEseparate(int i, int j){

   if (i>j)
      return KEseparate[j + (i*(i+1))/2];

   return KEseparate[i + (j*(j+1))/2];

}


/**
 * Setter of the kinetic energy matrix element (i|T|j)
 * @param i the first orbital
 * @param j the second orbital
 * @param v the new value of the mx element
 */
void MxElem::setKEseparate(int i, int j, double v){

   if (i>j)
      KEseparate[j + (i*(i+1))/2] = v;
   else
      KEseparate[i + (j*(j+1))/2] = v;

}


/**
 * Getter of the overlap matrix element (i|j)
 * @param i the first orbital
 * @param j the second orbital
 * @return (i|j)
 */
double MxElem::gSoverlap(int i, int j){

   if (i>j)
      return Soverlap[j + (i*(i+1))/2];

   return Soverlap[i + (j*(j+1))/2];

}


/**
 * Setter of the overlap matrix element (i|j)
 * @param i the first orbital
 * @param j the second orbital
 * @param v the new value of the mx element
 */
void MxElem::setSoverlap(int i, int j, double v){

   if (i>j)
      Soverlap[j + (i*(i+1))/2] = v;
   else
      Soverlap[i + (j*(j+1))/2] = v;

}


/**
 * Getter of the two body matrix element (ij|V|kl)
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @return (ij|V|kl)
 */
double MxElem::gVelem(int i, int j, int k, int l){

   // i smallest
   if ((i<=j) && (i<=k) && (i<=l)){
      if (j<=l)
         return gVelemOK(i,j,k,l);
      else
         return gVelemOK(i,l,k,j);
   }
   //j smallest
   if ((j<=k) && (j<=l)){
      if (i<=k)
         return gVelemOK(j,i,l,k);
      else
         return gVelemOK(j,k,l,i);
   }
   //k smallest
   if (k<=l){
      if (j<=l)
         return gVelemOK(k,j,i,l);
      else
         return gVelemOK(k,l,i,j);
   }
   // l smallest
   if (i<=k)
      return gVelemOK(l,i,j,k);

   return gVelemOK(l,k,j,i);

}


/**
 * Setter of the two body matrix element (ij|V|kl)
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @param v the new value of the mx element
 */
void MxElem::setVelem(int i, int j, int k, int l, double v){

   // i smallest
   if ((i<=j) && (i<=k) && (i<=l)){
      if (j<=l)
         setVelemOK(i,j,k,l,v);
      else
         setVelemOK(i,l,k,j,v);
   } else {
      //j smallest
      if ((j<=k) && (j<=l)){
         if (i<=k)
            setVelemOK(j,i,l,k,v);
         else
            setVelemOK(j,k,l,i,v);
      } else {
         //k smallest
         if (k<=l){
            if (j<=l)
               setVelemOK(k,j,i,l,v);
            else
               setVelemOK(k,l,i,j,v);
         } else {
            // l smallest
            if (i<=k)
               setVelemOK(l,i,j,k,v);
            else
               setVelemOK(l,k,j,i,v);
         }
      }
   }

}


/**
 * Getter of the two body matrix element (ij|V|kl) if i<j<l and i<k
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @return (ij|V|kl)
 */
double MxElem::gVelemOK(int i, int j, int k, int l){

   return Velem[i][j-i][k-i][l-j];

}


/**
 * Setter of the two body matrix element (ij|V|kl) if i<j<l and i<k
 * @param i the first orbital
 * @param j the second orbital
 * @param k the third orbital
 * @param l the fourth orbital
 * @param v the new value of the mx element
 */
void MxElem::setVelemOK(int i, int j, int k, int l, double v){

   //cout << "NOrbTot = " << NOrbTot << " and (" << i << "," << j << "," << k << "," << l << ")" << endl;
   Velem[i][j-i][k-i][l-j] = v;

}


/**
 * Save the matrix elements
 */
void MxElem::Save(){

   ofstream fileout("mxelem.bin", ios::binary | ios::trunc);
   fileout.seekp(0);
   fileout.write((char*)Telem,sizeof(*(Telem))*((NOrbTot*(NOrbTot+1))/2));
   fileout.write((char*)KEseparate,sizeof(*(KEseparate))*((NOrbTot*(NOrbTot+1))/2));
   fileout.write((char*)Soverlap,sizeof(*(Soverlap))*((NOrbTot*(NOrbTot+1))/2));

   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         for (int k=i; k<NOrbTot; k++)
            fileout.write((char*)Velem[i][j-i][k-i],sizeof(*(Velem[i][j-i][k-i]))*(NOrbTot-j+1));

   fileout.close();

}


/**
 * Load the matrix elements
 */
void MxElem::Load(){

   ifstream filein("mxelem.bin", ios::binary | ios::in);
   filein.seekg(0, ios::beg);
   filein.read((char*)Telem,sizeof(*(Telem))*((NOrbTot*(NOrbTot+1))/2));
   filein.read((char*)KEseparate,sizeof(*(KEseparate))*((NOrbTot*(NOrbTot+1))/2));
   filein.read((char*)Soverlap,sizeof(*(Soverlap))*((NOrbTot*(NOrbTot+1))/2));

   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         for (int k=i; k<NOrbTot; k++)
            filein.read((char*)Velem[i][j-i][k-i],sizeof(*(Velem[i][j-i][k-i]))*(NOrbTot-j+1));

   filein.close();

}


/**
 * Remove the stored matrix elements from disk
 */
void MxElem::RemoveFromDisk(){

   system("rm mxelem.bin");

}


/**
 * Calculate the nuclear-nuclear potential energy of a problem
 * @param problem the problem to be solved
 */
double MxElem::NuclPotEn(input & problem){

   double energy = 0;
   int Ncores = problem.gNcores();

   for (int i=0; i<Ncores; i++){

      R Ri(*(problem.gvector(i)));
      int Zi = problem.gcore(i);

      for (int j=i+1; j<Ncores; j++){

         energy += Zi * problem.gcore(j) / sqrt(Ri.DistanceSquared(*(problem.gvector(j))));
      }
   }

   return energy;

}


/**
 * Construct the transformation matrix for canonical orthogonalisation
 * @param A pointer to the array where the transformation needs to be stored
 */
void MxElem::CanOrth(double * A){

   char jobz = 'V';
   char uplo = 'U';
   int N = NOrbTot;
   int lda = N;
   
   int lwork = 3*NOrbTot-1;
   double * work = new double[lwork];
   double * eigs = new double[N];

   int info;

   for (int i=0; i<N; i++)
      for (int j=i; j<N; j++)
         A[i+N*j] = gSoverlap(i,j);

   dsyev_(&jobz, &uplo, &N, A, &lda, eigs, work, &lwork, &info);

   int incx = 1;

   for (int i=0; i<N; i++){
      eigs[i] = 1.0/sqrt(eigs[i]);
      dscal_(&N, eigs+i, A+N*i, &incx);
   }

   delete [] eigs;
   delete [] work;

}


/**
 * Construct the inverse of the overlap matrix
 * @param Sinv pointer to the array where the inverse of the overlapmatrix should be stored
 */
void MxElem::InverseOverlap(double * Sinv){

   double * V = new double[NOrbTot*NOrbTot];
   CanOrth(V);

   // Sinv = V*V^T 
   char trans = 'T';
   char notrans = 'N';
   int m = NOrbTot;
   int k = NOrbTot;
   int n = NOrbTot;
   double alpha = 1.0;
   double beta = 0.0;
   int lda = NOrbTot;
   int ldb = NOrbTot;
   int ldc = NOrbTot;
   dgemm_(&notrans, &trans, &m, &n, &k, &alpha, V, &lda, V, &ldb, &beta, Sinv, &ldc);

   delete [] V;

}


/**
 * Construct the Lodwin transformation
 * @param Slodwin pointer to the array where the Lodwin transformation should be stored
 */
void MxElem::MakeLodwinTfo(double * Slodwin){

   double * V = new double[NOrbTot*NOrbTot];
   
   // Fill V with the overlap matrix
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         V[i+NOrbTot*j] = gSoverlap(i,j);
   
   // Eigenvectors -> V, eigenvalues -> eigs
   char jobz = 'V';
   char uplo = 'U';
   int N = NOrbTot;
   int lda = N;
   int lwork = 3*NOrbTot-1;
   double * work = new double[lwork];
   double * eigs = new double[N];
   int info;
   dsyev_(&jobz, &uplo, &N, V, &lda, eigs, work, &lwork, &info);

   // V = V*D^(-0.25)
   int incx = 1;
   for (int i=0; i<N; i++){
      eigs[i] = pow(eigs[i],-0.25);
      dscal_(&N, eigs+i, V+N*i, &incx);
   }

   delete [] eigs;
   delete [] work;

   // Slodwin = V*V^T 
   char trans = 'T';
   char notrans = 'N';
   int m = NOrbTot;
   int k = NOrbTot;
   int n = NOrbTot;
   double alpha = 1.0;
   double beta = 0.0;
   lda = NOrbTot;
   int ldb = NOrbTot;
   int ldc = NOrbTot;
   dgemm_(&notrans, &trans, &m, &n, &k, &alpha, V, &lda, V, &ldb, &beta, Slodwin, &ldc);

   delete [] V;

}


/**
 * Function that performs the Lodwin transformation
 */
void MxElem::DoLodwinTfo(){

   int NOrb2 = NOrbTot*NOrbTot;
   int NOrb3 = NOrbTot*NOrb2;
   int NOrb4 = NOrb2*NOrb2;
   
   //Step 1. Fetch the Lodwin transformation. Notice that it doesn't have to be symmetric: possibly first blockdiagonaltfo and then full one.
   double * Slodwin = new double[NOrb2];
   MakeLodwinTfo(Slodwin);
   
   PrintFull(Slodwin,NOrbTot);
   
   //Step 2. Set the overlapmatrix to 1.
   for (int i=0; i<NOrbTot; i++){
      setSoverlap(i,i,1.0);
      for (int j=i+1; j<NOrbTot; j++)
         setSoverlap(i,j,0.0);
   }
   
   //Step 3. Transform One Body Mx and KEseparate
   //Step 3.1. Fetch One Body Mx
   double * OneBody = new double[NOrb2];
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         OneBody[i+j*NOrbTot] = gTelem(i,j);
   
   //Step 3.2. OneBody*Slodwin -> mem
   double * mem = new double[NOrb2];
   
   char side = 'L';
   char uplo = 'U';

   int m = NOrbTot;
   int n = NOrbTot;
   int k = NOrbTot;
   double alpha = 1.0;
   double beta = 0.0;
   int lda = NOrbTot;
   int ldb = NOrbTot;
   int ldc = NOrbTot;
   dsymm_(&side, &uplo, &m, &n, &alpha, OneBody, &lda, Slodwin, &ldb, &beta, mem, &ldc);
   
   //Step 3.3. Slodwin^T*mem -> OneBody
   char trans = 'T';
   char notrans = 'N';
   dgemm_(&trans, &notrans, &m, &n, &k, &alpha, Slodwin, &lda, mem, &ldb, &beta, OneBody, &ldc);
   
   //Step 3.4. OneBody -> Telem
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         setTelem(i,j,OneBody[i+NOrbTot*j]);
         
   //Step 3.5. Fetch KEseparate
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         OneBody[i+j*NOrbTot] = gKEseparate(i,j);
   
   //Step 3.6. Do the tfo's
   dsymm_(&side, &uplo, &m, &n, &alpha, OneBody, &lda, Slodwin, &ldb, &beta, mem, &ldc);
   dgemm_(&trans, &notrans, &m, &n, &k, &alpha, Slodwin, &lda, mem, &ldb, &beta, OneBody, &ldc);
   
   //Step 3.7. OneBody -> KEseparate
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         setKEseparate(i,j,OneBody[i+NOrbTot*j]);
       
   delete [] OneBody;
   delete [] mem;
   
   //Step 4. Transform Two Body Mx
   //Step 4.1. Initialise a temporary storage
   double * mem1 = new double[NOrb4];
   double * mem2 = new double[NOrb4];

   //Step 4.2. Construct mem1
   for (int i=0; i<NOrbTot; i++)
      for (int j=0; j<NOrbTot; j++)
         for (int k=0; k<NOrbTot; k++)
            for (int l=0; l<NOrbTot; l++)
               mem1[i+NOrbTot*j+NOrb2*k+NOrb3*l] = //(ijkl)
               mem1[k+NOrbTot*j+NOrb2*i+NOrb3*l] = //(kjil)
               mem1[i+NOrbTot*l+NOrb2*k+NOrb3*j] = //(ilkj)
               mem1[k+NOrbTot*l+NOrb2*i+NOrb3*j] = //(klij)
               mem1[j+NOrbTot*i+NOrb2*l+NOrb3*k] = //(jilk)
               mem1[l+NOrbTot*i+NOrb2*j+NOrb3*k] = //(lijk)
               mem1[j+NOrbTot*k+NOrb2*l+NOrb3*i] = //(jkli)
               mem1[k+NOrbTot*l+NOrb2*j+NOrb3*i] = gVelem(i,j,k,l); //(klji)

   //Step 4.3. Do Slodwin^T * mem1 -> mem2 or (ij|V|kl) -> (aj|V|kl).
   m = NOrbTot;
   n = NOrb3;
   k = NOrbTot;
   lda = m;
   ldb = k;
   ldc = m;
   dgemm_(&trans, &notrans, &m, &n, &k, &alpha, Slodwin, &lda, mem1, &ldb, &beta, mem2, &ldc);
   
   //Step 4.4. Do mem2 * Slodwin -> mem1 or (aj|V|kl) -> (aj|V|kd).
   m = NOrb3;
   n = NOrbTot;
   k = NOrbTot;
   lda = m;
   ldb = k;
   ldc = m;
   dgemm_(&notrans, &trans, &m, &n, &k, &alpha, mem2, &lda, Slodwin, &ldb, &beta, mem1, &ldc);
   
   //Step 4.5. Do (mem1+N^3*i) * Slodwin -> (mem2 + N^3*i) or (aj|V|kd) -> (aj|V|cd).
   m = NOrb2;
   n = NOrbTot;
   k = NOrbTot;
   lda = m;
   ldb = k;
   ldc = m;
   for (int i=0; i<NOrbTot; i++)
      dgemm_(&notrans, &notrans, &m, &n, &k, &alpha, mem1+NOrb3*i, &lda, Slodwin, &ldb, &beta, mem2+NOrb3*i, &ldc);
      
   //Step 4.6. Do (mem2+N^2*i+N^3*j) * Slodwin -> (mem1+N^2*i+N^3*j) or (aj|V|cd) -> (ab|V|cd).
   m = NOrbTot;
   n = NOrbTot;
   k = NOrbTot;
   lda = m;
   ldb = k;
   ldc = m;
   for (int i=0; i<NOrbTot; i++)
      for (int j=0; j<NOrbTot; j++)
         dgemm_(&notrans, &notrans, &m, &n, &k, &alpha, mem2+NOrb2*i+NOrb3*j, &lda, Slodwin, &ldb, &beta, mem1+NOrb2*i+NOrb3*j, &ldc);
         
   //Step 4.7. Copy back
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         for (int k=i; k<NOrbTot; k++)
            for (int l=j; l<NOrbTot; l++)
               setVelem(i,j,k,l,mem1[i+NOrbTot*j+NOrb2*k+NOrb3*l]);
   
   delete [] mem1;
   delete [] mem2;
   
   //Step 5. Delete the remaining allocation: The Lodwin Tfo
   delete [] Slodwin;
   

}


/**
 * Initialise the matrix elements
 * @param readin the problem to be solved
 */
void MxElem::Init(input & readin){

   MxElemFiller filler(readin);
   int count1, count2, count3, count4, start, start2, start3;
   double OneBodyElement;
   
   int symmaxis = 0;
   if ((readin.gRotationSymm()=='x')||(readin.gRotationSymm()=='X')) symmaxis = 1;
   else if ((readin.gRotationSymm()=='y')||(readin.gRotationSymm()=='Y')) symmaxis = 2;
   else if ((readin.gRotationSymm()=='z')||(readin.gRotationSymm()=='Z')) symmaxis = 3;

   count1 = -1;

   for (int i=0; i<readin.gNcores(); i++){
      Gauss * first = readin.gGaussInfo(i);

      for (int k=0; k<(first->gNtypes()); k++){
         int L1 = GetLofType(first->gtype(k));

         // loop over different [n1x,n1y,n1z] contributions with n1x+n1y+n1z=L1;
         for (int n1x=L1; n1x>=0; n1x--){
            for (int n1y=L1-n1x; n1y>=0; n1y--){
               int n1z = L1-n1x-n1y;
               count1++;

               count2 = -1;

               for (int j=i; j<readin.gNcores(); j++){
                  Gauss * second = readin.gGaussInfo(j);
                  start = 0;
                  if (i==j) start=k;

                  for (int l=start; l<(second->gNtypes()); l++){
                     int L2 = GetLofType(second->gtype(l));

                     // loop over different [n2x,n2y,n2z] contributions with n2x+n2y+n2z=L2;
                     for (int n2x=L2; n2x>=0; n2x--){
                        for (int n2y=L2-n2x; n2y>=0; n2y--){
                           int n2z = L2-n2x-n2y;

                           if ((i==j)&&(l==k)){
                              if (n2x>n1x){
                                 n2x=n1x;
                                 n2y=n1y;
                                 n2z=n1z;
                              } else {
                                 if (n2x==n1x){
                                    if (n2y>n1y){
                                       n2y=n1y;
                                       n2z=n1z;
                                    }
                                 }
                              }
                           }

                           count2++;
                           //If they're centered on the same atom and n1i+n2i odd for i=x,y or z -> Overlap & KE 0.0
                           //If they're centered on different atoms and n1i+n2i odd for i!=RotationSymm -> Overlap & KE 0.0
                           if (i==j){
                              if ((((n1x+n2x)%2)!=0)||(((n1y+n2y)%2)!=0)||(((n1z+n2z)%2)!=0)){
                                 setSoverlap(count1,count1+count2,0.0);
                                 setKEseparate(count1,count1+count2,0.0);
                                 OneBodyElement = 0.0;
                              } else {
                                 setSoverlap(count1,count1+count2,filler.Overlap(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z));
                                 OneBodyElement = filler.KE(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z);
                                 setKEseparate(count1,count1+count2,OneBodyElement);
                              }
                           } else {
                              
                              if (  ((symmaxis==1)&&((((n1y+n2y)%2)!=0)||(((n1z+n2z)%2)!=0))) || 
                                    ((symmaxis==2)&&((((n1x+n2x)%2)!=0)||(((n1z+n2z)%2)!=0))) ||
                                    ((symmaxis==3)&&((((n1y+n2y)%2)!=0)||(((n1x+n2x)%2)!=0))) ){
                                 setSoverlap(count1,count1+count2,0.0);
                                 setKEseparate(count1,count1+count2,0.0);
                                 OneBodyElement = 0.0;
                              } else {
                                 setSoverlap(count1,count1+count2,filler.Overlap(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z));
                                 OneBodyElement = filler.KE(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z);
                                 setKEseparate(count1,count1+count2,OneBodyElement);
                              }
                              
                              
                           }
                           
                           for (int Ncore=0; Ncore<readin.gNcores(); Ncore++)
                              OneBodyElement += filler.ElNucl(i,k,n1x,n1y,n1z,j,l,n2x,n2y,n2z,Ncore);

                           setTelem(count1,count1+count2,OneBodyElement);

                           count3 = -1;

                           for (int m=i; m<readin.gNcores(); m++){
                              Gauss * third = readin.gGaussInfo(m);
                              start2 = 0;
                              if (i==m) start2=k;

                              for (int n=start2; n<(third->gNtypes()); n++){
                                 int L3 = GetLofType(third->gtype(n));
                                 for (int n3x=L3; n3x>=0; n3x--){
                                    for (int n3y=L3-n3x; n3y>=0; n3y--){
                                       int n3z = L3-n3x-n3y;

                                       if ((i==m)&&(n==k)){
                                          if (n3x>n1x){
                                             n3x=n1x;
                                             n3y=n1y;
                                             n3z=n1z;
                                          } else {
                                             if (n3x==n1x){
                                                if (n3y>n1y){
                                                   n3y=n1y;
                                                   n3z=n1z;
                                                }
                                             }
                                          }
                                       }

                                       count3++;

                                       count4 = -1;
                                       for (int o=j; o<readin.gNcores(); o++){
                                          Gauss * fourth = readin.gGaussInfo(o);
                                          start3 = 0;
                                          if (j==o) start3=l;

                                          for (int p=start3; p<(fourth->gNtypes()); p++){
                                             int L4 = GetLofType(fourth->gtype(p));
                                             for (int n4x=L4; n4x>=0; n4x--){
                                                for (int n4y=L4-n4x; n4y>=0; n4y--){
                                                   int n4z = L4-n4x-n4y;

                                                   if ((j==o)&&(p==l)){
                                                      if (n4x>n2x){
                                                         n4x=n2x;
                                                         n4y=n2y;
                                                         n4z=n2z;
                                                      } else {
                                                         if (n4x==n2x){
                                                            if (n4y>n2y){
                                                               n4y=n2y;
                                                               n4z=n2z;
                                                            }
                                                         }
                                                      }
                                                   }

                                                   count4++;

                                                   //alpha=(i,k,1) beta=(j,l,2) gamma=(m,n,3) delta=(o,p,4)
                                                   //(alpha beta|V|gamma delta) = (alpha gamma beta delta) in MxElemFiller!!
                                                   setVelem(count1, count1+count2, count1+count3, count1+count2+count4, filler.ElEl(i,k,n1x,n1y,n1z,m,n,n3x,n3y,n3z,j,l,n2x,n2y,n2z,o,p,n4x,n4y,n4z));
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

}


/**
 * Print the overlap matrix
 */
void MxElem::PrintS(){

   for (int i=0; i<NOrbTot; i++){
      for (int j=0; j<i; j++)
         cout << "X\t";
      for (int j=i; j<NOrbTot; j++)
         cout << gSoverlap(i,j) << "\t";
      cout << endl;
   }

}


/**
 * Print the one body matrix
 */
void MxElem::PrintOneBody(){

   for (int i=0; i<NOrbTot; i++){
      for (int j=0; j<i; j++)
         cout << "X\t";
      for (int j=i; j<NOrbTot; j++)
         cout << gTelem(i,j) << "\t";
      cout << endl;
   }

}


/**
 * Print the two body matrix
 */
void MxElem::PrintTwoBody(){

   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         for (int k=i; k<NOrbTot; k++)
            for (int l=j; l<NOrbTot; l++)
               cout << "Velem[" << i << "," << j << "," << k << "," << l << "] = " << gVelem(i,j,k,l) << endl;

}


/**
 * Print all matrix elements larger than a given cutoff
 * @param cutoff the cutoff above which the element should be printed
 */
void MxElem::PrintLarger(double cutoff){

   cout << "Overlap matrix elements" << endl;
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         if (fabs(gSoverlap(i,j))>cutoff)
            cout << "S[" << i << "," << j << "] = " << gSoverlap(i,j) << endl;
   cout << " " << endl;

   cout << "One body matrix elements" << endl;
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         if (fabs(gTelem(i,j))>cutoff)
            cout << "T[" << i << "," << j << "] = " << gTelem(i,j) << endl;
   cout << " " << endl;

   cout << "Two body matrix elements" << endl;
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         for (int k=i; k<NOrbTot; k++)
            for (int l=j; l<NOrbTot; l++)
               if (fabs(gVelem(i,j,k,l))>cutoff)
                  cout << "V[" << i << "," << j << "," << k << "," << l << "] = " << gVelem(i,j,k,l) << endl;
   cout << " " << endl;

}


/**
 * Print all matrix elements smaller than a given cutoff
 * @param cutoff the cutoff below which an element should be printed
 */
void MxElem::PrintSmaller(double cutoff){

   cout << "Overlap matrix elements" << endl;
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         if (fabs(gSoverlap(i,j))<cutoff)
            if (gSoverlap(i,j)!=0.0)
               cout << "S[" << i << "," << j << "] = " << gSoverlap(i,j) << endl;
   cout << " " << endl;

   cout << "One body matrix elements" << endl;
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         if (fabs(gTelem(i,j))<cutoff)
            if (gTelem(i,j)!=0.0)
               cout << "T[" << i << "," << j << "] = " << gTelem(i,j) << endl;
   cout << " " << endl;

   cout << "Two body matrix elements" << endl;
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         for (int k=i; k<NOrbTot; k++)
            for (int l=j; l<NOrbTot; l++)
               if (fabs(gVelem(i,j,k,l))<cutoff)
                  if (gVelem(i,j,k,l)!=0.0)
                     cout << "V[" << i << "," << j << "," << k << "," << l << "] = " << gVelem(i,j,k,l) << endl;
   cout << " " << endl;

}


/**
 * Set all overlap matrix elements smaller than a cutoff to zero
 * @param cutoff the cutoff for setting the mx elements to zero
 */
void MxElem::SetOverlapSmallToZero(double cutoff){

   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         if (fabs(gSoverlap(i,j))<cutoff)
            if (gSoverlap(i,j)!=0.0)
               setSoverlap(i,j,0.0);

}


/**
 * Set all one body matrix elements smaller than a cutoff to zero
 * @param cutoff the cutoff for setting the mx elements to zero
 */
void MxElem::SetTelemSmallToZero(double cutoff){

   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         if (fabs(gTelem(i,j))<cutoff)
            if (gTelem(i,j)!=0.0)
               setTelem(i,j,0.0);

}


/**
 * Set all two body matrix elements smaller than a cutoff to zero
 * @param cutoff the cutoff for setting the mx elements to zero
 */
void MxElem::SetVelemSmallToZero(double cutoff){

   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++)
         for (int k=i; k<NOrbTot; k++)
            for (int l=j; l<NOrbTot; l++)
               if (fabs(gVelem(i,j,k,l))<cutoff)
                  if (gVelem(i,j,k,l)!=0.0)
                     setVelem(i,j,k,l,0.0);

}


/**
 * Print the right upper triangle of a matrix A with linear dimension N
 * @param A pointer to the matrix to be printed
 * @param N the linear dimension of the matrix
 */
void MxElem::PrintUpper(double * A, int N){

   cout << "PrintUpper (dim " << N << ")" << endl;
   for (int i=0; i<N; i++){
      for (int j=0; j<i; j++)
         cout << "X\t";
      for (int j=i; j<N; j++)
         cout << A[i+N*j] << "\t";
      cout << endl;
   }

}


/**
 * Print a matrix A with linear dimension N
 * @param A pointer to the matrix to be printed
 * @param N the linear dimension of the matrix
 */
void MxElem::PrintFull(double * A, int N){

   cout << "PrintFull (dim " << N << ")" << endl;
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++)
         cout << A[i+N*j] << "\t";
      cout << endl;
   }

}


/**
 * Print an array A with linear dimension N
 * @param A pointer to the array to be printed
 * @param N the linear dimension of the matrix
 */
void MxElem::PrintArray(double * A, int N){

   cout << "PrintArray (dim " << N << ")" << endl;
   for (int i=0; i<N; i++)
      cout << A[i] << "\t";
   cout << endl;

}





