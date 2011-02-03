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
      
   File information HF.h and HF.cpp
   
      A class to compute the unrestricted Hartree Fock solution.
      More proof of concept & tester of the class MxElem than an
      efficient and always correctly converging algorithm.
    
*/

#include "preamble.h"
#include "input.h"
#include "MxElem.h"
#include "HF.h"


/** 
 * Constructor for the HF class
 */
HF::HF(){

}


/** 
 * Copy constructor for the HF class
 */
HF::HF(HF & tocopy){

}


/** 
 * Standard destructor
 */
HF::~HF(){

}


/** 
 * Calculate the energy of the unrestricted HF solution
 * @param problem definition of the problem (read in from a setup file)
 * @param elements the overlap, one body and two body matrix elements of the problem
 * @param generalized boolean to decide whether to solve the generalized eigenvalue problem or to transform the problem
 * @return the HF energy
 */
double HF::CalcEnergy(input & problem, MxElem & elements, bool generalized){

   //Init the problem

   int Nelectrons = problem.NumberOfElectrons();
   int NOrbTot = elements.gNOrbTot();
   double mixingfactor = 0.0; // mixingfactor*Fock[previous] + (1-mixingfactor)*Fock[now]
   cout.precision(12);

   //Init part 1: Declare some constants for the problem solver

   int itype = 1; // A x = lambda B x
   char jobz = 'V'; //eigenvalues and eigenvectors
   char uplo = 'U'; //upper part of symmetric matrices is stored
   int N = 2*NOrbTot; //linear length of the square matrices
   int lda = N;
   int ldb = N;
   int lwork = 1 + 6*N + 2*N*N;
   int liwork = 3 + 5*N;
   int info;

   if (N<Nelectrons){
      cout << "Not enough orbitals" << endl;
      assert(N>=Nelectrons);
   }

   //Init part 2: Declare the necessary arrays

   double * S = new double[N*N]; //S contains the overlap mx (mx "B" of the problem) on entry. On exit its values are destroyed.
   double * F = new double[N*N]; //F contains the Fock operator (mx "A" of the problem) on entry. On exit it contains the eigenvectors, normalized as V^TSV = I.
   double * F_old = new double[N*N];
   double * eigs = new double[N];
   double * work = new double[lwork];
   int * iwork = new int[liwork];

   //Init part 3 : Obtain the canonical transformation (S) for the first guess and possibly for the TFO mx.
   double * V = new double[NOrbTot*NOrbTot];
   //elements.CanOrth(V);

   for (int i=0; i<NOrbTot*NOrbTot; i++)
      V[i] = 0.0;

   //Init part 4 : Find out the condition number of the overlap matrix and if its too bad, reshape the problem.
   for (int i=0; i<NOrbTot; i++)
      for (int j=i; j<NOrbTot; j++){
         S[2*i + N*2*j] = S[2*i+1 + N*(2*j+1)] = elements.gSoverlap(i,j);
         S[2*i+1 + N*2*j] = S[2*i + N*(2*j+1)] = 0.0;
      }
   //cout << "Condition number S : " << GetConditionNumber(S,N) << endl;
   double * TFO;
   double * G;
   if (!generalized){

      TFO = new double[N*N]; // TFO = V * D^(-1/2) => Can Orth
      for (int i=0; i<NOrbTot; i++)
         for (int j=0; j<NOrbTot; j++){
            //first of a pair is always up, second down: so ordering: (0, up) (0, down) (1, up) ....
            TFO[2*i + N*2*j] = TFO[2*i+1 + N*(2*j+1)] = V[i+NOrbTot*j];
            TFO[2*i + N*(2*j+1)] = TFO[2*i+1 + N*2*j] = 0.0;
         }
      G = new double[N*N]; // G = TFO^T * F * TFO

   }

   //Init part 5 : Get a first guess -> Can ortho wave functions
   for (int i=0; i<NOrbTot; i++)
      for (int j=0; j<NOrbTot; j++){
         //first of a pair is always up, second down: so ordering: (0, up) (0, down) (1, up) ....
         F[2*i + N*2*j] = F[2*i+1 + N*(2*j+1)] = V[i+NOrbTot*j];
         F[2*i + N*(2*j+1)] = F[2*i+1 + N*2*j] = 0.0;
      }   

   delete [] V;

   //Init part 6 : Declare some constants for matrix multiplication

   char trans = 'T';
   char notrans = 'N';
   int m = N;
   int kNel = Nelectrons;
   int kbis = N;
   int n = N;
   double alpha = 1.0;
   double beta = 0.0;
   int lda_ = N;
   int ldb_ = N;
   int ldc_ = N;

   //Init part 7 : Make sum lowest occ (j) of V[m,j]*V[n,j] and store in S[m,n] --> SYMMETRIC MX
   
   cout <<  "Number of electrons in the system: " << Nelectrons << endl;
   
   dgemm_(&notrans, &trans, &m, &n, &kNel, &alpha, F, &lda_, F, &ldb_, &beta, S, &ldc_);
   
   //Iterate:

   double convcrit = 1e-10;
   double FockDifference = 1.0;
   double Enow = 0.0;
   double ENuclPot = elements.NuclPotEn(problem);
   int iteration = 0;
   double Zscaling = 1.0;

   while (FockDifference>=convcrit){

      iteration++;

      //Step 1: Build the Fock Operator (for F only upper half needs to be filled)

      BuildFock(F, S, elements, NOrbTot, Zscaling);
      
      if (iteration!=1) FockDifference = TwoNormOfDifference(F, F_old, N);

      //Step 2: Mix the Fock operator with the old one and make a copy to F_old
      for (int i=0; i<N; i++)
         for (int j=i; j<N; j++){

            if (iteration==1)
               F_old[i+N*j] = F[i+N*j];
            else
               F_old[i+N*j] = F[i+N*j] = mixingfactor*F_old[i+N*j] + (1-mixingfactor)*F[i+N*j];

         }

      //Step 3: Solve the problem:
      //If generalized problem: fill S and solve Fx = lambda Sx (x eigenvectors stored in F).
      //If not: build G = TFO^T F TFO, solve Gz = lambda z, calculate x = TFO z and store in F.
      if (generalized){
      
         //Step 3(OK) part1 : Refill the upper triangle half of S as the values of S are destroyed by the algorithm.
         for (int i=0; i<NOrbTot; i++)
            for (int j=i; j<NOrbTot; j++){
               S[2*i + N*2*j] = S[2*i+1 + N*(2*j+1)] = elements.gSoverlap(i,j);
               S[2*i+1 + N*2*j] = S[2*i + N*(2*j+1)] = 0.0;
            }

         //Step 3(OK) part2: Solve the generalized eigenvalue problem (for both eigenvalues and eigenvectors).
         //cout << "Condition number F : " << GetConditionNumber(F,N) << endl;
         dsygvd_(&itype, &jobz, &uplo, &N, F, &lda, S, &ldb, eigs, work, &lwork, iwork, &liwork, &info);

      } else {
         //Step 3(NOK) part1: Build G
         for (int i=0; i<N; i++)
            for (int j=i; j<N; j++)
               F[j+N*i] = F[i+N*j]; //Maybe use a symmetric * general matrix multiplier in lapack. Check...

         dgemm_(&notrans, &notrans, &m, &n, &kbis, &alpha, F, &lda_, TFO, &ldb_, &beta, S, &ldc_); //F * TFO -> S
         dgemm_(&trans, &notrans, &m, &n, &kbis, &alpha, TFO, &lda_, S, &ldb_, &beta, G, &ldc_); //TFO^T * S -> G

         //Step 3(NOK) part2: Solve Gz = lambda z
         //cout << "Condition number G : " << GetConditionNumber(G,N) << " [condition number F : " << GetConditionNumber(F,N) << "]." << endl;
         dsyevd_(&jobz, &uplo, &N, G, &lda_, eigs, work, &lwork, iwork, &liwork, &info);
         
         //cout << "Info = " << info << endl;

         //Step 3(NOK) part3: Calc x = TFO z and store in F.
         dgemm_(&notrans, &notrans, &m, &n, &kbis, &alpha, TFO, &lda_, G, &ldb_, &beta, F, &ldc_); //TFO * G -> F

      }

      //Step 4: Make sum lowest occ (j) of V[m,j]*V[n,j] and store in S[m,n]

      dgemm_(&notrans, &trans, &m, &n, &kNel, &alpha, F, &lda_, F, &ldb_, &beta, S, &ldc_);

      //Step 5: Enow to Eprev, Calculate the HF energy and store it at Enow.
      Enow = ENuclPot;
      for (int i=0; i<Nelectrons; i++)
         Enow += eigs[i];

      for (int k=0; k<NOrbTot; k++){
         for (int m=0; m<NOrbTot; m++){
            for (int l=0; l<NOrbTot; l++){
               for (int n=0; n<NOrbTot; n++){
                  Enow -= 0.5*elements.gVelem(k,l,m,n)*(S[2*k+N*2*m]*S[2*l+N*2*n]+S[2*k+1+N*(2*m+1)]*S[2*l+N*2*n]+S[2*k+N*2*m]*S[2*l+1+N*(2*n+1)]+S[2*k+1+N*(2*m+1)]*S[2*l+1+N*(2*n+1)]); //Direct correction
                  Enow += 0.5*elements.gVelem(k,l,m,n)*(S[2*k+N*2*n]*S[2*l+N*2*m]+S[2*k+1+N*2*n]*S[2*l+N*(2*m+1)]+S[2*k+N*(2*n+1)]*S[2*l+1+N*2*m]+S[2*k+1+N*(2*n+1)]*S[2*l+1+N*(2*m+1)]); //Exchange correction
               }
            }
         }
      }
      
      cout << "    Energy HF at iteration " << iteration << "  =  " << Enow << " and ||F - Fold||  =  " << FockDifference << endl;

      if ((FockDifference<convcrit)&&(Zscaling!=1.0)){
         cout << "The calculation with Zscaling = " << Zscaling << " has converged: Energy HF at iteration " << iteration << "  =  " << Enow << endl;
         Enow = 0.0;
         iteration = 0;
         FockDifference = 1.0;
         if (fabs(Zscaling-1.0)<=0.01)
            Zscaling = 1.0;
         else
            Zscaling = pow(Zscaling,0.75);
      }

   }
   
   cout << "The calculation has converged: Energy HF  =  " << Enow << endl;

   //Terminate part 1: Delete arrays

   delete [] S;
   delete [] F;
   delete [] F_old;
   delete [] eigs;
   delete [] work;
   delete [] iwork;
   if (!generalized){
      delete [] TFO;
      delete [] G;
   }

   //Terminate part 2: Return the HF energy
   return Enow;

}


/**
 * Function to calculate the 2-norm of the difference of two symmetric matrices, stored in upper triangular form
 * @param F1 pointer to the first matrix
 * @param F2 pointer to the second matrix
 * @param Nlin linear dimension of the matrices
 * @return the desired 2-norm
 */
double HF::TwoNormOfDifference(double * F1, double * F2, int Nlin){

   double norm2 = 0.0;
   double interm;
   for (int i=0; i<Nlin; i++){
      interm = F1[i+Nlin*i] - F2[i+Nlin*i];
      norm2 += interm*interm;
      for (int j=i+1; j<Nlin; j++){
         interm = F1[i+Nlin*j] - F2[i+Nlin*j];
         norm2 += 2.0*interm*interm;
      }
   }

   norm2 = sqrt(norm2);
   return norm2;

}


/**
 * Function to build the fock operator
 * @param F pointer to the array where the Fock operator should be stored
 * @param S pointer to the array where S[m,n] = sum j occupied V[m,j]V[n,j]
 * @param NOrbTot the total number of orbitals (spin doubling not included)
 * @param Zscaling scaling factor of the nuclear charges
 */
void HF::BuildFock(double * F, double * S, MxElem & elements, int NOrbTot, double Zscaling){

   int N = NOrbTot*2;

   //for F only upper half needs to be filled
   for (int l=0; l<NOrbTot; l++)
      for (int k=l; k<NOrbTot; k++){

         //One body elements are diagonal in the spin
         F[2*l + N*2*k] = F[2*l+1 + N*(2*k+1)] = elements.gKEseparate(l,k) + Zscaling*(elements.gTelem(l,k)-elements.gKEseparate(l,k));
         F[2*l + N*(2*k+1)] = F[2*l+1 + N*2*k] = 0.0;

         //Two body elements: Direct coulomb is diagonal in the spin
         for (int m=0; m<NOrbTot; m++){

            double addition = elements.gVelem(l,m,k,m)*(S[2*m + N*2*m] + S[2*m+1 + N*(2*m+1)]);
            F[2*l + N*2*k] += addition;
            F[2*l+1 + N*(2*k+1)] += addition;

            for (int n=m+1; n<NOrbTot; n++){

               addition = 2.0*elements.gVelem(l,m,k,n)*(S[2*m + N*2*n] + S[2*m+1 + N*(2*n+1)]);
               F[2*l + N*2*k] += addition;
               F[2*l+1 + N*(2*k+1)] += addition;

            }
         }
      }

   //Two body elements: Exchange not necessarily diagonal
   for (int l=0; l<N; l++){
      int spinl = l%2; //0 is up, 1 is down
      for (int k=l; k<N; k++){
         int spink = k%2; //0 is up, 1 is down
         for (int m=0; m<NOrbTot; m++) //use 2*m+spink
            for (int n=0; n<NOrbTot; n++) //use 2*n+spinl
               F[l+N*k] -= elements.gVelem(l/2,m,n,k/2)*S[2*m+spink + N*(2*n+spinl)];
      }
   }

}


/**
 * Function to get the condition number of a square matrix
 * @param A pointer to the matrix elements
 * @param N linear dimensions of the matrix
 */
double HF::GetConditionNumber(double * A, int N){

   double * B = new double[N*N];

   //copy A to B
   int dim = N*N;
   int incA = 1;
   int incB = 1;
   dcopy_(&dim, A, &incA, B, &incB);

   //DSYTRF computes the factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method.
   char uplo = 'U';
   int n = N;
   int lda = N;
   int * ipiv = new int[N];
   int lwork = 4*N;
   double * work = new double[lwork];
   int info = -1;

   dsytrf_(&uplo, &n, B, &lda, ipiv, work, &lwork, &info);

   //calculate the one-norm of B
   char norm = '1';
   double onenorm = dlansy_(&norm, &uplo, &n, B, &lda, work);

   //calculate the condition number
   int * iwork = new int[N];
   double rcond = -1.0;
 
   dsycon_(&uplo, &n, B, &lda, ipiv, &onenorm, &rcond, work, iwork, &info);

   //deallocate memory
   delete [] iwork;
   delete [] work;
   delete [] ipiv;
   delete [] B;

   return 1.0/rcond;

}


