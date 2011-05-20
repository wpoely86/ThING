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
      
   File information Diag.h and Diag.cpp
   
      A class to compute the exact solution. More proof of concept &
      tester of the class MxElem than an efficient algorithm. Note
      that the overlap matrix is assumed to be the identity (as can
      be reached by e.g. a Lodwin transformation).
    
*/

#include<iostream>
#include<gmpxx.h>
#include "input.h"
#include "lapack.h"
#include "MxElem.h"
#include "Diag.h"

using std::cout;
using std::endl;

using namespace ThING;

/** 
 * Constructor for the Diag class
 */
Diag::Diag(){

}


/** 
 * Copy constructor for the Diag class
 */
Diag::Diag(Diag & tocopy){

}


/** 
 * Standard destructor
 */
Diag::~Diag(){

}


/**
 * Calculate the exact energy
 * @param problem definition of the problem (read in from a setup file)
 * @param elements the overlap, one body and two body matrix elements of the problem
 * @return the exact energy
 */
double Diag::CalcEnergy(input & problem, MxElem & elements){

   int Nelectrons = problem.NumberOfElectrons();
   int Norbitals = elements.gNOrbTot();
   mpz_class bin;
   unsigned long int upper = 2*Norbitals;
   unsigned long int lower = Nelectrons;
   mpz_bin_uiui(bin.get_mpz_t(),upper,lower);
   int Nstates = (int)bin.get_si();
   
   double * Hamiltonian = new double[Nstates*Nstates];
   for (int i=0; i<Nstates*Nstates; i++){
      Hamiltonian[i] = 0.0;
   }

   int * rowstate = new int[Nelectrons];
   int * colstate = new int[Nelectrons];
   for (int i=0; i<Nelectrons; i++){
      rowstate[i] = i;
      colstate[i] = i;
   }
   
   for (int colcnt = 0; colcnt<Nstates; colcnt++){
      if (colcnt!=0){
         Next(colstate, Norbitals, Nelectrons);
      }
      for (int rowcnt = 0; rowcnt<Nstates; rowcnt++){
         if (rowcnt!=0){
            Next(rowstate, Norbitals, Nelectrons);
         }
         // Calculation of the matrix element
         // 1. Check whether spin projection is preserved
         int sL = 0;
         int sR = 0;
         for (int i=0; i<Nelectrons; i++){
            sL += rowstate[i]%2;
            sR += colstate[i]%2;
         }
         if (sL==sR){
            // 2. Get the number of different indices and the first two we can find.
            // 2a. Find the left indices (row).
            int Ndifferent = 0;
            int i1, i2, j1, j2;
            
            for (int i=0; i<Nelectrons; i++){
               bool same = false;
               for (int j=0; j<Nelectrons; j++){
                  if (rowstate[i]==colstate[j]){
                     same = true;
                     j = Nelectrons;
                  }
               }
               if (!same){
                  Ndifferent++;
                  if (Ndifferent==1) i1 = i;
                  if (Ndifferent==2) i2 = i;
               }
            }
            
            // 2b. Find the right indices (col).
            Ndifferent = 0;
            for (int j=0; j<Nelectrons; j++){
               bool same = false;
               for (int i=0; i<Nelectrons; i++){
                  if (rowstate[i]==colstate[j]){
                     same = true;
                     i = Nelectrons;
                  }
               }
               if (!same){
                  Ndifferent++;
                  if (Ndifferent==1) j1 = j;
                  if (Ndifferent==2) j2 = j;
               }
            }
            
            // 3. One body terms
            if (Ndifferent == 1){
               Hamiltonian[rowcnt + Nstates*colcnt] += elements.gTelem(rowstate[i1]/2,colstate[j1]/2)*((i1-j1)%2==0?1:-1);
            }
            if (Ndifferent == 0){
               for (int i=0; i<Nelectrons; i++){
                  Hamiltonian[rowcnt + Nstates*colcnt] += elements.gTelem(rowstate[i]/2,rowstate[i]/2);
               }
            }
            
            // 4. Two body terms
            if (Ndifferent == 0){
               for (int k=0; k<Nelectrons-1; k++){
                  for (int l=k+1; l<Nelectrons; l++){
                     Hamiltonian[rowcnt+Nstates*colcnt] += elements.gVelem(rowstate[k]/2,rowstate[l]/2,rowstate[k]/2,rowstate[l]/2);
                     if ((rowstate[k]-rowstate[l])%2==0){
                        Hamiltonian[rowcnt+Nstates*colcnt] -= elements.gVelem(rowstate[k]/2,rowstate[l]/2,rowstate[l]/2,rowstate[k]/2);
                     }
                  }
               }
            }
            
            if (Ndifferent == 2){
               if ((rowstate[i1]-colstate[j1])%2==0){
                  int phase = ((i1-j1)%2==0?1:-1)*((i2-j2)%2==0?1:-1);
                  if (!(((j1>j2)&&(i1>i2))||((j1<j2)&&(i1<i2)))){
                     phase *= -1;
                  }
                  Hamiltonian[rowcnt+Nstates*colcnt] += phase*elements.gVelem(rowstate[i1]/2,rowstate[i2]/2,colstate[j1]/2,colstate[j2]/2);
               }
               
               if ((rowstate[i1]-colstate[j2])%2==0){
                  int phase = ((i1-j2)%2==0?1:-1)*((i2-j1)%2==0?1:-1);
                  if (((j1>j2)&&(i1>i2))||((j1<j2)&&(i1<i2))){
                     phase *= -1;
                  }
                  Hamiltonian[rowcnt+Nstates*colcnt] += phase*elements.gVelem(rowstate[i1]/2,rowstate[i2]/2,colstate[j2]/2,colstate[j1]/2);
               }
            }
            
            if (Ndifferent == 1){
            
               for (int i=0; i<Nelectrons; i++){
                  if (i!=i1){
                     Hamiltonian[rowcnt+Nstates*colcnt] += elements.gVelem(rowstate[i1]/2,rowstate[i]/2,colstate[j1]/2,rowstate[i]/2)*((i1-j1)%2==0?1:-1);
                     if ((rowstate[i]-rowstate[i1])%2==0){
                        Hamiltonian[rowcnt+Nstates*colcnt] += elements.gVelem(rowstate[i1]/2,rowstate[i]/2,rowstate[i]/2,colstate[j1]/2)*((i1-j1)%2==0?-1:1);
                     }
                  }
               }
            
            }

         }
      }
      for (int i=0; i<Nelectrons; i++){
         rowstate[i] = i;
      }
   }
   
   delete [] rowstate;
   delete [] colstate;
   
   char jobz = 'N';
   char uplo = 'U';
   int N = Nstates;
   int lda = N;
   double * W = new double[N];
   int info = 0;
   int lwork = 3*N-1;
   double * work = new double[lwork];
   
   dsyev_(&jobz, &uplo, &N, Hamiltonian, &lda, W, work, &lwork, &info);
   double energy = W[0];
   
   if (info!=0)
      cout << "Diag::CalcEnergy(...) The info parameter after diagonalisation is not zero." << endl;
   
   delete [] work;
   delete [] W;
   delete [] Hamiltonian;
   
   energy += elements.NuclPotEn(problem);
   
   cout << "The exact ground state energy  =  " << energy << endl;
   
   return energy;
   
}


/**
 * Give the next state in Nelectrons particle space.
 * @param state on entry the current state, on exit the next state
 * @param Norbitals the number of spatial orbitals of the problem
 * @param Nelectrons the number of electrons
 */
void Diag::Next(int * state, int Norbitals, int Nelectrons){

   int index = Nelectrons-1;
   bool stop = false;
   
   while(!stop){
   
      if (state[index]!=2*Norbitals+index-Nelectrons){
         stop = true;
         state[index]++;
         for (int j=index+1; j<Nelectrons; j++){
            state[j] = state[index] + (j-index);
         }
      } else {
         index--;
      }
      
      if (index<0){
         stop = true;
         cout << "Diag::Next() index smaller than zero." << endl;
      }
   
   }
   
}

 
