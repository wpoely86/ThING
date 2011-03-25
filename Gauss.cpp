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
      
   File information Gauss.h and Gauss.cpp
   
      Gauss is a container class to store the information for the
      gaussian contractions of all the orbitals of a single atom.
    
*/

#include "preamble.h"
#include "Gauss.h"

using namespace ThING;

/**
 * Constructor for the Gauss class
 * @param Ntypes the number of different contractions for the atom under consideration
 */
Gauss::Gauss(int NumberOfTypes){

   Ntypes = NumberOfTypes;
   Ncontr = new int[Ntypes];
   type = new char[Ntypes];
   alpha = new double*[Ntypes];
   prefactors = new double*[Ntypes];
   init = false;

}


/** 
 * Copy constructor for the Gauss class
 * @param Gauss_c the Gauss object to be copied
 */
Gauss::Gauss(Gauss & Gauss_c){

   Ntypes = Gauss_c.gNtypes();
   Ncontr = new int[Ntypes];
   type = new char[Ntypes];
   alpha = new double*[Ntypes];
   prefactors = new double*[Ntypes];
   init = Gauss_c.ginit();

   if (init) {
      for (int cnt=0; cnt<Ntypes; cnt++){
         Ncontr[cnt] = Gauss_c.gNcontr(cnt);
         type[cnt] = Gauss_c.gtype(cnt);
         alpha[cnt] = new double[Ncontr[cnt]];
         prefactors[cnt] = new double[Ncontr[cnt]];
         for (int cnt2=0; cnt2<Ncontr[cnt]; cnt2++){
            alpha[cnt][cnt2] = Gauss_c.galpha(cnt,cnt2);
            prefactors[cnt][cnt2] = Gauss_c.gprefactors(cnt,cnt2);
         }
      }
   }

}


/**
 * Standard destructor to deallocate the memory. It assumes that if one array alpha[i]; prefactors[i] has been allocated, they all have.
 */
Gauss::~Gauss(){

   if (init) {
      for (int cnt=0; cnt<Ntypes; cnt++){
         delete [] alpha[cnt];
         delete [] prefactors[cnt];
      }
   }

   delete [] Ncontr;
   delete [] type;
   delete [] alpha;
   delete [] prefactors;

}


/**
 * Getter function for the variable Ntypes
 * @return Ntypes the number of different contractions for the atom under consideration
 */
int Gauss::gNtypes(){

   return Ntypes;

}


/**
 * Getter function for the variable init
 * @return init whether on or more of the alpha[i]; prefactors[i] arrays have been allocated
 */
bool Gauss::ginit(){

   return init;

}


/**
 * Getter function for the number of gaussians of the i-th contraction
 * @param i the number of the contraction
 * @return Ncontr[i] the number of gaussians of the i-th contraction
 */
int Gauss::gNcontr(int i){

   return Ncontr[i];

}


/**
 * Getter function for the type of orbital the i-th contraction approximates
 * @param i the number of the contraction
 * @return type[i] the type of orbital the i-th contraction approximates
 */
char Gauss::gtype(int i){

   return type[i];

}


/**
 * Getter function for the j-th exponent of the i-th contraction
 * @param i the number of the contraction
 * @param j the number of the exponent
 * @return alpha[i][j] the j-th exponent of the i-th contraction
 */
double Gauss::galpha(int i, int j){

   return alpha[i][j];

}


/**
 * Getter function for the j-th coefficient of the i-th contraction
 * @param i the number of the contraction
 * @param j the number of the coefficient
 * @return prefeactors[i][j] the j-th coefficient of the i-th contraction
 */
double Gauss::gprefactors(int i, int j){

   return prefactors[i][j];

}


/**
 * Setter function for all the information of the i-th contraction
 * @param i the number of the contraction
 * @param Ngaussian the number of gaussians in the i-th contraction
 * @param alphas the exponents for the i-th contraction
 * @param coeff the coefficients for the i-th contraction
 * @param orbitaltype the type of orbital the i-th contraction represents
 */
void Gauss::set(int i, int Ngaussian, double * alphas, double * coeff, char orbitaltype){

   init = true;

   Ncontr[i] = Ngaussian;
   alpha[i] = new double[Ngaussian];
   prefactors[i] = new double[Ngaussian];
   for (int cnt=0; cnt<Ngaussian; cnt++){
      alpha[i][cnt] = alphas[cnt];
      prefactors[i][cnt] = coeff[cnt];
   }
   type[i] = orbitaltype;

}

namespace ThING
{
/**
 * Ostream overloader
 * @param output the ostream to write the info to
 * @param krusty the Gauss of which the info needs to be written
 */
ostream &operator<<(ostream &output, Gauss & krusty){

   output << "# types = " << krusty.gNtypes() << endl;
   for (int cnt=0; cnt<krusty.gNtypes(); cnt++){
      output << "    nr. " << cnt+1 << " is of type " << krusty.gtype(cnt) << " and has " << krusty.gNcontr(cnt) << " contractions." << endl;
      for (int cnt2=0; cnt2<krusty.gNcontr(cnt); cnt2++){
         output << "        contr nr. " << cnt2+1 << " has exponent " << krusty.galpha(cnt,cnt2) << " and prefactor " << krusty.gprefactors(cnt,cnt2) << endl;
      }
   }
   return output;

}
}




