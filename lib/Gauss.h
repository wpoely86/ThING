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

#ifndef GAUSS_H
#define GAUSS_H

namespace ThING
{

class Gauss{

   friend std::ostream &operator<<(std::ostream &output, Gauss &);

   public:

      //Constructor
      Gauss(int);

      //Copy constructor
      Gauss(Gauss &);

      //Destructor
      virtual ~Gauss();

      //Getters
      int gNtypes();
      bool ginit();
      int gNcontr(int);
      char gtype(int);
      double galpha(int, int);
      double gprefactors(int, int);

      //Setters
      void set(int, int, double *, double *, char);

   private:

      //!Number of types of orbitals (px, py and pz are still considered as 1 at this moment)
      int Ntypes;

      //!Number of contractions for the orbital
      int * Ncontr;

      //!Orbital type of the gaussian contraction: s,p,d,f...
      char * type;

      //!Exponents of the gaussian contraction
      double ** alpha;

      //!Prefactors of the gaussian contraction
      double ** prefactors;

      //!Whether one or more of the alpha[iii] arrays and prefactors[iii] arrays have been allocated
      bool init;

};

}

#endif

