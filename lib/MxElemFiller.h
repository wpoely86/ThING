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
      
   File information MxElemFiller.h and MxElemFiller.cpp
   
      All electron integral recursions as noted in
      "Efficient recursive computation of molecular integrals
      over Cartesian Gaussian Functions", by S. Obara and A. Saika,
      in J. Chem. Phys. 84 (7), 1986.

      The Fm(U) function recursion can be easily proven from its definition.

      Efficiency can be improved most in this class. But, as Gaussian sometimes states:
      
         THERE IS MORE TO LIFE THAN INCREASING ITS SPEED.
                                        -- GANDHI
    
*/

#ifndef MXELEMFILLER_H
#define MXELEMFILLER_H

#include "input.h"
#include "R.h"

namespace ThING
{

class MxElemFiller{

   public:

      //Constructor
      MxElemFiller(input &);

      //Copy constructor
      MxElemFiller(MxElemFiller &);

      //Destructor
      virtual ~MxElemFiller();

      //Getter
      input * ginfo();

      //Calculators
      double Overlap(int, int, int, int, int, int, int, int, int, int);
      double KE(int, int, int, int, int, int, int, int, int, int);
      double ElNucl(int, int, int, int, int, int, int, int, int, int, int);
      double ElEl(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int);

   private:    

      //!From here one can calculate everything recursive
      input * info;

      //Internal calculators
      double OverlapRecursion(double, double, R &, R &, int, int, int, int, int, int);
      double KERecursion(double, double, R &, R &, int, int, int, int, int, int);
      double ElNucRecursion(double, double, R &, R &, R &, int, int, int, int, int, int, int);
      double ElElRecursion(int, double, R &, int, int, int, double, R &, int, int, int, double, R &, int, int, int, double, R &, int, int, int);
      double FactHelper(int);
      double NormCst(double,int,int,int);
      double FunctionF(int, double);


};

}

#endif
