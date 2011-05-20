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

#ifndef HF_H
#define HF_H

#include "input.h"
#include "MxElem.h"

namespace ThING
{

class HF{

   public:

      //Constructor
      HF();

      //Copy constructor
      HF(HF &);

      //Destructor
      virtual ~HF();

      //Get the condition number of a real symmetric matrix
      double GetConditionNumber(double *, int);

      //Calculate the HF g.s. energy
      double CalcEnergy(input &, MxElem &);

   private:

      //Build the Fock operator
      void BuildFock(double *, double *, MxElem &, int, double);
      
      //Calculate the 2-norm of the difference of two arrays
      double TwoNormOfDifference(double *, double *, int);

};

}

#endif

