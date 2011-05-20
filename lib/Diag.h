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

#ifndef DIAG_H
#define DIAG_H

#include "input.h"
#include "MxElem.h"

namespace ThING
{

class Diag{

   public:

      //Constructor
      Diag();

      //Copy constructor
      Diag(Diag &);

      //Destructor
      virtual ~Diag();

      //Calculate the exact ground state energy
      double CalcEnergy(input &, MxElem &);

   private:
   
      //Give the next state
      void Next(int *, int, int);

};

}

#endif

