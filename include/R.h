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
      
   File information R.h and R.cpp
   
      R is a container class for a vector. Basic operations are available.
    
*/

#ifndef R_H
#define R_H

#include "preamble.h"

class R{

   friend ostream &operator<<(ostream &output, R &);

   public:

      //Constructor
      R(double, double, double);

      //Copy constructor
      R(R &);

      //Destructor
      virtual ~R();

      //Getters
      double gxco();
      double gyco();
      double gzco();

      //Setter
      void set(double, double, double);

      //Distance
      double DistanceSquared(R &);

      //Product
      double product(R &);

   private:

      //!x-co
      double xco;

      //!y-co
      double yco;

      //!z-co
      double zco;
      
};

#endif

