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
      
   File information input.h and input.cpp
   
      Handling class to do all the dirty read-in jobs: setup-file,
      elements mapper and the basisset.
      
      Comment 1 : It contains all the necessary info for the MxElem
                  class to do its job.
      Comment 2 : Unfortunately there are still fragments of bad
                  programming, but, as you all will check, they are
                  of no significance to the CPU time.
      Comment 3 : This is actually the first important one : The
                  basissets have to be in the Gaussian 94 format and
                  they can be downloaded from 
                  <https://bse.pnl.gov/bse/portal>.
      Comment 4 : The second important one : the setupfile must be
                  specifically formatted. Examples are given in the 
                  folder "examples/".
    
*/

#ifndef INPUT_H
#define INPUT_H

#include "preamble.h"
#include "Gauss.h"
#include "R.h"

namespace ThING
{

class input{

   public:

      //Constructor
      input(string);

      //Copy constructor
      input(input &);

      //Destructor
      virtual ~input();

      //Getters
      char gRotationSymm();
      int gCharge();
      int gNcores();
      string gbasisset();
      int gcore(int);
      R * gvector(int);
      Gauss * gGaussInfo(int);
      int NumberOfElectrons();

      //Printer
      void debugprint();

   private:

      //!To know the proton number of an element: elements[Z-1] = "name".
      string * elements;

      //!Charge = number of protons - number of electrons
      int Charge;
      
      //!RotationSymm = n (no), x, y or z rotation axis
      char RotationSymm;

      //!Basisset
      string basisset;

      //!Number of cores
      int Ncores;

      //!Cores (atomic number)
      int * cores;

      //!Positions of the cores
      R ** vectors;

      //!Per core: a Gauss object that stores the basic info to construct the desired Gaussian contractions
      Gauss ** GaussInfo;

      //Helper functions
      void initelements();
      void readinsetupfile(string);
      void fillgaussinfo();
      int getZ(string);

};

}

#endif

