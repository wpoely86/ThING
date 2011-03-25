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
      
   File information main.cpp
   
      A file (in this case "start.stp") containing the chemical
      configuration is read in by the class input. It also reads in
      the required gaussian contraction parameters. Other examples
      can be found in the folder "examples/". A MxElem object is
      created and all the matrix elements are calculated by the
      function Init. The Hartree Fock energy is calculated.
      
   General information
   
      Notice that you should have installed LAPACK, BLAS and GMP.
      I have tested it with ATLAS and GMP on Ubuntu 10.10.
    
*/

#include "preamble.h"
#include "input.h"
#include "MxElem.h"
#include "HF.h"

void printGPL(){

   cout << "   Program information\n\n      ThING Is Not Gaussian is a rudimentary program to calculate \n      matrix elements of gaussian contractions using recursion\n      relations.\n\n   Copyright notice\n\n      Copyright 2011, 2011 Sebastian Wouters\n      <sebastianwouters@gmail.com>\n\n   Copyright permission\n\n      This file is part of ThING Is Not Gaussian.\n\n      ThING Is Not Gaussian is free software: you can redistribute\n      it and/or modify it under the terms of the GNU General \n      Public License as published by the Free Software Foundation,\n      either version 3 of the License, or (at your option) any\n      later version.\n\n      ThING Is Not Gaussian is distributed in the hope that it will\n      be useful, but WITHOUT ANY WARRANTY; without even the implied\n      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n      See the GNU General Public License for more details.\n\n      You should have received a copy of the GNU General Public\n      License along with ThING Is Not Gaussian. If not, see \n      <http://www.gnu.org/licenses/>.\n" << endl;

}

using namespace ThING;

int main(void){

   printGPL();
   
   cout.precision(12);

   input readin("start.stp");

   MxElem setup(readin);
   setup.Init(readin);
   
   setup.PrintS();
   setup.PrintOneBody();
   setup.PrintTwoBody();
   
   //setup.Save();
   //setup.Load();
   //setup.RemoveFromDisk();
   
   HF HFSolver;
   HFSolver.CalcEnergy(readin,setup);
   
   setup.DoLodwinTfo();
   
   setup.PrintS();
   setup.PrintOneBody();
   setup.PrintTwoBody();
   
   HFSolver.CalcEnergy(readin,setup);

   return 0;

}

