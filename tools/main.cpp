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

#include <iostream>
#include <getopt.h>
#include <input.h>
#include <MxElem.h>
#include <HF.h>
#include <Diag.h>
#include "../config.h"

using namespace std;

using namespace ThING;

int main(int argc,char **argv){

   struct option long_options[] =
   {
      {"input",  required_argument, 0, 'i'},
      {"hfock",  no_argument, 0, 'f'},
      {"diag", no_argument, 0, 'd'},
      {"version",  no_argument, 0, 'v'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   bool hartreefock = false;
   bool diagonalisation = false;
   string filename = "";

   int i,j;
   while( (j = getopt_long (argc, argv, "hvi:fd", long_options, &i)) != -1)
       switch(j)
       {
	   case 'h':
	   case '?':
	       cout << "Usage: " << PACKAGE_NAME << " [OPTIONS]\n"
		   "ThING can calculate gaussian type matrix elements and use them\n"
		   "to do a Hartree-Fock or exact diagonalisation calculation.\n"
		   "\n"
		   "    -i, --input=file    Give a input file\n"
		   "    -f, --hfock         Use Hartree-Fock\n"
		   "    -d, --diag          Use exact diagonalisation\n"
		   "    -h, --help          Display this help\n"
		   "    -v, --version       Display the program version\n"
		   "\n"
		   "Either the -f or the -d option should be given.\n"
		   "\n"
		   << PACKAGE_STRING << "\n"
		   "Report bugs to " << PACKAGE_BUGREPORT << "\n";
	       return 0;
	       break;
	   case 'v':
	       cout << "   Program information\n\n"
		       "     ThING Is Not Gaussian is a rudimentary program to calculate\n"
		       "     matrix elements of gaussian contractions using recursion\n"
		       "     relations.\n"
		       "\n"
		       "   Copyright notice\n"
		       "\n"
		       "     Copyright 2011, 2011 Sebastian Wouters\n"
		       "     <sebastianwouters@gmail.com>\n"
		       "\n"
		       "   Copyright permission\n"
		       "\n"
		       "     This file is part of ThING Is Not Gaussian.\n"
		       "\n"
		       "     ThING Is Not Gaussian is free software: you can redistribute\n"
		       "     it and/or modify it under the terms of the GNU General\n"
		       "     Public License as published by the Free Software Foundation,\n"
		       "     either version 3 of the License, or (at your option) any\n"
		       "     later version.\n"
		       "\n"
		       "     ThING Is Not Gaussian is distributed in the hope that it will\n"
		       "     be useful, but WITHOUT ANY WARRANTY; without even the implied\n"
		       "     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
		       "     See the GNU General Public License for more details.\n"
		       "\n"
		       "     You should have received a copy of the GNU General Public\n"
		       "     License along with ThING Is Not Gaussian. If not, see \n"
		       "     <http://www.gnu.org/licenses/>.\n\n";
		   return 0;
	       break;
	   case 'i':
	       filename = optarg;
	       break;
	   case 'f':
	       hartreefock = true;
	       break;
	   case 'd':
	       diagonalisation = true;
	       break;
       }

   if( filename.empty() )
   {
       cout << "I need an inputfile to do any calculation. The format of the\n"
	       "inputfiles can be found in the examples directory" << endl;
       return 0;
   }

   if( !hartreefock && !diagonalisation )
   {
       cout << "You need to tell me what you want: Hartree-Fock or exact diagonalisation." << endl;
       cout << "Use the --help option to get more information." << endl;
       return 0;
   }

   cout.precision(12);

   input readin(filename);
   MxElem setup(readin);
   setup.Init(readin);

   if( hartreefock )
   {
       HF HFSolver;
       HFSolver.CalcEnergy(readin,setup);
   }

   if( diagonalisation )
   {
       setup.DoLodwinTfo();

       Diag DiagSolver;
       DiagSolver.CalcEnergy(readin,setup);
   }

   return 0;

}

