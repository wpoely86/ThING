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
      
   File information MxElem.h and MxElem.cpp
   
      Class to fetch & store the matrix elements of a chemical
      problem. It contains:
         - the overlap matrix elements
         - the one body matrix elements
         - the kinetic energy matrix elements (so that nuclear charge
           rescaling can be used in HF if desired)
         - the teo body matrix elements

      Basic operations concerning the overlap matrix (Lodwin tfo,
      canonical tfo, inverse overlap matrix) are available
      "to a certain extent".
    
*/

#ifndef MXELEM_H
#define MXELEM_H

#include "preamble.h"
#include "input.h"

namespace ThING
{

class MxElem{

   public:

      //Constructor
      MxElem(input &);

      //Copy constructor
      MxElem(MxElem &);

      //Destructor
      virtual ~MxElem();

      //Getters
      int gNOrbTot();
      double gKEseparate(int,int);
      double gTelem(int,int);
      double gSoverlap(int,int);
      double gVelem(int,int,int,int);
      double gVelemOK(int,int,int,int);

      //Binary storage: One should be aware to use exactly the same setup file of course.
      void Save();
      void Load();
      void RemoveFromDisk();

      //Init the gaussians
      void Init(input &);

      //Different orthogonalisations
      void CanOrth(double *);

      //Inverse overlap matrix
      void InverseOverlap(double *);
      
      //Construct the Lodwin Transformation
      void MakeLodwinTfo(double *);
      
      //Perform a Lodwin Transformation
      void DoLodwinTfo();

      //Printers
      void PrintUpper(double *, int);
      void PrintFull(double *, int);
      void PrintArray(double *, int);
      void PrintS();
      void PrintOneBody();
      void PrintTwoBody();
      void PrintLarger(double);
      void PrintSmaller(double);

      //Set certain matrix elements smaller than a cutoff to zero
      void SetOverlapSmallToZero(double);
      void SetVelemSmallToZero(double);
      void SetTelemSmallToZero(double);
      
      //Calculate the nuclear nuclear repulsion energy of the problem
      double NuclPotEn(input &);

   private:    

      //! Total number of orbitals
      int NOrbTot;

      //! 1body: Telem[index] = (i|T|j) with j>=i and index = i + j*(j+1)/2. Minimum storage for upper triangular filling.
      double * Telem;
      
      //! 1body: KEseparate[index] = (i|T|j) with j>=i and index = i + j*(j+1)/2. Minimum storage for upper triangular filling.
      double * KEseparate;

      //! 2body: Velem[i][x][y][z] = (ij|V|kl) with k>=i; j>=i; l>=j. [ i : 0->N-1 ][ x = (j-i) : 0->N-i ][ y = (k-i) : 0->N-i ][ z = l-j : 0->N-j : 0->N-i-x]
      double **** Velem;

      //! Overlap: S[index] = (i|j) with j>=i and index = i + j*(j+1)/2. Minimum storage for upper triangular filling.
      double * Soverlap;

      //Calculate total number of basis functions
      int CalcTotalNumberOfOrbitals(input &);

      //Get the orbital angular momentum corresponding to a subshell (char s,p,d,f...)
      int GetLofType(char);

      //Setters
      void setKEseparate(int,int,double);
      void setTelem(int,int,double);
      void setSoverlap(int,int,double);
      void setVelem(int,int,int,int,double);
      void setVelemOK(int,int,int,int,double);      

      //allocate
      void allocate();

};

}

#endif
