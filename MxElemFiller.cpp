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

#include "preamble.h"
#include "input.h"
#include "Gauss.h"
#include "R.h"
#include "MxElemFiller.h"


/**
 * Constructor for the MxElemFiller class
 * @param readin the problem to be solved
 */
MxElemFiller::MxElemFiller(input & readin){

   info = new input(readin);

}


/** 
 * Copy constructor for the MxElemFiller class
 * @param tocopy the MxElemFiller object to be copied
 */
MxElemFiller::MxElemFiller(MxElemFiller & tocopy){

   info = new input(*tocopy.ginfo());

}


/**
 * Standard destructor to deallocate the memory.
 */
MxElemFiller::~MxElemFiller(){

   delete info;

}


/**
 * Getter function for the problem info
 * @return info the problem info
 */
input * MxElemFiller::ginfo(){

   return info;

}


/**
 * Calculate the overlap matrix element for
 * @param atom_1 the first atom
 * @param orbit_1 the orbital of the first atom
 * @param n1x the power of (x-R1x) in front of the first Gaussian basis function
 * @param n1y the power of (y-R1y) in front of the first Gaussian basis function
 * @param n1z the power of (z-R1z) in front of the first Gaussian basis function
 * @param atom_2 the second atom
 * @param orbit_2 the orbital of the second atom
 * @param n2x the power of (x-R2x) in front of the second Gaussian basis function
 * @param n2y the power of (y-R2y) in front of the second Gaussian basis function
 * @param n2z the power of (z-R2z) in front of the second Gaussian basis function
 * @return the overlap matrix element
 */
double MxElemFiller::Overlap(int atom_1, int orbit_1, int n1x, int n1y, int n1z, int atom_2, int orbit_2, int n2x, int n2y, int n2z){

   R R1(*(info->gvector(atom_1)));
   R R2(*(info->gvector(atom_2)));

   Gauss * A1 = info->gGaussInfo(atom_1);
   Gauss * A2 = info->gGaussInfo(atom_2);

   int Ncontr1 = A1->gNcontr(orbit_1);
   int Ncontr2 = A2->gNcontr(orbit_2);

   double theOverlap = 0.0;

   double alpha1, alpha2, intermediate;

   for (int i=0; i<Ncontr1; i++)
      for (int j=0; j<Ncontr2; j++){

         alpha1 = A1->galpha(orbit_1,i);
         alpha2 = A2->galpha(orbit_2,j);
         intermediate = OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z);
         theOverlap += (A1->gprefactors(orbit_1,i))*(A2->gprefactors(orbit_2,j))*NormCst(alpha1,n1x,n1y,n1z)*NormCst(alpha2,n2x,n2y,n2z)*intermediate;

      }

   return theOverlap;

}


/**
 * Recursion relation for the calculation of the overlap matrix element
 * @param alpha1 first gaussian power parameter
 * @param alpha2 second gaussian power parameter
 * @param R1 the center of the first gaussian
 * @param R2 the center of the second gaussian
 * @param n1x the power of (x-R1x) in front of the first Gaussian basis function
 * @param n1y the power of (y-R1y) in front of the first Gaussian basis function
 * @param n1z the power of (z-R1z) in front of the first Gaussian basis function
 * @param n2x the power of (x-R2x) in front of the second Gaussian basis function
 * @param n2y the power of (y-R2y) in front of the second Gaussian basis function
 * @param n2z the power of (z-R2z) in front of the second Gaussian basis function
 * @return the overlap matrix element of a step in the recursion
 */
double MxElemFiller::OverlapRecursion(double alpha1, double alpha2, R & R1, R & R2, int n1x, int n1y, int n1z, int n2x, int n2y, int n2z){

   if ((n1x<0)||(n1y<0)||(n1z<0)||(n2x<0)||(n2y<0)||(n2z<0))
      return 0.0;

   double zeta = alpha1+alpha2;

   if (n1x==0){

      if (n2x==0){

         if (n1y==0){

            if (n2y==0){

               if (n1z==0){

                  if (n2z==0){
                     //all zero
                     return pow(M_PI/zeta,1.5) * exp(- alpha1*alpha2/(alpha1+alpha2)*R1.DistanceSquared(R2));
                  }

                  //all but n2z zero
                  return ((alpha1*(R1.gzco()) + alpha2*(R2.gzco()))/(alpha1+alpha2) - R2.gzco())*OverlapRecursion(alpha1, alpha2, R1, R2, n1x, n1y, n1z, n2x, n2y, n2z-1) + 0.5/zeta*(n2z-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z-2);
               }

               //n1x==n2x==n1y==n2y==0, but n1z not zero or smaller
               return ((alpha1*(R1.gzco()) + alpha2*(R2.gzco()))/(alpha1+alpha2) - R1.gzco())*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z-1,n2x,n2y,n2z) + 0.5/zeta*((n1z-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z-2,n2x,n2y,n2z) + n2z*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z-1,n2x,n2y,n2z-1));

            }

            //n1x==n2x==n1y==0, but n2y not zero or smaller
            return ((alpha1*(R1.gyco()) + alpha2*(R2.gyco()))/(alpha1+alpha2) - R2.gyco())*OverlapRecursion(alpha1, alpha2, R1, R2, n1x, n1y, n1z, n2x, n2y-1, n2z) + 0.5/zeta*(n2y-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y-2,n2z);

         }

         //n1x==n2x==0, but n1y not zero or smaller
         return ((alpha1*(R1.gyco()) + alpha2*(R2.gyco()))/(alpha1+alpha2) - R1.gyco())*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y-1,n1z,n2x,n2y,n2z) + 0.5/zeta*((n1y-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y-2,n1z,n2x,n2y,n2z) + n2y*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y-1,n1z,n2x,n2y-1,n2z)); 

      }

      //n1x==0, but n2x not zero or smaller
      return ((alpha1*(R1.gxco()) + alpha2*(R2.gxco()))/(alpha1+alpha2) - R2.gxco())*OverlapRecursion(alpha1, alpha2, R1, R2, n1x, n1y, n1z, n2x-1, n2y, n2z) + 0.5/zeta*(n2x-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x-2,n2y,n2z);

   }

   //n1x not zero or smaller
   return ((alpha1*(R1.gxco()) + alpha2*(R2.gxco()))/(alpha1+alpha2) - R1.gxco())*OverlapRecursion(alpha1,alpha2,R1,R2,n1x-1,n1y,n1z,n2x,n2y,n2z) + 0.5/zeta*((n1x-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x-2,n1y,n1z,n2x,n2y,n2z) + n2x*OverlapRecursion(alpha1,alpha2,R1,R2,n1x-1,n1y,n1z,n2x-1,n2y,n2z)); 

}


/**
 * Calculate the kinetic energy matrix element for
 * @param atom_1 the first atom
 * @param orbit_1 the orbital of the first atom
 * @param n1x the power of (x-R1x) in front of the first Gaussian basis function
 * @param n1y the power of (y-R1y) in front of the first Gaussian basis function
 * @param n1z the power of (z-R1z) in front of the first Gaussian basis function
 * @param atom_2 the second atom
 * @param orbit_2 the orbital of the second atom
 * @param n2x the power of (x-R2x) in front of the second Gaussian basis function
 * @param n2y the power of (y-R2y) in front of the second Gaussian basis function
 * @param n2z the power of (z-R2z) in front of the second Gaussian basis function
 * @return the kinetic energy matrix element
 */
double MxElemFiller::KE(int atom_1, int orbit_1, int n1x, int n1y, int n1z, int atom_2, int orbit_2, int n2x, int n2y, int n2z){

   R R1(*(info->gvector(atom_1)));
   R R2(*(info->gvector(atom_2)));

   Gauss * A1 = info->gGaussInfo(atom_1);
   Gauss * A2 = info->gGaussInfo(atom_2);

   int Ncontr1 = A1->gNcontr(orbit_1);
   int Ncontr2 = A2->gNcontr(orbit_2);

   double alpha1, alpha2, intermediate;

   double theKineticEnergy = 0.0;

   for (int i=0; i<Ncontr1; i++)
      for (int j=0; j<Ncontr2; j++){

         alpha1 = A1->galpha(orbit_1,i);
         alpha2 = A2->galpha(orbit_2,j);
         intermediate = KERecursion(alpha1, alpha2, R1, R2, n1x, n1y, n1z, n2x, n2y, n2z);
         theKineticEnergy += (A1->gprefactors(orbit_1,i))*(A2->gprefactors(orbit_2,j))*NormCst(alpha1,n1x,n1y,n1z)*NormCst(alpha2,n2x,n2y,n2z)*intermediate;

      }

   return theKineticEnergy;

}


/**
 * Recursion relation for the calculation of the kinetic energy matrix element
 * @param alpha1 first gaussian power parameter
 * @param alpha2 second gaussian power parameter
 * @param R1 the center of the first gaussian
 * @param R2 the center of the second gaussian
 * @param n1x the power of (x-R1x) in front of the first Gaussian basis function
 * @param n1y the power of (y-R1y) in front of the first Gaussian basis function
 * @param n1z the power of (z-R1z) in front of the first Gaussian basis function
 * @param n2x the power of (x-R2x) in front of the second Gaussian basis function
 * @param n2y the power of (y-R2y) in front of the second Gaussian basis function
 * @param n2z the power of (z-R2z) in front of the second Gaussian basis function
 * @return the kinetic energy matrix element of a step in the recursion
 */
double MxElemFiller::KERecursion(double alpha1, double alpha2, R & R1, R & R2, int n1x, int n1y, int n1z, int n2x, int n2y, int n2z){

   if ((n1x<0)||(n1y<0)||(n1z<0)||(n2x<0)||(n2y<0)||(n2z<0))
      return 0.0;

   double zeta = alpha1+alpha2;
   double xi = alpha1*alpha2/(alpha1+alpha2);

   if (n1x==0){

      if (n2x==0){

         if (n1y==0){

            if (n2y==0){

               if (n1z==0){

                  if (n2z==0){
                     //all zero
                     double distSq = R1.DistanceSquared(R2);
                     return xi * (3.0 - 2.0*xi*distSq) * pow(M_PI/zeta,1.5) * exp(- xi*distSq);
                  }

                  //all but n2z zero
                  return ((alpha1*(R1.gzco()) + alpha2*(R2.gzco()))/(alpha1+alpha2) - R2.gzco())*KERecursion(alpha1, alpha2, R1, R2, n1x, n1y, n1z, n2x, n2y, n2z-1) + 0.5/zeta*(n2z-1)*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z-2) + 2*xi*(OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z) - 0.5/alpha2*(n2z-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z-2));
               }

               //n1x==n2x==n1y==n2y==0, but n1z not zero or smaller
               return ((alpha1*(R1.gzco()) + alpha2*(R2.gzco()))/(alpha1+alpha2) - R1.gzco())*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z-1,n2x,n2y,n2z) + 0.5/zeta*((n1z-1)*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z-2,n2x,n2y,n2z) + n2z*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z-1,n2x,n2y,n2z-1)) + 2*xi*(OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z) - 0.5/alpha1*(n1z-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z-2,n2x,n2y,n2z));

            }

            //n1x==n2x==n1y==0, but n2y not zero or smaller
            return ((alpha1*(R1.gyco()) + alpha2*(R2.gyco()))/(alpha1+alpha2) - R2.gyco())*KERecursion(alpha1, alpha2, R1, R2, n1x, n1y, n1z, n2x, n2y-1, n2z) + 0.5/zeta*(n2y-1)*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y-2,n2z) + 2*xi*(OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z) - 0.5/alpha2*(n2y-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y-2,n2z));

         }

         //n1x==n2x==0, but n1y not zero or smaller
         return ((alpha1*(R1.gyco()) + alpha2*(R2.gyco()))/(alpha1+alpha2) - R1.gyco())*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y-1,n1z,n2x,n2y,n2z) + 0.5/zeta*((n1y-1)*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y-2,n1z,n2x,n2y,n2z) + n2y*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y-1,n1z,n2x,n2y-1,n2z)) + 2*xi*(OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z) - 0.5/alpha1*(n1y-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y-2,n1z,n2x,n2y,n2z));

      }

      //n1x==0, but n2x not zero or smaller
      return ((alpha1*(R1.gxco()) + alpha2*(R2.gxco()))/(alpha1+alpha2) - R2.gxco())*KERecursion(alpha1, alpha2, R1, R2, n1x, n1y, n1z, n2x-1, n2y, n2z) + 0.5/zeta*(n2x-1)*KERecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x-2,n2y,n2z) + 2*xi*(OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z) - 0.5/alpha2*(n2x-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x-2,n2y,n2z));

   }

   //n1x not zero or smaller
   return ((alpha1*(R1.gxco()) + alpha2*(R2.gxco()))/(alpha1+alpha2) - R1.gxco())*KERecursion(alpha1,alpha2,R1,R2,n1x-1,n1y,n1z,n2x,n2y,n2z) + 0.5/zeta*((n1x-1)*KERecursion(alpha1,alpha2,R1,R2,n1x-2,n1y,n1z,n2x,n2y,n2z) + n2x*KERecursion(alpha1,alpha2,R1,R2,n1x-1,n1y,n1z,n2x-1,n2y,n2z)) + 2*xi*(OverlapRecursion(alpha1,alpha2,R1,R2,n1x,n1y,n1z,n2x,n2y,n2z) - 0.5/alpha1*(n1x-1)*OverlapRecursion(alpha1,alpha2,R1,R2,n1x-2,n1y,n1z,n2x,n2y,n2z));

}


/**
 * Calculate the electron-nuclear interaction matrix element for
 * @param atom_1 the first atom
 * @param orbit_1 the orbital of the first atom
 * @param n1x the power of (x-R1x) in front of the first Gaussian basis function
 * @param n1y the power of (y-R1y) in front of the first Gaussian basis function
 * @param n1z the power of (z-R1z) in front of the first Gaussian basis function
 * @param atom_2 the second atom
 * @param orbit_2 the orbital of the second atom
 * @param n2x the power of (x-R2x) in front of the second Gaussian basis function
 * @param n2y the power of (y-R2y) in front of the second Gaussian basis function
 * @param n2z the power of (z-R2z) in front of the second Gaussian basis function
 * @param atom_3 center of the nucleus
 * @return the electron-nuclear interaction matrix element
 */
double MxElemFiller::ElNucl(int atom_1, int orbit_1, int n1x, int n1y, int n1z, int atom_2, int orbit_2, int n2x, int n2y, int n2z, int atom_3){

   R RA(*(info->gvector(atom_1)));
   R RB(*(info->gvector(atom_2)));
   R RC(*(info->gvector(atom_3)));

   Gauss * A1 = info->gGaussInfo(atom_1);
   Gauss * A2 = info->gGaussInfo(atom_2);

   int Z3 = info->gcore(atom_3);

   int Ncontr1 = A1->gNcontr(orbit_1);
   int Ncontr2 = A2->gNcontr(orbit_2);

   double alpha1, alpha2, intermediate;

   double theElNuclEn = 0.0;

   for (int i=0; i<Ncontr1; i++)
      for (int j=0; j<Ncontr2; j++){

         alpha1 = A1->galpha(orbit_1,i);
         alpha2 = A2->galpha(orbit_2,j);
         intermediate = ElNucRecursion(alpha1, alpha2, RA, RB, RC, 0, n1x, n1y, n1z, n2x, n2y, n2z);
         theElNuclEn += (A1->gprefactors(orbit_1,i))*(A2->gprefactors(orbit_2,j))*NormCst(alpha1,n1x,n1y,n1z)*NormCst(alpha2,n2x,n2y,n2z)*intermediate;

      }

   return -Z3*theElNuclEn;

}


/**
 * Recursion relation for the calculation of the electron-nuclear interaction matrix element
 * @param alpha1 first gaussian power parameter
 * @param alpha2 second gaussian power parameter
 * @param RA the center of the first gaussian
 * @param RB the center of the second gaussian
 * @param RC the location of the nucleus
 * @param m parameter in the recursion relation (see article)
 * @param n1x the power of (x-R1x) in front of the first Gaussian basis function
 * @param n1y the power of (y-R1y) in front of the first Gaussian basis function
 * @param n1z the power of (z-R1z) in front of the first Gaussian basis function
 * @param n2x the power of (x-R2x) in front of the second Gaussian basis function
 * @param n2y the power of (y-R2y) in front of the second Gaussian basis function
 * @param n2z the power of (z-R2z) in front of the second Gaussian basis function
 * @return the electron-nuclear interaction matrix element of a step in the recursion
 */
double MxElemFiller::ElNucRecursion(double alpha1, double alpha2, R & RA, R & RB, R & RC, int m, int n1x, int n1y, int n1z, int n2x, int n2y, int n2z){

   if ((n1x<0)||(n1y<0)||(n1z<0)||(n2x<0)||(n2y<0)||(n2z<0))
      return 0.0;

   double zeta = alpha1+alpha2;
   R P((alpha1*RA.gxco()+alpha2*RB.gxco())/zeta,(alpha1*RA.gyco()+alpha2*RB.gyco())/zeta,(alpha1*RA.gzco()+alpha2*RB.gzco())/zeta);

   if (n1x==0){

      if (n2x==0){

         if (n1y==0){

            if (n2y==0){

               if (n1z==0){

                  if (n2z==0){
                     //cout << "F(" << m << "," << zeta*P.DistanceSquared(RC) << ")  =  " << FunctionF(m,zeta*P.DistanceSquared(RC)) << endl;
                     //all zero
                     return 2.0 * FunctionF(m,zeta*P.DistanceSquared(RC)) * M_PI/zeta * exp(- alpha1*alpha2/zeta*RA.DistanceSquared(RB));
                  }

                  //all but n2z zero
                  return (P.gzco() - RB.gzco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z,n2x,n2y,n2z-1) - (P.gzco() - RC.gzco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z,n2x,n2y,n2z-1) + 0.5/zeta*(n2z-1)*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z,n2x,n2y,n2z-2) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z,n2x,n2y,n2z-2));
               }

               //n1x==n2x==n1y==n2y==0, but n1z not zero or smaller
               return (P.gzco() - RA.gzco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z-1,n2x,n2y,n2z) - (P.gzco() - RC.gzco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z-1,n2x,n2y,n2z) + 0.5/zeta*(n1z-1)*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z-2,n2x,n2y,n2z) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z-2,n2x,n2y,n2z)) + 0.5/zeta*n2z*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z-1,n2x,n2y,n2z-1) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z-1,n2x,n2y,n2z-1));

            }

            //n1x==n2x==n1y==0, but n2y not zero or smaller
            return (P.gyco() - RB.gyco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z,n2x,n2y-1,n2z) - (P.gyco() - RC.gyco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z,n2x,n2y-1,n2z) + 0.5/zeta*(n2y-1)*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z,n2x,n2y-2,n2z) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z,n2x,n2y-2,n2z));

         }

         //n1x==n2x==0, but n1y not zero or smaller
         return (P.gyco() - RA.gyco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y-1,n1z,n2x,n2y,n2z) - (P.gyco() - RC.gyco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y-1,n1z,n2x,n2y,n2z) + 0.5/zeta*(n1y-1)*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y-2,n1z,n2x,n2y,n2z) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y-2,n1z,n2x,n2y,n2z)) + 0.5/zeta*n2y*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y-1,n1z,n2x,n2y-1,n2z) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y-1,n1z,n2x,n2y-1,n2z));

      }

      //n1x==0, but n2x not zero or smaller
      return (P.gxco() - RB.gxco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z,n2x-1,n2y,n2z) - (P.gxco() - RC.gxco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z,n2x-1,n2y,n2z) + 0.5/zeta*(n2x-1)*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x,n1y,n1z,n2x-2,n2y,n2z) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x,n1y,n1z,n2x-2,n2y,n2z));

   }

   //n1x not zero or smaller
   return (P.gxco() - RA.gxco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x-1,n1y,n1z,n2x,n2y,n2z) - (P.gxco() - RC.gxco())*ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x-1,n1y,n1z,n2x,n2y,n2z) + 0.5/zeta*(n1x-1)*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x-2,n1y,n1z,n2x,n2y,n2z) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x-2,n1y,n1z,n2x,n2y,n2z)) + 0.5/zeta*n2x*(ElNucRecursion(alpha1,alpha2,RA,RB,RC,m,n1x-1,n1y,n1z,n2x-1,n2y,n2z) - ElNucRecursion(alpha1,alpha2,RA,RB,RC,m+1,n1x-1,n1y,n1z,n2x-1,n2y,n2z));

}


/**
 * Calculate the electron-electron interaction matrix element for
 * @param atom_1 the first atom
 * @param orbit_1 the orbital of the first atom
 * @param n1x the power of (x-R1x) in front of the first Gaussian basis function
 * @param n1y the power of (y-R1y) in front of the first Gaussian basis function
 * @param n1z the power of (z-R1z) in front of the first Gaussian basis function
 * @param atom_2 the second atom
 * @param orbit_2 the orbital of the second atom
 * @param n2x the power of (x-R2x) in front of the second Gaussian basis function
 * @param n2y the power of (y-R2y) in front of the second Gaussian basis function
 * @param n2z the power of (z-R2z) in front of the second Gaussian basis function
 * @param atom_3 the third atom
 * @param orbit_3 the orbital of the third atom
 * @param n3x the power of (x-R3x) in front of the third Gaussian basis function
 * @param n3y the power of (y-R3y) in front of the third Gaussian basis function
 * @param n3z the power of (z-R3z) in front of the third Gaussian basis function
 * @param atom_4 the fourth atom
 * @param orbit_4 the orbital of the fourth atom
 * @param n4x the power of (x-R4x) in front of the fourth Gaussian basis function
 * @param n4y the power of (y-R4y) in front of the fourth Gaussian basis function
 * @param n4z the power of (z-R4z) in front of the fourth Gaussian basis function
 * @return the electron-electron interaction matrix element
 */
double MxElemFiller::ElEl(int atom_1, int orbit_1, int n1x, int n1y, int n1z, int atom_2, int orbit_2, int n2x, int n2y, int n2z, int atom_3, int orbit_3, int n3x, int n3y, int n3z, int atom_4, int orbit_4, int n4x, int n4y, int n4z){

   R RA(*(info->gvector(atom_1)));
   R RB(*(info->gvector(atom_2)));
   R RC(*(info->gvector(atom_3)));
   R RD(*(info->gvector(atom_4)));

   Gauss * Ai = info->gGaussInfo(atom_1);
   Gauss * Bi = info->gGaussInfo(atom_2);
   Gauss * Ci = info->gGaussInfo(atom_3);
   Gauss * Di = info->gGaussInfo(atom_4);

   int NcontrA = Ai->gNcontr(orbit_1);
   int NcontrB = Bi->gNcontr(orbit_2);
   int NcontrC = Ci->gNcontr(orbit_3);
   int NcontrD = Di->gNcontr(orbit_4);

   double alphaA, alphaB, alphaC, alphaD, intermediate;

   double theElElEn = 0.0;

   for (int i=0; i<NcontrA; i++)
      for (int j=0; j<NcontrB; j++)
         for (int k=0; k<NcontrC; k++)
            for (int l=0; l<NcontrD; l++){

               alphaA = Ai->galpha(orbit_1,i);
               alphaB = Bi->galpha(orbit_2,j);
               alphaC = Ci->galpha(orbit_3,k);
               alphaD = Di->galpha(orbit_4,l);

               intermediate = ElElRecursion(0, alphaA, RA, n1x, n1y, n1z, alphaB, RB, n2x, n2y, n2z, alphaC, RC, n3x, n3y, n3z, alphaD, RD, n4x, n4y, n4z);

               theElElEn += (Ai->gprefactors(orbit_1,i))*(Bi->gprefactors(orbit_2,j))*(Ci->gprefactors(orbit_3,k))*(Di->gprefactors(orbit_4,l))*NormCst(alphaA,n1x,n1y,n1z)*NormCst(alphaB,n2x,n2y,n2z)*NormCst(alphaC,n3x,n3y,n3z)*NormCst(alphaD,n4x,n4y,n4z)*intermediate;

            }      

   return theElElEn;

}


/**
 * Recursion relation for the calculation of the electron-electron interaction matrix element
 * @param m parameter in the recursion relation (see article)
 * @param alphaA first gaussian power parameter
 * @param RA the center of the first gaussian
 * @param nAx the power of (x-RAx) in front of the first Gaussian basis function
 * @param nAy the power of (y-RAy) in front of the first Gaussian basis function
 * @param nAz the power of (z-RAz) in front of the first Gaussian basis function
 * @param alphaB second gaussian power parameter
 * @param RB the center of the second gaussian
 * @param nBx the power of (x-RBx) in front of the second Gaussian basis function
 * @param nBy the power of (y-RBy) in front of the second Gaussian basis function
 * @param nBz the power of (z-RBz) in front of the second Gaussian basis function
 * @param alphaC third gaussian power parameter
 * @param RC the center of the third gaussian
 * @param nCx the power of (x-RCx) in front of the third Gaussian basis function
 * @param nCy the power of (y-RCy) in front of the third Gaussian basis function
 * @param nCz the power of (z-RCz) in front of the third Gaussian basis function
 * @param alphaD fourth gaussian power parameter
 * @param RD the center of the fourth gaussian
 * @param nDx the power of (x-RDx) in front of the fourth Gaussian basis function
 * @param nDy the power of (y-RDy) in front of the fourth Gaussian basis function
 * @param nDz the power of (z-RDz) in front of the fourth Gaussian basis function
 * @return the electron-electron interaction matrix element of a step in the recursion
 */
double MxElemFiller::ElElRecursion(int m, double alphaA, R & RA, int nAx, int nAy, int nAz, double alphaB, R & RB, int nBx, int nBy, int nBz, double alphaC, R & RC, int nCx, int nCy, int nCz, double alphaD, R & RD, int nDx, int nDy, int nDz){

   if ((nAx<0)||(nAy<0)||(nAz<0)||(nBx<0)||(nBy<0)||(nBz<0)||(nCx<0)||(nCy<0)||(nCz<0)||(nDx<0)||(nDy<0)||(nDz<0))
      return 0.0;

   double zeta = alphaA+alphaB;
   double eta = alphaC+alphaD;
   double rho = zeta*eta/(zeta+eta);

   R P((alphaA*RA.gxco()+alphaB*RB.gxco())/zeta,(alphaA*RA.gyco()+alphaB*RB.gyco())/zeta,(alphaA*RA.gzco()+alphaB*RB.gzco())/zeta);
   R Q((alphaC*RC.gxco()+alphaD*RD.gxco())/eta, (alphaC*RC.gyco()+alphaD*RD.gyco())/eta, (alphaC*RC.gzco()+alphaD*RD.gzco())/eta);
   R W((zeta*P.gxco()+eta*Q.gxco())/(zeta+eta),(zeta*P.gyco()+eta*Q.gyco())/(zeta+eta),(zeta*P.gzco()+eta*Q.gzco())/(zeta+eta));

   if (nAx==0){
      if (nBx==0){
         if (nCx==0){
            if (nDx==0){
               if (nAy==0){
                  if (nBy==0){
                     if (nCy==0){
                        if (nDy==0){
                           if (nAz==0){
                              if (nBz==0){
                                 if (nCz==0){
                                    if (nDz==0){

                                       //all coefficients zero
                                       return 2.0/(sqrt(zeta + eta)*zeta*eta)*FunctionF(m,rho*P.DistanceSquared(Q))*pow(M_PI,2.5)*exp(-alphaA*alphaB/(alphaA+alphaB)*RA.DistanceSquared(RB))*exp(-alphaC*alphaD/(alphaC+alphaD)*RC.DistanceSquared(RD));

                                    }

                                    //lower nDz, nAz=nBz=nCz=0
                                    return (Q.gzco()-RD.gzco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz-1) + (W.gzco()-Q.gzco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz-1) + 0.5/eta*(nDz-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz-2)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz-2)); //because nAz=nBz=nCz=0

                                 }

                                 //lower nCz, nAz=nBz=0
                                 return (Q.gzco()-RC.gzco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz-1,alphaD,RD,nDx, nDy, nDz) + (W.gzco()-Q.gzco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz-1,alphaD,RD,nDx, nDy, nDz) + 0.5/eta*(nCz-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz-2,alphaD,RD,nDx, nDy, nDz)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz-2,alphaD,RD,nDx, nDy, nDz)) + 0.5/eta*nDz*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz-1,alphaD,RD,nDx, nDy, nDz-1)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz-1,alphaD,RD,nDx, nDy, nDz-1)); //because nAz=nBz=0
      
                              }

                              //lower nBz, nAz=0
                              return (P.gzco()-RB.gzco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz-1,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gzco()-P.gzco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz-1,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/zeta*(nBz-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz-2,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)-rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz-2,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/(zeta+eta)*nCz*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz-1,alphaC,RC,nCx,nCy,nCz-1,alphaD,RD,nDx,nDy,nDz) + 0.5/(zeta+eta)*nDz*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz-1,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz-1); // because nAz=0

                           }

                           //lower nAz
                           return (P.gzco()-RA.gzco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz-1,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gzco()-P.gzco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz-1,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/zeta*(nAz-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz-2,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) - rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz-2,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/zeta*nBz*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz-1,alphaB,RB,nBx,nBy,nBz-1,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) - rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz-1,alphaB,RB,nBx,nBy,nBz-1,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/(zeta+eta)*nCz*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz-1,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz-1,alphaD,RD,nDx,nDy,nDz) + 0.5/(zeta+eta)*nDz*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz-1,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz-1);

                        }

                        //lower nDy nAy=nBy=nCy=0
                        return (Q.gyco()-RD.gyco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy-1, nDz) + (W.gyco()-Q.gyco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy-1, nDz) + 0.5/eta*(nDy-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy-2, nDz)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy-2, nDz)); //because nAy=nBy=nCy=0

                     }

                     //lower nCy, nAy=nBy=0
                     return (Q.gyco()-RC.gyco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy-1,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gyco()-Q.gyco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy-1,nCz,alphaD,RD,nDx, nDy, nDz) + 0.5/eta*(nCy-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy-2,nCz,alphaD,RD,nDx, nDy, nDz)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy-2,nCz,alphaD,RD,nDx, nDy, nDz)) + 0.5/eta*nDy*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy-1,nCz,alphaD,RD,nDx, nDy-1, nDz)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy-1,nCz,alphaD,RD,nDx, nDy-1, nDz)); //because nAy=nBy=0

                  }

                  //lower nBy, nAy=0
                  return (P.gyco()-RB.gyco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy-1,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gyco()-P.gyco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy-1,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/zeta*(nBy-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy-2,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)-rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy-2,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/(zeta+eta)*nCy*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy-1,nBz,alphaC,RC,nCx,nCy-1,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/(zeta+eta)*nDy*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy-1,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy-1,nDz); // because nAy=0

               }

               //lower nAy
               return (P.gyco()-RA.gyco())*ElElRecursion(m,alphaA,RA,nAx,nAy-1,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gyco()-P.gyco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy-1,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/zeta*(nAy-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy-2,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) - rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx,nAy-2,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/zeta*nBy*(ElElRecursion(m,alphaA,RA,nAx,nAy-1,nAz,alphaB,RB,nBx,nBy-1,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) - rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx,nAy-1,nAz,alphaB,RB,nBx,nBy-1,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/(zeta+eta)*nCy*ElElRecursion(m+1,alphaA,RA,nAx,nAy-1,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy-1,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/(zeta+eta)*nDy*ElElRecursion(m+1,alphaA,RA,nAx,nAy-1,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy-1,nDz);

            }

            //nAx=nBx=nCx=0, nDx not
            return (Q.gxco()-RD.gxco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx-1,nDy,nDz) + (W.gxco()-Q.gxco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx-1,nDy,nDz) + 0.5/eta*(nDx-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx-2, nDy, nDz)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx-2, nDy, nDz)); //because nAx=nBx=nCx=0

         }

         //nAx=nBx=0, nCx not
         return (Q.gxco()-RC.gxco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx-1,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gxco()-Q.gxco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx-1,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + 0.5/eta*(nCx-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx-2,nCy,nCz,alphaD,RD,nDx, nDy, nDz)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx-2,nCy,nCz,alphaD,RD,nDx, nDy, nDz)) + 0.5/eta*nDx*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx-1,nCy,nCz,alphaD,RD,nDx-1, nDy, nDz)-rho/eta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx-1,nCy,nCz,alphaD,RD,nDx-1, nDy, nDz)); //because nAx=nBx=0

      }

      //nAx=0, nBx not
      return (P.gxco()-RB.gxco())*ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx-1,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gxco()-P.gxco())*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx-1,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/zeta*(nBx-1)*(ElElRecursion(m,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx-2,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)-rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx-2,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/(zeta+eta)*nCx*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx-1,nBy,nBz,alphaC,RC,nCx-1,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/(zeta+eta)*nDx*ElElRecursion(m+1,alphaA,RA,nAx,nAy,nAz,alphaB,RB,nBx-1,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx-1,nDy,nDz); // because nAx=0

   }

   //nAx not zero
   return (P.gxco()-RA.gxco())*ElElRecursion(m,alphaA,RA,nAx-1,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx, nDy, nDz) + (W.gxco()-P.gxco())*ElElRecursion(m+1,alphaA,RA,nAx-1,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/zeta*(nAx-1)*(ElElRecursion(m,alphaA,RA,nAx-2,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) - rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx-2,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/zeta*nBx*(ElElRecursion(m,alphaA,RA,nAx-1,nAy,nAz,alphaB,RB,nBx-1,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz) - rho/zeta*ElElRecursion(m+1,alphaA,RA,nAx-1,nAy,nAz,alphaB,RB,nBx-1,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx,nDy,nDz)) + 0.5/(zeta+eta)*nCx*ElElRecursion(m+1,alphaA,RA,nAx-1,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx-1,nCy,nCz,alphaD,RD,nDx,nDy,nDz) + 0.5/(zeta+eta)*nDx*ElElRecursion(m+1,alphaA,RA,nAx-1,nAy,nAz,alphaB,RB,nBx,nBy,nBz,alphaC,RC,nCx,nCy,nCz,alphaD,RD,nDx-1,nDy,nDz);

}


/**
 * Calculator of Fm(U) function (see article)
 * @param m the m parameter of Fm(U)
 * @param U the U parameter of Fm(U)
 * @return Fm(U)
 */
double MxElemFiller::FunctionF(int m, double U){

   if (U<=1e-12)
      return 1.0/(2.0*m+1.0) - U/(2.0*m+3.0) + U*U/(2.0*(2.0*m+5.0));

   if (m<0) {
      cout << "m<0 in FunctionF @ MxElemFiller." << endl;
      assert(m>=0);
   }

   if (m==0)
      return 0.5*sqrt(M_PI/U)*erf(sqrt(U));

   return 0.5*((2*m-1)*FunctionF(m-1,U)-exp(-U))/U;

}


/**
 * Calculate the normalisation constant of a Gaussian
 * @param alpha the gaussian power parameter
 * @param nx the power of the term (x-Rx) in front of the gaussian
 * @param ny the power of the term (y-Ry) in front of the gaussian
 * @param nz the power of the term (z-Rz) in front of the gaussian
 * @return the normalisation constant
 */
double MxElemFiller::NormCst(double alpha, int nx, int ny, int nz){

   return pow(2.0*alpha/M_PI,0.75) * pow(4.0*alpha,0.5*(nx+ny+nz)) * FactHelper(nx) * FactHelper(ny) * FactHelper(nz);

}


/**
 * Calculate 1/sqrt((2*ni-1)!!)
 * @param ni integer input
 * @return 1/sqrt((2*ni-1)!!)
 */
double MxElemFiller::FactHelper(int ni){

   if (ni==0)
      return 1.0;

   unsigned long int n = 2*ni-1;
   mpz_class fact;
   mpz_fac_ui(fact.get_mpz_t(),n);
   n = fact.get_ui();
   mpz_fac_ui(fact.get_mpz_t(),n);

   signed long int expon;

   double result = mpz_get_d_2exp(&expon,fact.get_mpz_t());
   return pow(1.0/result,0.5)*pow(0.5,0.5*double(expon));

}



