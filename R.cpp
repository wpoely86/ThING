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

#include "preamble.h"
#include "R.h"

/** 
 * Constructor for the R class
 * @param xco the x coordinate
 * @param yco the y coordinate
 * @param zco the z coordinate
 */
R::R(double xco, double yco, double zco){

   this->xco = xco;
   this->yco = yco;
   this->zco = zco;

}


/** 
 * Copy contructor for the R class
 * @param R_c the R to be copied
 */
R::R(R & R_c){

   this->xco = R_c.gxco();
   this->yco = R_c.gyco();
   this->zco = R_c.gzco();

}


/**
 * Standard destructor
 */
R::~R(){

}


/**
 * Ostream overloader
 * @param output the ostream to write the info to
 * @param R_p the R of which the coordinates need to be written
 */
ostream &operator<<(ostream &output, R & R_p){

   output << "[ " << R_p.gxco() << " ; " << R_p.gyco() << " ; " << R_p.gzco() << " ]";
   return output;

}


/**
 * Getter function for the x coordinate
 * @return xco the x coordinate
 */
double R::gxco(){

   return this->xco;

}


/**
 * Getter function for the y coordinate
 * @return yco the y coordinate
 */
double R::gyco(){

   return this->yco;

}


/**
 * Getter function for the z coordinate
 * @return zco the z coordinate
 */
double R::gzco(){

   return this->zco;

}


/**
 * Setter function for the coordinates
 * @param xco x coordinate
 * @param yco y coordinate
 * @param zco z coordinate
 */
void R::set(double xco, double yco, double zco){

   this->xco = xco;
   this->yco = yco;
   this->zco = zco;

}


/**
 * Calculate the distance between this R and the R Rd
 * @param Rd the other R
 * @return the distance between this R and the R Rd
 */
double R::DistanceSquared(R & Rd){

   return (xco-Rd.gxco())*(xco-Rd.gxco()) + (yco-Rd.gyco())*(yco-Rd.gyco()) + (zco-Rd.gzco())*(zco-Rd.gzco());

}


/**
 * Calculate the product between this R and the R Rd
 * @param Rd the other R
 * @return the product between this R and the R Rd
 */
double R::product(R & Rd){

   return xco * Rd.gxco() + yco * Rd.gyco() + zco * Rd.gzco();

}





