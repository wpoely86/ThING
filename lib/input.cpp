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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <assert.h>
#include "input.h"
#include "Gauss.h"

using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::ifstream;
using std::string;
using std::vector;

using namespace ThING;

int MENDELJEV = 118;


/**
 * Constructor for the input class
 * @param setupfile filename of the setupfile of the job
 */
input::input(string setupfile){

   initelements();
   readinsetupfile(setupfile);
   fillgaussinfo();

}


/**
 * Copy constructor for the input class
 * @param lisa the input to be copied
 */
input::input(input & lisa){

   initelements();
   Charge = lisa.gCharge();
   RotationSymm = lisa.gRotationSymm();
   basisset = lisa.gbasisset();
   Ncores = lisa.gNcores();
   cores = new int[Ncores];
   vectors = new R*[Ncores];
   GaussInfo = new Gauss*[Ncores];
   for (int cnt=0; cnt<Ncores; cnt++){
      cores[cnt] = lisa.gcore(cnt);
      vectors[cnt] = new R(*(lisa.gvector(cnt)));
      GaussInfo[cnt] = new Gauss(*(lisa.gGaussInfo(cnt)));
   }

}


/**
 * Standard destructor to deallocate the memory
 */
input::~input(){

   delete [] elements;
   delete [] cores;
   for (int cnt=0; cnt<Ncores; cnt++){
      delete vectors[cnt];
      delete GaussInfo[cnt];
   }
   delete [] vectors;
   delete [] GaussInfo;

}


/**
 * Getter function for the global charge
 * @return the charge
 */
int input::gCharge(){

   return Charge;

}


/**
 * Getter function for the rotation symmetry
 * @return the rotation symmetry
 */
char input::gRotationSymm(){

   return RotationSymm;

}


/**
 * Getter function for the number of cores in the problem
 * @return the number of cores in the problem
 */
int input::gNcores(){

   return Ncores;

}


/**
 * Getter function for the basis set name
 * @return the basis set name
 */
string input::gbasisset(){

   return basisset;

}


/**
 * Getter function for charge of the i-th core
 * @param i the number of the core
 * @return the charge of the i-th core
 */
int input::gcore(int i){

   return cores[i];

}


/**
 * Getter function for the position vector of the i-th core
 * @param i the number of the core
 * @return the position vector of the i-th core
 */
R* input::gvector(int i){

   return vectors[i];

}


/**
 * Getter function for the Gauss information object of the i-th core
 * @param i the number of the core
 * @return the Gauss information object of the i-th core
 */
Gauss* input::gGaussInfo(int cnt){

   return GaussInfo[cnt];

}


/**
 * Function to read in the elements mapper : elements[Z-1] = "name"
 */
void input::initelements(){

   elements = new string[MENDELJEV];
   string temp;
   int counter = 0;

   ifstream elem("basissets/elements.bs", ios::in);
   while(getline(elem,temp)) {
      int place = temp.find("\t") + 1;
      elements[counter] = temp.substr(place,temp.size()-place);
      counter++;
   }
   elem.close();

}


/**
 * Function that establishes the reverse elements mapper : getZ("name")=Z
 */
int input::getZ(string elem){

   for (int cnt=0; cnt<MENDELJEV; cnt++)
      if (!elem.compare(elements[cnt]))
         return cnt+1;

   return -1;

}

int input::NumberOfElectrons(){

   int count = 0;
   for (int i=0; i<Ncores; i++)
      count += cores[i];

   count -= Charge;
   return count;

}

/**
 * Function to read in the setup file
 * @param setupfile filename of the setup file
 */
void input::readinsetupfile(string setupfile){

   string temp, temp2, temp3;
   int pos1, pos2, pos3, length;

   ifstream inputfile(setupfile.c_str());

   //1: the basisset
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   basisset = temp.substr(pos1,length);

   //2: the charge
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   Charge = atoi((temp.substr(pos1,length)).c_str());
   
   //3: the rotation symmetry
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   RotationSymm = temp.at(pos1);

   //4: white line
   getline(inputfile,temp);

   //5: predefined distances
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   int numberofpred = atoi((temp.substr(pos1,length)).c_str());
   string * variables = new string[numberofpred];
   double * values = new double[numberofpred];

   //6: fill the arrays
   for (int cnt=0; cnt<numberofpred; cnt++){
      getline(inputfile,temp);
      pos1 = temp.find("=");
      variables[cnt] = temp.substr(0,pos1);
      pos1++;
      length = temp.size()-pos1;
      values[cnt] = atof((temp.substr(pos1,length)).c_str());
   }

   //7: white line
   getline(inputfile,temp);

   //8: number of cores;
   getline(inputfile,temp);
   pos1 = temp.find("=")+1;
   length = temp.size()-pos1;
   Ncores = atoi((temp.substr(pos1,length)).c_str());
   cores = new int[Ncores];
   vectors = new R*[Ncores];
      
   double * cotemp = new double[3];

   //9: read in cores
   for (int cnt=0; cnt<Ncores; cnt++){
      getline(inputfile,temp);

      //9.1 Core name -> Z
      pos1 = temp.find(" ")+1;
      pos2 = temp.find(" ",pos1);
      temp2 = temp.substr(pos1,pos2-pos1);
      cores[cnt] = getZ(temp2);

      //9.2 Core coordinates -> R
      for (int cnt2=0; cnt2<3; cnt2++){
         pos1 = pos2+1;
         if (cnt2<2) 
            pos2 = temp.find(" ",pos1);
         else
            pos2 = temp.size();
         temp2 = temp.substr(pos1,pos2-pos1);

         //Predefined distance used?
         pos3 = temp2.find("*");
         if (pos3==(int)string::npos){
            cotemp[cnt2] = atof(temp2.c_str());
         } else {
            temp3 = temp2.substr(pos3+1,temp2.size()-pos3-1);
            bool found = false;
            int counter = 0;
            while ((!found)&&(counter<numberofpred)){
               if (!temp3.compare(variables[counter]))
                  found = true;
               counter++;
            }
            cotemp[cnt2] = values[counter-1] * atof((temp2.substr(0,pos3)).c_str());
         }
      }

      vectors[cnt] = new R(cotemp[0],cotemp[1],cotemp[2]);
   }

   delete [] cotemp;
   delete [] variables;
   delete [] values;

   inputfile.close();

}


/**
 * Function to create the GaussInfo array and fill it.
 */
void input::fillgaussinfo(){

   GaussInfo = new Gauss*[Ncores];

   string marge = basisset;
   while (marge.find("*")!=string::npos){
      marge.replace(marge.find("*"),1,"star");
   }

   marge = "basissets/" + marge + ".bs";
   ifstream basiss(marge.c_str(), ios::in);
   int maggy = 0;
   int numberofgauss, Z, cnt, count;

   bool itchy = true;

   //jump to start of the input of the first atom
   do {
      getline(basiss,marge);
      if (marge.find("***")!=string::npos)
         itchy = false;
   } while (itchy);

   bool scratchy = false;

   //loop over the atoms in the basisset and stop if (i) maggy equals Ncores or (ii) the end of the file has been reached
   do {

      //End of file? Scratchy kills the loop.
      if (!getline(basiss,marge)) scratchy = true;
      else {
         //which atom?
         Z = getZ(marge.substr(0,marge.find(" ")));

         //atom required? if yes: cnt contains the index of the core in our problem
         cnt = 0;
         itchy = true;
         while ((itchy)&&(cnt<Ncores)) {
            if (Z==cores[cnt]) itchy = false;
            else cnt++;
         }

         //if required: readin
         if (!itchy){

            //dump everything in an ugly vector and count i.t.m. the number of different contractions
            vector<string> milhouse;
            itchy = true;
            count = 0;
            do {
               getline(basiss,marge);
               if (marge.find("***")!=string::npos) itchy = false;
               else {
                  milhouse.push_back(marge);
                  if ((marge.substr(0,2)).compare("SP")==0) count += 2;
                  else
                     if ((marge.substr(0,1)).compare(" ")) count++;
               }
            } while (itchy);

            //make a Gauss object to store the contractions
            GaussInfo[cnt] = new Gauss(count);
            maggy++;

            //loop over the vector and store everything in the Gauss object
            count = 0;
            for (int cnt2=0; cnt2<(int)milhouse.size(); cnt2++){
               marge = milhouse[cnt2];
               if ((marge.substr(0,2)).compare("SP")==0) {
                  numberofgauss = atoi((marge.substr(4,3)).c_str());
                  double * alphas = new double[numberofgauss];
                  double * fors = new double[numberofgauss];
                  double * forp = new double[numberofgauss];
                  
                  for (int cnt3=0; cnt3<numberofgauss; cnt3++){
                     marge = milhouse[cnt2+1+cnt3];

                     string n1 = marge.substr(0,27);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     alphas[cnt3] = atof(n1.c_str());

                     n1 = marge.substr(27,23);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     fors[cnt3] = atof(n1.c_str());

                     n1 = marge.substr(50,marge.size()-50);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     forp[cnt3] = atof(n1.c_str());
                  }

                  GaussInfo[cnt]->set(count,numberofgauss,alphas,fors,'S');
                  GaussInfo[cnt]->set(count+1,numberofgauss,alphas,forp,'P');
                  count+=2;

                  cnt2+=numberofgauss;

                  delete [] alphas;
                  delete [] fors;
                  delete [] forp;

               } else {
                  numberofgauss = atoi((marge.substr(4,3)).c_str());
                  double * alphas = new double[numberofgauss];
                  double * coeff = new double[numberofgauss];
                  
                  for (int cnt3=0; cnt3<numberofgauss; cnt3++){
                     marge = milhouse[cnt2+1+cnt3];

                     string n1 = marge.substr(0,27);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     alphas[cnt3] = atof(n1.c_str());

                     n1 = marge.substr(27,marge.size()-50);
                     if (n1.find("D")!=string::npos) n1.replace(n1.find("D"),1,"E");
                     coeff[cnt3] = atof(n1.c_str());
                  }

                  GaussInfo[cnt]->set(count,numberofgauss,alphas,coeff,milhouse[cnt2].at(0));
                  count++;

                  cnt2+=numberofgauss;

                  delete [] alphas;
                  delete [] coeff;
               }

            }

            //check if other cores in the problem have the same number of protons. if yes: copy the gauss object.
            count = cnt;
            cnt++;
            while (cnt<Ncores) {
               if (Z==cores[cnt]){
                  GaussInfo[cnt] = new Gauss(*(GaussInfo[count]));
                  maggy++;
               }
               cnt++;
            }
            
         //if the atom is not required move on to the next atom.
         } else {
            itchy = true;
            do {
               if (marge.find("***")!=string::npos) itchy = false;
            } while ((itchy)&&(getline(basiss,marge)));
         }
      }

   } while ((maggy<Ncores)&&(!scratchy));

   //terminate the program if the basisset doesn't contain all the required atoms
   if (maggy<Ncores){
      cout << "Basisset doesn't contain all the required atoms!" << endl;
      assert(maggy>=Ncores);
   }

   basiss.close();

}


/**
 * Printer of this input object
 */
void input::debugprint(){

   cout << "Basisset : " << basisset << endl;
   cout << "Charge : " << Charge << endl;
   cout << "Rotation symmetry : " << RotationSymm << endl;
   cout << "Cores, vectors, orbitals and their contractions : " << endl;
   for (int cnt=0; cnt<Ncores; cnt++){
      cout << "Z = " << cores[cnt] << " and R = " << *(vectors[cnt]) << endl;
      cout << *(GaussInfo[cnt]) << endl;
      cout << " " << endl;
   }

}




