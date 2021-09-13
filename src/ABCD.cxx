//  ABCfit++ -- General constrained kinematic fit using ABC-parametrisation
//  Copyright (C) 2021
//    - Jørgen Beck Hansen, <beck@nbi.ku.dk>
//    - Julie Munch Torndal, <julie@torndal.com>
//
//  This file is part of ABCfit++.
//
//  ABCfit++ is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  ABCfit++ is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "ABCD.h"
#include "MatrixAlgebra.h"
#include <iostream>
#include <cmath>

namespace ABCFit{
  //translate from this represenation to internal representation
  //internal = PxPyPzM
  Coordinates ABCD::Transform(Coordinates p) {
    if (p.size()!=GetNumberOfCoordinates()) {
      std::cout<< "coordinate vector wrong dimension" << std::endl;
      return Coordinates();
    }
    Coordinates result(GetNumberOfCoordinates(),0.0);
    std::copy(p.begin()+4,p.end(),result.begin()+4); // preserve unit vector definition!
    //result[0] = p[0]*result[4]+p[1]*result[7]+p[2]*result[9];
    //result[1] = p[0]*result[5]+p[1]*result[8]+p[2]*result[10];
    //result[2] = p[0]*result[6]+p[2]*result[11];
    result[0] = p[0]*result[4]+p[1]*result[9]+p[2]*result[7];
    result[1] = p[0]*result[5]+p[1]*result[10]+p[2]*result[8];
    result[2] = p[0]*result[6]+p[1]*result[11];
    //result[3] = p[3];
    if (m_logscale) 
      result[3] = exp(p[3])*p[12];
    else            
      result[3] = p[3]*p[12];
    if (this->GetNumberOfParameters()==3) result[3]=sqrt(p[0]*p[0]+(p[1]*p[1]+p[2]*p[2])/(result[4]*result[4]+result[5]*result[5]+result[6]*result[6]))*p[12];
    return result;
  }
  //translate from internal coordinate representation to this represenation
  Coordinates ABCD::InvTransform(Coordinates p, bool RecalculateUnitVectors) {
    Coordinates result(GetNumberOfCoordinates(),0.0);
    double momentum = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    if (p.size()==GetNumberOfCoordinates() && !RecalculateUnitVectors) {
      std::copy(p.begin()+4,p.end(),result.begin()+4); // preserve unit vector definition!
      result[0]=(p[0]*result[4]+p[1]*result[5]+p[2]*result[6])/(result[4]*result[4]+result[5]*result[5]+result[6]*result[6]);
      result[2]=p[0]*result[7]+p[1]*result[8];
      result[1]=p[0]*result[9]+p[1]*result[10]+p[2]*result[11];
      if (p[12]<MassTolerance ) {
	if (!m_logscale) result[3]=1.0;
	return result;
      }
      if (m_logscale)
	result[3] = log(p[3]/p[12]);
      else
	result[3] = p[3]/p[12];
      return result;
    }
    //Pa vector (obs not normalised)
    result[4]=p[0]; 
    result[5]=p[1];
    result[6]=p[2];
    //Pb unit vector
    result[7]=-p[1];
    result[8]=p[0];
    double norm = sqrt(result[7]*result[7]+result[8]*result[8]);
    if (norm < MassTolerance ) {
      result[7]=1.0;
      result[8]=0.0;
    }
    else {
      result[7]=result[7]/norm;
      result[8]=result[8]/norm;
    }
    //Pc unit vector = Pa x Pb (normalised by construction)
    result[9] = result[6]*result[8]/momentum;
    result[10] = -result[6]*result[7]/momentum;
    result[11] = -(result[4]*result[8]-result[5]*result[7])/momentum;
    //ABCD
    result[0]=1.0;
    result[1]=0.0;
    result[2]=0.0;
    if (!m_logscale) result[3]=1.0;
    result[12]=p[3]; //Save input mass for mass scaling option
    return result;
  }
  //Derivative of internal representation w.r.t. this representation 
  std::vector<Coordinates> ABCD::Derivative(Coordinates p) {
    if (p.size()!=GetNumberOfCoordinates()) std::cout<< "coordinate vector wrong dimension" << std::endl;
    std::vector<Coordinates> result;
    if (this->GetNumberOfParameters()==3) {
      result.push_back({p[4],p[9],p[7]});
      result.push_back({p[5],p[10],p[8]});
      result.push_back({p[6],p[11],0.0});
      double temp =p[4]*p[4]+p[5]*p[5]+p[6]*p[6];
      double norm = sqrt(p[0]*p[0]+(p[1]*p[1]+p[2]*p[2])/temp)*p[12];
      if(abs(p[12])>MassTolerance ) result.push_back({p[0]/norm,p[1]/norm/temp,p[2]/norm/temp});
      else result.push_back({0.0,0.0,0.0});
      return result;
    }
    result.push_back({p[4],p[9],p[7],0.0});
    result.push_back({p[5],p[10],p[8],0.0});
    result.push_back({p[6],p[11],0.0,0.0});
    if (m_logscale) 
      result.push_back({0.0,0.0,0.0,exp(p[3])*p[12]});
    else
      result.push_back({0.0,0.0,0.0,p[12]});
    return result;
  }
  //Derivative of this representation w.r.t. internal representation
  std::vector<Coordinates> ABCD::InvDerivative(Coordinates p) {
    if (p.size()!=GetNumberOfCoordinates()) std::cout<< "coordinate vector wrong dimension" << std::endl;
    std::vector<Coordinates> result;
    result.push_back({p[4],p[5],p[6],0.0});
    result.push_back({p[9],p[10],p[11],0.0});
    result.push_back({p[7],p[8],0.0,0.0});
    if (this->GetNumberOfParameters()==3) return result;
    if(abs(p[12])<MassTolerance) {
      result.push_back({0.0,0.0,0.0,0.0});
      return result;
    }
    if (m_logscale) 
      result.push_back({0.0,0.0,0.0,1.0/p[3]});
    else
      result.push_back({0.0,0.0,0.0,1.0/p[12]});
    return result;
  }
  //assumes Coordinate Representation in internal representation
  Coordinates ABCD::GetDefaultExpectation(Coordinates p) {
    if (this->GetNumberOfParameters()==3) return Coordinates ({1.0,0.0,0.0});
    return Coordinates ({1.0,0.0,0.0,p[3]});
  }
  //assumes Coordinate Representation in internal representation
  Matrix ABCD::GetDefaultExpectationCovMatrix(Coordinates p) {
    Matrix result;
    double varE = 0.15*0.15/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]); //sampling term 15%
    if (this->GetNumberOfParameters()==3) {
      result = Matrix(3,{0.0,0.0,0.0});
      result[0][0] = varE;
      result[1][1] = 0.09;
      result[2][2] = 0.09; //Lambda_QCD^2 quark pT
      return result;
    }
    result = Matrix(4,{0.0,0.0,0.0,0.0});
    result[0][0] = varE;
    result[1][1] = 0.09;
    result[2][2] = 0.09; //Lambda_QCD^2 quark pT
    if(abs(p[3])>MassTolerance) result[3][3] = 0.15*0.15/p[3]; //assuming basic scaling
    return result;
  }

}
