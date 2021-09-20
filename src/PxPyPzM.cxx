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

#include "PxPyPzM.h"
#include <iostream>
#include <cmath>

namespace ABCFit{
  //translate from this represenation to internal representation                                                                                              f(p
  //internal = PxPyPzM
  Coordinates PxPyPzM::Transform(Coordinates p) {
    if (p.size()!=4) std::cout<< "coordinate vector wrong dimension" << std::endl;
    if (p.size()!=4) return Coordinates();
    return p;
  }
  //translate from internal coordinate representation to this represenation
  Coordinates PxPyPzM::InvTransform(Coordinates p) {
    return p;
  }
  //Derivative of internal representation w.r.t. this representation 
  std::vector<Coordinates> PxPyPzM::Derivative(Coordinates p) {
    std::vector<Coordinates> result;
    Coordinates UnitVec (4,0.0);
    for (unsigned int i=0; i<UnitVec.size();i++) {
      UnitVec[i]=1.0;
      result.push_back(UnitVec);
      UnitVec[i]=0.0;
    }
    return result;
  }
  //Derivative of this representation w.r.t. internal representation
  std::vector<Coordinates> PxPyPzM::InvDerivative(Coordinates p) {
    std::vector<Coordinates> result;
    Coordinates UnitVec (4,0.0);
    for (unsigned int i=0; i<UnitVec.size();i++) {
      UnitVec[i]=1.0;
      result.push_back(UnitVec);
      UnitVec[i]=0.0;
    }
    return result;
  }

  Coordinates PxPyPzM::GetDefaultExpectation(Coordinates p) {
    return p;
  }

  Matrix PxPyPzM::GetDefaultExpectationCovMatrix(Coordinates p) {
    Matrix result(4,{0.0,0.0,0.0,0.0});
    double varE = 0.15*0.15*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
    result[0][0] = varE;
    result[1][1] = varE;
    result[2][2] = varE;
    result[3][3] = varE;
    return result;
  }

}
