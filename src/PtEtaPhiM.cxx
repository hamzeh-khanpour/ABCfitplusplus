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

#include "PtEtaPhiM.h"
#include <iostream>
#include <cmath>
namespace ABCFit{
  //translate from this represenation to internal representation                                                                                              f(p
  //internal = PxPyPzM
  Coordinates PtEtaPhiM::Transform(Coordinates p) {
    if (p.size()!=4) std::cout<< "coordinate vector wrong dimension" << std::endl;
    if (p.size()!=4) return Coordinates();
    Coordinates result{p[0]*cos(p[2]), p[0]*sin(p[2]), p[0]*sinh(p[1]), p[3]};
    return result;
  }
  //translate from internal coordinate representation to this represenation
  Coordinates PtEtaPhiM::InvTransform(Coordinates p) {
    Coordinates result{sqrt(p[0]*p[0]+p[1]*p[1]), asinh(p[2]/sqrt(p[0]*p[0]+p[1]*p[1])), atan2(p[1],p[0]) ,p[3]};
    return result;
  }
  //Derivative of internal representation w.r.t. this representation 
  std::vector<Coordinates> PtEtaPhiM::Derivative(Coordinates p) {
    std::vector<Coordinates> result;   
    result.push_back({cos(p[2]), sin(p[2]), sinh(p[1]), 0.0});
    result.push_back({0.0, 0.0, p[0]*cosh(p[1]), 0.0});
    result.push_back({-p[0]*sin(p[2]), p[0]*cos(p[2]), 0.0, 0.0});
    result.push_back({0.0, 0.0, 0.0, 1.0});
    return result;
  }

  //Derivative of this representation w.r.t. internal representation
  std::vector<Coordinates> PtEtaPhiM::InvDerivative(Coordinates p) {
    std::vector<Coordinates> result;
    //not sure whether mass is explicitly dependent on {px,py,pz} i.e M=sqrt(E^2-p^2) hence dM/dp_i=-p_i/sqrt(E^2-p^2)=-p_i/M???
    result.push_back({p[0]/sqrt(p[0]*p[0]+p[1]*p[1]), -p[2]*p[0]/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/(p[0]*p[0]+p[1]*p[1]), -p[1]/(p[0]*p[0]+p[1]*p[1]), 0.0});
    result.push_back({p[1]/sqrt(p[0]*p[0]+p[1]*p[1]), -p[2]*p[1]/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/(p[0]*p[0]+p[1]*p[1]), p[0]/(p[0]*p[0]+p[1]*p[1]), 0.0});
    result.push_back({0.0, 1/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]), 0.0, 0.0});
    result.push_back({0.0, 0.0, 0.0, 1.0});
    return result;
  }

  //assumes Coordinate Representation in internal representation
  Coordinates PtEtaPhiM::GetDefaultExpectation(Coordinates p) {
    return p;
  }
  //assumes Coordinate Representation in internal representation
  Matrix PtEtaPhiM::GetDefaultExpectationCovMatrix(Coordinates p) {
    Matrix result(4,{0.0,0.0,0.0,0.0});
    double varpT = 0.01*(p[0]*p[0]+p[1]*p[1]); //10 % for pT
    result[0][0] = varpT/2.0;
    result[1][1] = varpT/2.0;
    result[2][2] = varpT/sin(atan2(p[1],p[0]));
    result[3][3] = 0.0;
    return result;
  }

}
