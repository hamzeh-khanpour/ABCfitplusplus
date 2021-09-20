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

#include "InvMassConstraint.h"
#include "MatrixAlgebra.h"
#include <cmath>

namespace ABCFit{

  double InvMassConstraint::ConstraintFunction() { 
    std::vector<Coordinates> Particles = GetParticleCoordinates();
    Coordinates sumvector(4,0.0);
    for (auto p: Particles) sumvector=MatrixAlgebra::VectorSum(sumvector,p);
    return sqrt(std::max(sumvector[3]*sumvector[3]-(sumvector[0]*sumvector[0]+sumvector[1]*sumvector[1]+sumvector[2]*sumvector[2]),0.0));
  }

  Coordinates InvMassConstraint::ConstraintFunctionDerivative() {
    std::vector<Coordinates> Particles = GetParticleCoordinates();
    Coordinates result(4*Particles.size(),0.0);//Chosen parametrisation has four coordinates per particle 
    Coordinates sumvector(4,0.0);
    for (auto p: Particles) sumvector=MatrixAlgebra::VectorSum(sumvector,p);
    double m = sqrt(std::max(sumvector[3]*sumvector[3]-(sumvector[0]*sumvector[0]+sumvector[1]*sumvector[1]+sumvector[2]*sumvector[2]),0.0));
    for (unsigned int i=0; i<result.size(); i+=4) {
      result[i+0]=-sumvector[0]/m;
      result[i+1]=-sumvector[1]/m;
      result[i+2]=-sumvector[2]/m;
      result[i+3]=sumvector[3]/m;
    }
    return result;
  }

}
