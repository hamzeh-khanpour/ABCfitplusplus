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

#include "SumPyConstraint.h"

namespace ABCFit{

  double SumPyConstraint::ConstraintFunction() { 
    std::vector<Coordinates> Particles = GetParticleCoordinates();
    double sum=0.0;
    for (auto p: Particles) sum+=p.at(1);
    return sum;
  }
  Coordinates SumPyConstraint::ConstraintFunctionDerivative() {
    std::vector<Coordinates> Particles = GetParticleCoordinates();
    Coordinates result(4*Particles.size(),0.0);//Chosen parametrisation has four coordinates per particle 
    for (unsigned int i=1; i<result.size(); i+=4) result[i]=1.0;
    return result;
  }
}
