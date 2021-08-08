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

#include "constraint.h"
#include "MatrixAlgebra.h"
#include <iostream>

namespace ABCFit{

  Constraint::Constraint(CoordRepr* CoordType, std::vector<ParticleObject*> ListOfParticles, double ConstraintValue) : m_CoordType(CoordType), m_ListOfParticles(ListOfParticles), m_ConstraintValue(ConstraintValue), m_ConstraintValueDist(NULL) {};

  Constraint::Constraint(CoordRepr* CoordType, std::vector<ParticleObject*> ListOfParticles,  ProbDistFunc* ConstraintValueDist) : m_CoordType(CoordType), m_ListOfParticles(ListOfParticles), m_ConstraintValue(0.0), m_ConstraintValueDist(ConstraintValueDist) {};

  std::vector<Coordinates> Constraint::GetParticleCoordinates() {
    std::vector<Coordinates> result;
    for (auto* p : m_ListOfParticles) {
      Coordinates temp = p->GetCoordinates(m_CoordType);
      result.push_back(Coordinates(temp.begin(),temp.begin()+m_CoordType->GetNumberOfCoordinates()));
     }
    return result;
  }

  Coordinates Constraint::ConstraintInternalDerivative() {
    Coordinates result;
    Coordinates vector=ConstraintFunctionDerivative(); 
    //std::cout << "ConstraintInternalDerivative = "; MatrixAlgebra::print(vector); 
    for (int i =0; i<m_ListOfParticles.size(); i++){
	Coordinates singleparticlederiv = MatrixAlgebra::VectorMatrixMultiplication(Coordinates(vector.begin()+4*i,vector.begin()+4+4*i),m_CoordType->InvDerivative(m_ListOfParticles[i]->GetCoordinates()));
	//std::cout << "Inv deriv = "; MatrixAlgebra::print(m_CoordType->InvDerivative(m_ListOfParticles[i]->GetCoordinates()));
      result.insert(result.end(), singleparticlederiv.begin(), singleparticlederiv.end());
    }
    return result;
  }

}
