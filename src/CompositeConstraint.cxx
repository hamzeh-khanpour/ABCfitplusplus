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

#include "CompositeConstraint.h"
#include "MatrixAlgebra.h"
#include <algorithm>
#include <vector>

namespace ABCFit{

  CompositeConstraint::CompositeConstraint(ListOfConstraints Constraints, double ConstraintValue) : m_ListOfConstraints(Constraints), m_ConstraintValue(ConstraintValue), m_ConstraintValueDist(NULL) {if (!CheckInput()) m_ListOfConstraints=ListOfConstraints();};

  CompositeConstraint::CompositeConstraint(ListOfConstraints Constraints, ProbDistFunc* ConstraintValueDist) : m_ListOfConstraints(Constraints), m_ConstraintValue(0.0), m_ConstraintValueDist(ConstraintValueDist) {if (!CheckInput()) m_ListOfConstraints=ListOfConstraints();};

  double CompositeConstraint::ConstraintFunction() {
    double result = 0.0;
    for (auto c : m_ListOfConstraints) result=result+c.second*c.first->ConstraintFunction();
    return result;
  }

  Coordinates CompositeConstraint::ConstraintFunctionDerivative() {
    Coordinates result;
    for (auto c : m_ListOfConstraints) {
      Coordinates temp = MatrixAlgebra::ScaleVector(c.first->ConstraintFunctionDerivative(),c.second);
      result.insert(result.end(), temp.begin(), temp.end());
    }
    unsigned int removed = 0;
    for (unsigned int i=0; i<m_Index.size(); i++) {
      if (m_Index[i]!=i){
	for (unsigned int j=0; j<4; j++) {
	  result[(m_Index[i]-removed)*4+j]+=result[(i-removed)*4];
	  result.erase(result.begin()+(i-removed)*4);
	}
      }
    } 
    return result;
  }
 
  Coordinates CompositeConstraint::ConstraintInternalDerivative() {
    Coordinates result;
    for (auto c : m_ListOfConstraints) {
      Coordinates temp = MatrixAlgebra::ScaleVector(c.first->ConstraintInternalDerivative(),c.second);
      result.insert(result.end(), temp.begin(), temp.end());
    }
    unsigned int removed = 0;
    for (unsigned int i=0; i<m_Index.size(); i++) {
      if (m_Index[i]!=i){
	for (unsigned int j=0; j<4; j++) {
          result[(m_Index[i]-removed)*4+j]+=result[(i-removed)*4];
          result.erase(result.begin()+(i-removed)*4);
	}
      }
    }
    return result;
  }

  bool CompositeConstraint::CheckInput() {
    std::string name = m_ListOfConstraints.size()>0 ? m_ListOfConstraints.begin()->first->name() : "undefined";
    ListOfParticleObjects AllParticles; //List of unique particles
    unsigned int count = 0;
    m_Index.clear();
    for (auto c : m_ListOfConstraints) {
      if (c.first->name()!=name) return false;
      for (auto p : c.first->GetListOfParticleObjects()){
	if (std::find(AllParticles.begin(),AllParticles.end(), p)!=AllParticles.end()) {
	  m_Index.push_back(std::find(AllParticles.begin(),AllParticles.end(), p)-AllParticles.begin());
	}
	else {
	  AllParticles.push_back(p);
	  m_Index.push_back(count);
	}
	count++;
      }
    }
    return true;
  }

}
