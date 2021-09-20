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

#ifndef _COMPOSITECONSTRAINT_
#define _COMPOSITECONSTRAINT_
#include "constraint.h"

namespace ABCFit {

class CompositeConstraint {
public:
  CompositeConstraint(ListOfConstraints Constraints,
                      double ConstraintValue = 0.0);
  CompositeConstraint(ListOfConstraints Constraints,
                      ProbDistFunc *ConstraintValueDist);
  CompositeConstraint(Constraint *constraint) {
    m_ListOfConstraints.push_back({constraint, 1.0});
    m_ConstraintValue = constraint->GetConstraintValue();
    m_ConstraintValueDist = constraint->GetProbDistFunc();
  };
  double ConstraintFunction();
  Coordinates ConstraintFunctionDerivative();
  Coordinates ConstraintInternalDerivative();
  bool HasProbDistFunc() {
    return m_ConstraintValueDist;
  }; // Set which constraint to use
  ListOfConstraints GetListOfConstraints() {
    return m_ListOfConstraints;
  }; // returns pointers to constraints
  double GetConstraintValue() { return m_ConstraintValue; };
  ProbDistFunc *GetProbDistFunc() {
    return m_ConstraintValueDist;
  }; // to check that ProbDistFunc is defined
private:
  ListOfConstraints m_ListOfConstraints;
  double m_ConstraintValue;
  ProbDistFunc *m_ConstraintValueDist;
  bool CheckInput();
  std::vector<unsigned int> m_Index;
  // particle allocator
};
} // namespace ABCFit
#endif
