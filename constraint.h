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

#ifndef _CONSTRAINT_
#define _CONSTRAINT_
#include "CoordRepr.h"
#include "ParticleObject.h"
#include "ProbDistFunc.h"
#include <vector>
// Use constraints to determine initial momenta of unmeasured particles
namespace ABCFit {

class Constraint {
public:
  Constraint(CoordRepr *CoordType, ListOfParticleObjects ListOfParticles,
             double ConstraintValue = 0.0);
  Constraint(CoordRepr *CoordType, ListOfParticleObjects ListOfParticles,
             ProbDistFunc *ConstraintValueDist);
  virtual double ConstraintFunction() = 0;
  virtual Coordinates ConstraintFunctionDerivative() = 0;
  Coordinates ConstraintInternalDerivative();
  bool HasProbDistFunc() {
    return m_ConstraintValueDist;
  }; // Set which constraint to use
  CoordRepr *GetCoordType() { return m_CoordType; };
  std::vector<Coordinates> GetParticleCoordinates();
  ListOfParticleObjects GetListOfParticleObjects() {
    return m_ListOfParticles;
  };
  double GetConstraintValue() { return m_ConstraintValue; };
  ProbDistFunc *GetProbDistFunc() {
    return m_ConstraintValueDist;
  }; // to check that ProbDistFunc is defined
  virtual std::string name() const = 0;

private:
  CoordRepr *m_CoordType;
  ListOfParticleObjects m_ListOfParticles;
  double m_ConstraintValue;
  ProbDistFunc *m_ConstraintValueDist;
};
typedef std::vector<std::pair<Constraint *, double>>
    ListOfConstraints; // linear combination of constraints
} // namespace ABCFit
#endif
