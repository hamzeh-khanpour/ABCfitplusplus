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

#ifndef _FITKERNEL_
#define _FITKERNEL_
#include "CompositeConstraint.h"
#include "CoordRepr.h"
#include "ParticleObject.h"

namespace ABCFit {
enum statuscode {
  success = 0,
  failure = -1,
  MaxNumIter = 1,
  Chi2Overflow = -2,
  NonInvertibleCovMatrix = -3,
  InconsistentExpectationValues = -4,
  NonInvertibleMatrix = -5,
  InconsistentDerivativeMatrix = -6,
  Initiate = 2,
  InitiateFailure = -7
};
struct FitResult {
  FitResult(std::vector<ParticleObject> FittedParticles, double chi2, int Ndof,
            statuscode status, unsigned int Niter);
  void print(bool PrintParticles = false);
  std::vector<ParticleObject> FittedParticles;
  double chi2;
  int Ndof;
  statuscode status;
  unsigned int Niter;
};

class FitKernel {
public:
  FitKernel(std::vector<CompositeConstraint> ListOfConstraints,
            unsigned int maxiter = 100);
  FitResult Fit(unsigned int maxiter = 0);
  void SetMaxIterations(unsigned int maxiter = 100) {
    m_MaxIterations = maxiter;
  };
  Matrix GetParameterCovMatrix(); // Full covariance matrix of all internal chi2
                                  // parameters = measured + unmeasured
                                  // parameters in that order
  Matrix GetParticleCovMatrix(); // covariance matrix for fitted momenta in
                                 // internal representation (PxPyPxM)
  Coordinates GetMeasuredParameters() { return m_ParameterValues; };
  Coordinates GetUnmeasuredMeasuredParameters() {
    return m_UnmeasuredParameterValues;
  };
  bool Initialize();

private:
  std::vector<CompositeConstraint> m_ListOfConstraints;
  unsigned int m_MaxIterations;
  statuscode m_State;
  double m_DeltaChiConvergence = 0.0005;
  double m_ConstraintConvergence = 0.0001;
  ListOfParticleObjects m_AllParticles;
  std::vector<unsigned int> m_Particle2Index;
  unsigned int m_NumberOfConstraints;
  unsigned int m_NumberOfFreeParameters;
  unsigned int m_NumberOfParameters;
  unsigned int m_NumberOfConstraintParameters;
  unsigned int m_TotalNumberOfParticleParameters;
  std::vector<std::vector<unsigned int>> m_Constraint2Particle;
  std::vector<unsigned int> m_HasProbDistFunc;
  Coordinates m_ConstraintValues;
  Coordinates m_ParameterValues;
  Coordinates m_UnmeasuredParameterValues;
  Matrix m_dpdParameters;
  Matrix m_Amatrix;
  Matrix m_Bmatrix;
  Matrix m_Covariance;
  Matrix m_WAWB;
  Matrix m_vbtWb;
  Matrix m_WaInv;
  statuscode CalculateDerivativeMatrix();
};

} // namespace ABCFit
#endif
