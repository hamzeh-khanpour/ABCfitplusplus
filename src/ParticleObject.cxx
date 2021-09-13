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

#include "ParticleObject.h"
#include "MatrixAlgebra.h"
#include <iostream>
namespace ABCFit {

  ParticleObject::ParticleObject(Coordinates coordinates, CoordRepr* repr_input) : m_InputRepr(repr_input), m_Chi2Repr(internalRepresentation), m_Unmeasured(true) {
    setup(coordinates, repr_input);
  }

  ParticleObject::ParticleObject(Coordinates coordinates, CoordRepr* repr_input, CoordRepr* Chi2Repr) : m_InputRepr(repr_input), m_Chi2Repr(Chi2Repr), m_Unmeasured(false) {
    setup(coordinates, repr_input);
  }

  void ParticleObject::SetCoordinates(Coordinates coordinates) {
    Coordinates internal = m_InputRepr->Transform(coordinates);
    if (m_InputRepr == m_Chi2Repr)
      m_Coordinates = internal;
    else {
      m_Coordinates = m_Chi2Repr->InvTransform(internal); // calculate auxiliary coordinates for Chi2Repr (e.g. unit vectors)
      for (unsigned int i=0; i<m_NumberOfInternalCoordinates; i++) 
	m_Coordinates[i]=internal[i]; //First m_NumberOfInternalCoordinates coordinates are allocated for internal representation (currently 4-momentum), auxiliary coordinates should lie afterwards
    }

  }

  void ParticleObject::setup(Coordinates coordinates, CoordRepr* repr_input) {
    Coordinates internal = repr_input->Transform(coordinates);
    if (repr_input == m_Chi2Repr)
      m_Coordinates = internal;
    else {
      m_Coordinates = m_Chi2Repr->InvTransform(internal); // calculate auxiliary coordinates for Chi2Repr (e.g. unit vectors)
      for (unsigned int i=0; i<m_NumberOfInternalCoordinates; i++) 
	m_Coordinates[i]=internal[i]; //First m_NumberOfInternalCoordinates coordinates are allocated for internal representation (currently 4-momentum), auxiliary coordinates should lie afterwards
    }
    m_FixParameter = std::vector<bool> (m_Chi2Repr->GetNumberOfParameters(), false);    
    //NOTE/ATTENTION/WARNING: Since internal representation uses mass, we need to protect for zero masses:
    if (m_Chi2Repr->GetNumberOfParameters()>3)
      if (internal[3]<MassTolerance) m_FixParameter[3]=true; //assumes 4th parameter in chi2 representation is the mass parameter!!
  }

  void ParticleObject::UpdateCoordinates(Coordinates parameters) {
    Coordinates currentParameters = m_Chi2Repr->InvTransform(m_Coordinates);
    //std::cout << "currentParameters = "; MatrixAlgebra::print(currentParameters);
    unsigned int count = 0;
    for (unsigned int i=0; i< m_FixParameter.size(); i++)
      if (!m_FixParameter[i] && count<parameters.size()) currentParameters[i]=parameters[count++];
    //std::cout << "currentParameters after = "; MatrixAlgebra::print(currentParameters);
    m_Coordinates = m_Chi2Repr->Transform(currentParameters);
    //std::cout << "coodinates = "; MatrixAlgebra::print(m_Coordinates);
  }
}
