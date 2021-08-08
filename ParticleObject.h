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

#ifndef _PARTICLEOBJECT_
#define _PARTICLEOBJECT_
#include "CoordRepr.h"
#include "PxPyPzM.h"
#include <algorithm>

namespace ABCFit {

class ParticleObject {
public:
  ParticleObject(Coordinates coordinates, CoordRepr *repr_input,
                 CoordRepr *Chi2Repr); // Measured particles constribute to chi2
  ParticleObject(Coordinates coordinates,
                 CoordRepr *repr_input); // Unmeasured particles do not
  Coordinates GetCoordinates(CoordRepr *convertTo = internalRepresentation) {
    return convertTo == internalRepresentation
               ? m_Coordinates
               : convertTo->InvTransform(m_Coordinates);
  }; // Returns coordinates by default in internal representation
  CoordRepr *GetParametrisation() { return m_Chi2Repr; };
  void UpdateCoordinates(
      Coordinates parameters); // Update momenta using chi2Repr parameter values
  bool IsUnmeasured() { return m_Unmeasured; };
  unsigned int TotalNumberOfParameters() const {
    return m_FixParameter.size();
  }; // Total number of parameters
  void FixParameter(unsigned int ipar) { m_FixParameter[ipar] = true; };
  void ReleaseParameter(unsigned int ipar) { m_FixParameter[ipar] = false; };
  bool IsParameterFixed(unsigned int ipar) { return m_FixParameter[ipar]; };
  unsigned int NumberOfFreeParameters() {
    return std::count(m_FixParameter.begin(), m_FixParameter.end(), false);
  };

private:
  Coordinates m_Coordinates;
  CoordRepr *m_Chi2Repr; // coordrepresentation with used in chi2-calculation
  bool m_Unmeasured;
  std::vector<bool>
      m_FixParameter; // Allows to remove some coordinates from chi2
  void setup(Coordinates coordinates, CoordRepr *repr_input);
};
typedef std::vector<ParticleObject *> ListOfParticleObjects;
} // namespace ABCFit
#endif
