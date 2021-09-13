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

#include "CoordRepr.h"
#include "PxPyPzM.h"
#include "PxPyPzE.h"
#include "PtEtaPhiM.h"
#include "ABCD.h"
#include "MatrixAlgebra.h"
#include "ParticleObject.h"
//#include "ProbDistFunc.h"
#include "GaussianPDF.h"
#include "constraint.h"
#include "SumPxConstraint.h"
#include "SumPyConstraint.h"
#include "SumPzConstraint.h"
#include "SumEConstraint.h"
#include "InvMassConstraint.h"
#include "CompositeConstraint.h"
#include "ABCFit.h"

#include <iostream>
#include <typeinfo>
#include <cmath>

using namespace ABCFit;
//***************************************************************************<*
//User specifications of particle parametrisations can be controlled by function pointers
//In the following there is a list of examples to illustrate syntax and usage
//***************************************************************************<*
//Generic example parametrisation and corresponding covariance matrix
//ABCD parametrisation 
//Argument: Particle coordinates in internal representation
//Returns: ABCD parameters
Coordinates TestParameterisation(Coordinates p) {
  Coordinates result({0.0,0.0,0.0,0.0});
  result[0] = 1.0;
  result[1] = 0.0;
  result[2] = 0.0; 
  result[3] = p[3]; 
  return result;
}
//ABCD parametrisation 
//Argument: Particle coordinates in internal representation
//Returns: Diagonal covariance matrix for ABCD parameters 
Matrix TestParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  double sigmaE = 0.15/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]); //sampling term 15%                                                                         
  result[0][0] = sigmaE*sigmaE;
  result[1][1] = 0.09;
  result[2][2] = 0.09; //Lambda_QCD^2 quark pT                                                                                                                       
  result[3][3] = sigmaE*sigmaE;
  return result;
}
//***************************************************************************<*
//In the following we show a list of type dependent parametrisations used to illustrate different particle types inside the constrained fit
//***************************************************************************<*
//Example parametrisation for leptons from ttbar @ FCC-ee
//ABC parametrisation 
//Argument: Particle coordinates in internal representation
//Returns: Diagonal covariance matrix for ABC parameters (NOTE: mass term is assumed not to be used)
Matrix LeptonParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  result[0][0] = 0.003;
  result[1][1] = 0.003;
  result[2][2] = 0.003; 
  result[3][3] = 10.0; //mass value not used  
  return result;
}
//Example parametrisation for b-jets from ttbar @ FCC-ee
//ABC parametrisation 
//Argument: Particle coordinates in internal representation
//Returns: Diagonal covariance matrix for ABC parameters (NOTE: mass term is assumed not to be used)
Matrix BjetParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  result[0][0] = 0.07;
  result[1][1] = 1.2;
  result[2][2] = 1.2; 
  result[3][3] = 10.0; //mass value not used
  return result;
}
//Example parametrisation for light flavour jets from ttbar @ FCC-ee
//ABC parametrisation 
//Argument: Particle coordinates in internal representation
//Returns: Diagonal covariance matrix for ABC parameters (NOTE: mass term is assumed not to be used)
Matrix LFjetParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  result[0][0] = 0.07;
  result[1][1] = 1.4;
  result[2][2] = 1.4; 
  result[3][3] = 10.0; //mass value not used
  return result;
}
//Example parametrisation for jets (b-jet and light-flavour jets) from ttbar @ FCC-ee
//ABC parametrisation 
//Argument: Particle coordinates in internal representation
//Returns: ABC parameters (NOTE: mass term is assumed not to be used)
Coordinates JetParameterisation(Coordinates p) {
  Coordinates result({0.0,0.0,0.0,0.0});
  result[0] = 1.03;
  result[1] = 0.0;
  result[2] = 0.0; 
  result[3] = 10.0; //mass value not used
  return result;
}

//***************************************************************************<*

//***************************************************************************<*
int main() {
  //A few basic examples on usage of Coordinate representation (CoordRepr)
  //Coordinates for test particle
  Coordinates test{1.0,2.0,3.0,4.0};

  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Example of CoordRepr with PxPyPzE and PtEtaPhiM" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;

  //Input coordinates
  std::cout << "Input coordinates = "; MatrixAlgebra::print(test);
  //Testing for class PxPyPzE
  Coordinates outPxPyPzE1 = PxPyPzE().Transform(test);
  Coordinates outPxPyPzE2 = PxPyPzE().InvTransform(test);
  std::cout << "--------- PxPyPzE class ---------:" << std::endl;
  std::cout << "Coordinates transformation from internal representation (input coordinates) to PxPyPzE: "; MatrixAlgebra::print(outPxPyPzE1);
  std::cout << "Coordinates transformation from PxPyPzE (input coordinates) to internal representation: "; MatrixAlgebra::print(outPxPyPzE2);
  std::cout << "Derivative of internal representation w.r.t. PxPyPzE representation : "; MatrixAlgebra::print(PxPyPzE().Derivative(test));
  std::cout << "Derivative of PxPyPzE representation w.r.t. internal representation : "; MatrixAlgebra::print(PxPyPzE().Derivative(test));

  //Testing for class PtEtaPhiM
  Coordinates outPtEtaPhiM1 = PtEtaPhiM().Transform(test);
  Coordinates outPtEtaPhiM2 = PtEtaPhiM().InvTransform(test);
  std::cout << "--------- PtEtaPhiM class ---------:" << std::endl;
  std::cout << "Coordinates transformation from internal representation (input coordinates) to PtEtaPhiM: "; MatrixAlgebra::print(outPtEtaPhiM1);
  std::cout << "Coordinates transformation from PtEtaPhiM (input coordinates) to internal representation: "; MatrixAlgebra::print(outPtEtaPhiM2);
  std::cout << "Derivative of internal representation w.r.t. PtEtaPhiM representation : "; MatrixAlgebra::print(PtEtaPhiM().Derivative(test));
  std::cout << "Derivative of PtEtaPhiM representation w.r.t. internal representation : "; MatrixAlgebra::print(PtEtaPhiM().Derivative(test));

  //One can transform between representation via internal representation
  PxPyPzE repPxPyPzE;
  PtEtaPhiM repPtEtaPhiM;
  Coordinates out2transforms1 = ((CoordRepr*)&repPxPyPzE)->Transform(test,(CoordRepr*)&repPtEtaPhiM);
  Coordinates out2transforms2 = ((CoordRepr*)&repPtEtaPhiM)->Transform(test,(CoordRepr*)&repPxPyPzE);
  std::cout << "Coordinates transformation from PxPyPzE (input coordinates) to PtEtaPhiM via internal representation: ";  MatrixAlgebra::print(out2transforms1);
  std::cout << "Coordinates transformation from PtEtaPhiM (input coordinates) to PxPyPzE via internal representation: ";  MatrixAlgebra::print(out2transforms2);
  std::cout << "Derivative of PtEtaPhiM representation w.r.t. PxPyPzE representation : "; MatrixAlgebra::print(((CoordRepr*)&repPxPyPzE)->Derivative(test,(CoordRepr*)&repPtEtaPhiM));
  std::cout << "Derivative of PxPyPzE representation w.r.t. PtEtaPhiM representation : "; MatrixAlgebra::print(((CoordRepr*)&repPtEtaPhiM)->Derivative(test,(CoordRepr*)&repPxPyPzE));
  
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "ProbDistFunc" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "ParticleObject" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  //Creating Particle Object with coordinates test in PxPyPzE representation and with chi2 contribution in PtEtaPhiM representation
  ParticleObject POtest(test, pxpypze, ptetaphim);
  Coordinates POcoord=POtest.GetCoordinates(); //Coordinates in internal representation
  std::cout << "Particle object coordinates in internal representation = "; MatrixAlgebra::print(POcoord);
  std::cout << "using parametrisation representation: " << POtest.GetParametrisation()->name() << std::endl;

//***************************************************************************<*
// Complete WW kinematic fit example illustrating usage of constraints and fit
//***************************************************************************<*
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Complete WW kinematic fit example @ 183 GeV" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;

  std::vector<Coordinates> WWList = { //Particle coordinates in PxPyPzE representation for WW->12 34 particle event
    {2.67377, 13.5278, -14.7699, 21.5475}, //lepton
    {-36.0753, -33.4628, 32.5065, 64.3082}, //neutrino: determined by using 4-momentum conservation constraints
    {35.4921, -19.6482, -32.833, 52.1897}, //jet
    {-2.09058, 39.5834, 15.0964, 42.416}}; //jet
  /*C PARTICLE FIVE -- YES! WHY NOT? WHEN ABCFIT CAN DO IT!
C THIS COULD A THREE-JET LIKE W DECAY....
P_REC(1,5)=-2.13801       ! PX
P_REC(2,5)=18.2512        ! PY
P_REC(3,5)=7.74521        ! PZ
P_REC(4,5)=22.207         ! E   
  */

  std::vector<ParticleObject*> WWParticleList;
  for (int i=0; i<WWList.size(); i++) {
    ParticleObject* po = new ParticleObject(WWList[i], pxpypze, abcd); //NOTE: all particles are assumed measured including the neutrino
    WWParticleList.push_back(po); 
    po->FixParameter(3); //Fix particle mass parameter
  }
  //Setup of usual 4-momentum conservation constraints: SumPx, SumPy, SumPz, SumE
  //Following constraints define particle systems (all particles) to be used
  SumPxConstraint WWPxConstraint(WWParticleList); 
  SumPyConstraint WWPyConstraint(WWParticleList); 
  SumPzConstraint WWPzConstraint(WWParticleList); 
  SumEConstraint  WWEConstraint(WWParticleList); 
  //Define particle systems for later W mass constraints 
  InvMassConstraint WWInvMass1({WWParticleList[0],WWParticleList[1]});
  InvMassConstraint WWInvMass2({WWParticleList[2],WWParticleList[3]});
  //Define particle systems for later energy sum constraints
  SumEConstraint  WWEConstraint1({WWParticleList[0],WWParticleList[1]});
  SumEConstraint  WWEConstraint2({WWParticleList[2],WWParticleList[3]});

  //A fit is defined by: 
  //    -> a sequence of composite constraints
  //        -> linear combinations of constraints 
  //            -> list of particles
  FitKernel WWfit(std::vector<CompositeConstraint> {
      CompositeConstraint ({{&WWPxConstraint,1.0}}),   // this lines means "1.0*Sum(Px) = 0.0" (constraint value defaults to zero) 
	CompositeConstraint ({{&WWPyConstraint,1.0}}), // this lines means "1.0*Sum(Py) = 0.0" (constraint value defaults to zero) 
	CompositeConstraint ({{&WWPzConstraint,1.0}}), // this lines means "1.0*Sum(Pz) = 0.0" (constraint value defaults to zero) 
	CompositeConstraint ({{&WWEConstraint,1.0}},183.0) // this lines means "1.0*Sum(E) = 183.0"
	},10); //Number of iterations
  FitResult WWTest = WWfit.Fit();
  std:: cout << "4-momentum conservation fit output : " << std::endl;
  WWTest.print(true);
  std::cout << "---------------------------------------------" << std::endl;
  //Constraints can be applied in steps
  //Additional constraints can be added to the fit by a new fit using the SAME particles (automatically updated in previous fit)
  //This can be used to gradually increase constraint complexity in case of convergence issues
  //Adding equal mass constraint
  FitKernel WWfitMASS(std::vector<CompositeConstraint> {
      CompositeConstraint ({{&WWPxConstraint,1.0}}), 
	CompositeConstraint ({{&WWPyConstraint,1.0}}), 
	CompositeConstraint ({{&WWPzConstraint,1.0}}), 
	CompositeConstraint ({{&WWEConstraint,1.0}},183.0), 
	CompositeConstraint ({{&WWInvMass1,1.0},{&WWInvMass2,-1.0}}) //Additional equal mass constraint: "1.0*InvMass1 -1.0*InvMass2 = 0.0"
	},10);
  FitResult WWTestMASS = WWfitMASS.Fit();
  std:: cout << "Equal mass fit output : " << std::endl;
  WWTestMASS.print(true);
  std::cout << "Mass of W1 = " << WWInvMass1.ConstraintFunction() << std::endl;
  std::cout << "Mass of W2 = " << WWInvMass2.ConstraintFunction() << std::endl;
  std::cout << "---------------------------------------------" << std::endl;

  //Now trying equal energy – should give same result
  //Let's start from scratch: Reset particle objects
  for (int i=0; i<WWParticleList.size(); i++) {
    WWParticleList[i]->SetCoordinates(WWList[i]); //Setting coordinates in input representation
  }

  FitKernel WWfitENERGY(std::vector<CompositeConstraint> {
      CompositeConstraint ({{&WWPxConstraint,1.0}}), 
	CompositeConstraint ({{&WWPyConstraint,1.0}}), 
	CompositeConstraint ({{&WWPzConstraint,1.0}}), 
	CompositeConstraint ({{&WWEConstraint,1.0}},183.0), 
	CompositeConstraint ({{&WWEConstraint1,1.0},{&WWEConstraint2,-1.0}}) //Equal energy constraint: "1.0*Sum(E(set1)) -1.0*Sum(E(set2)) = 0.0" 
	},10);
  FitResult WWTestENERGY = WWfitENERGY.Fit();
  std:: cout << "Equal energy fit output : " << std::endl;
  WWTestENERGY.print(true);
  std::cout << "Mass of W1 = " << WWInvMass1.ConstraintFunction() << std::endl;
  std::cout << "Mass of W2 = " << WWInvMass2.ConstraintFunction() << std::endl;
  std::cout << "---------------------------------------------" << std::endl;

  //Now trying with exact W mass constraints (one for each W)
  FitKernel WWfitWMASS(std::vector<CompositeConstraint> {
      CompositeConstraint ({{&WWPxConstraint,1.0}}), 
	CompositeConstraint ({{&WWPyConstraint,1.0}}), 
	CompositeConstraint ({{&WWPzConstraint,1.0}}), 
	CompositeConstraint ({{&WWEConstraint,1.0}},183.0), 
	CompositeConstraint ({{&WWInvMass1,1.0}},80.5), // Exact mass constraint: "1.0*InvMass1 = 80.5"
	CompositeConstraint ({{&WWInvMass2,1.0}},80.5)  // Exact mass constraint: "1.0*InvMass2 = 80.5"
	},10);
  FitResult WWTestWMASS = WWfitWMASS.Fit();
  std:: cout << "Exact W mass fit output : " << std::endl;
  WWTestWMASS.print(true);
  std::cout << "---------------------------------------------" << std::endl;
  //Now trying with Gaussian distributed mass constraints (one for each W)

  GaussianPDF WmassDist(80.5,2.7); //Gaussian pdf with mean 80.5 and sigma 2.7

  FitKernel WWgausfit(std::vector<CompositeConstraint> {
      CompositeConstraint ({{&WWPxConstraint,1.0}}), 
	CompositeConstraint ({{&WWPyConstraint,1.0}}), 
	CompositeConstraint ({{&WWPzConstraint,1.0}}), 
	CompositeConstraint ({{&WWEConstraint,1.0}},183.0), 
	CompositeConstraint ({{&WWInvMass1,1.0}},&WmassDist), // Constraint value distributed around 80.5 with sigma = 2.7
	CompositeConstraint ({{&WWInvMass2,1.0}},&WmassDist)  // Constraint value distributed around 80.5 with sigma = 2.7
	},10);
  FitResult WWTestgaus = WWgausfit.Fit();
  std:: cout << "Gaussian distributed W mass fit output : " << std::endl;
  WWTestgaus.print(true);
  //Printout mass values of constraints
  std::cout << "Mass of W1 = " << WWInvMass1.ConstraintFunction() << std::endl;
  std::cout << "Mass of W2 = " << WWInvMass2.ConstraintFunction() << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "End of WW example" << std::endl;
//***************************************************************************<*
// Complete ttbar kinematic fit example illustrating usage of constraints and fit
//***************************************************************************<*
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Complete ttbar kinematic fit example @ 365 GeV" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;

  //This example illustrates custom parameterisation setup for chi2-calculation
  CoordRepr* abcdlepton = new ABCD(); //lepton
  CoordRepr* abcdbjet = new ABCD();   //b-jet
  CoordRepr* abcdlfjet = new ABCD();  //light flavour jet
  //Setting covariance matrices
  abcdlepton->SetParametrisationCovMatrixFunction(&LeptonParameterisationCovMatrix);
  abcdbjet->SetParametrisationCovMatrixFunction(&BjetParameterisationCovMatrix);
  abcdlfjet->SetParametrisationCovMatrixFunction(&LFjetParameterisationCovMatrix);
  //Setting parametrisation for jets only (lepton is default (1.0,0.0,0.0,mass))
  abcdbjet->SetParametrisationFunction(&JetParameterisation);
  abcdlfjet->SetParametrisationFunction(&JetParameterisation);

  //Event for semileptonic top: in this event, the neutrino is treated as an unmeasured particle, however an initial value 4-momenta must be provided 
  //Unmeasured momenta should be set by using some of the constraints
  std::vector<Coordinates> TopList = {//Particle coordinates in PxPyPzE representation for ttbar->123 456 particle event  
    {14.3186, -25.5041, -7.69136, 30.243}, //lepton
    {21.6725, 23.0272, 51.9512, 60.8184},  //neutrino: determined by using 4-momentum conservation constraints
    {8.1591, 57.7045, -69.235, 91.7773},   //b-jet from leptonic top 
    {-30.4919, -51.0642, -12.7182, 61.4842}, //b-jet from hadronic top
    {12.2401, -28.7724, 27.6045, 46.4874}, //lf jet
    {-25.8985, 24.609, 10.0887, 37.7875}   //lf jet
  };
  //Particle type dependent parametrisations
  std::vector<ParticleObject*> TopParticleList;
  TopParticleList.push_back(new ParticleObject(TopList[0], pxpypze, abcdlepton));
  TopParticleList.push_back(new ParticleObject(TopList[1], pxpypze)); //NOTE unmeasured particles do not contribute to chi2
  TopParticleList.push_back(new ParticleObject(TopList[2], pxpypze, abcdbjet));
  TopParticleList.push_back(new ParticleObject(TopList[3], pxpypze, abcdbjet));
  TopParticleList.push_back(new ParticleObject(TopList[4], pxpypze, abcdlfjet));
  TopParticleList.push_back(new ParticleObject(TopList[5], pxpypze, abcdlfjet));
  
  for (int i=0; i<TopList.size(); i++) {
    TopParticleList[i]->FixParameter(3); //Particles masses are kept fixed (in this example input particle masses are all zero)
  }

  std:: cout << "Lepton ABCD expectation values = " ; MatrixAlgebra::print(abcdlepton->GetExpectation(test));
  std:: cout << "Lepton ABCD Covariance Matrix = " ;MatrixAlgebra::print(abcdlepton->GetExpectationCovMatrix(test));


  //Now we are ready for the fit
  //Setup of usual 4-momentum conservation constraints: SumPx, SumPy, SumPz, SumE
  //Following constraints define particle systems (all particles) to be used
  SumPxConstraint TopPxConstraint(TopParticleList); 
  SumPyConstraint TopPyConstraint(TopParticleList); 
  SumPzConstraint TopPzConstraint(TopParticleList); 
  SumEConstraint  TopEConstraint(TopParticleList); 
  //Define particle systems for later top mass constraints 
  InvMassConstraint TopInvMass1({TopParticleList[0],TopParticleList[1],TopParticleList[2]});
  InvMassConstraint TopInvMass2({TopParticleList[3],TopParticleList[4],TopParticleList[5]});
  //Define particle systems for later W mass constraints 
  InvMassConstraint WInvMass1({TopParticleList[0],TopParticleList[1]});
  InvMassConstraint WInvMass2({TopParticleList[4],TopParticleList[5]});

  GaussianPDF TopWmassDist(80.5,2.7); //Gaussian pdf with mean 80.5 and sigma 2.7  
  GaussianPDF TopmassLepDist(173.5,10.0); //Example Gaussian resolution taken from FCC-ee simulation
  GaussianPDF TopmassHadDist(173.5,9.0);  //Example Gaussian resolution taken from FCC-ee simulation

  FitKernel TopfitMASS(std::vector<CompositeConstraint> {CompositeConstraint ({{&TopPxConstraint,1.0}}), 
							 CompositeConstraint ({{&TopPyConstraint,1.0}}),
							 CompositeConstraint ({{&TopPzConstraint,1.0}}),
							 CompositeConstraint ({{&TopEConstraint,1.0}},365.0),
							 CompositeConstraint ({{&TopInvMass1,1.0}},&TopmassLepDist),
							 CompositeConstraint ({{&TopInvMass2,1.0}},&TopmassHadDist),
							 CompositeConstraint ({{&WInvMass1,1.0}},&TopWmassDist),
                                                         CompositeConstraint ({{&WInvMass2,1.0}},&TopWmassDist)
    }); 
  FitResult TopTestMASS = TopfitMASS.Fit();
  std:: cout << "Gaussian distributed Top and W masses fit output : " << std::endl;
  TopTestMASS.print(true);
  std::cout << "Mass of leptonic W = " << WInvMass1.ConstraintFunction() << std::endl;
  std::cout << "Mass of hadronic W = " << WInvMass2.ConstraintFunction() << std::endl;
  std::cout << "Mass of leptonic top = " << TopInvMass1.ConstraintFunction() << std::endl;
  std::cout << "Mass of hadronic top = " << TopInvMass2.ConstraintFunction() << std::endl;
  std::cout << "---------------------------------------------" << std::endl;  

  //Following lines illustrates how to access covariance matrix of parameters or particle momenta
  std::cout << "Following lines shows covariance matrix of parameters: "<< std::endl;
  std::cout << " ABCD cov matrix = "; MatrixAlgebra::print(TopfitMASS.GetParameterCovMatrix(),true);
  std::cout << " and particle momenta: "<< std::endl;
  std::cout << " Momenta cov matrix = "; MatrixAlgebra::print(TopfitMASS.GetParticleCovMatrix(),true);


  // clean up memory
  for(auto ptr : WWParticleList) delete ptr;
  for(auto ptr : TopParticleList) delete ptr;
  return 0;
}
