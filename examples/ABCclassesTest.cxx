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
#include "ABC.h"
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

//assumes Coordinate Representation in internal representation                                                                                                    
Matrix TestParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  double sigmaE = 0.15*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]); //sampling term 15%                                                                         
  result[0][0] = 10*sigmaE;
  result[1][1] = 0.3;
  result[2][2] = 0.3; //Lambda_QCD quark pT                                                                                                                       
  result[3][3] = 20*sigmaE;
  return result;
}
Matrix LeptonParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  result[0][0] = 0.003;
  result[1][1] = 0.003;
  result[2][2] = 0.003; 
  result[3][3] = 10.0; //mass value not used  
  return result;
}
Matrix BjetParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  result[0][0] = 0.07;
  result[1][1] = 1.2;
  result[2][2] = 1.2; 
  result[3][3] = 10.0; //mass value not used
  return result;
}
Matrix LFjetParameterisationCovMatrix(Coordinates p) {
  Matrix result(4,{0.0,0.0,0.0,0.0});
  result[0][0] = 0.07;
  result[1][1] = 1.4;
  result[2][2] = 1.4; 
  result[3][3] = 10.0; //mass value not used
  return result;
}
Coordinates JetParameterisation(Coordinates p) {
  Coordinates result({0.0,0.0,0.0,0.0});
  result[0] = 1.03;
  result[1] = 0.0;
  result[2] = 0.0; 
  result[3] = 10.0; //mass value not used
  return result;
}

int main() {
  //Coordinates for test particle
  Coordinates test{1.0,2.0,3.0,4.0};

  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "CoordRepr with PxPyPzM, PxPyPzE and PtEtaPhiM" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;

  //Testing for the internal class PxPyPzM
  Coordinates outPxPyPzM1 = PxPyPzM().Transform(test);
  Coordinates outPxPyPzM2 = PxPyPzM().InvTransform(test);

  std::cout << "--------- Internal class PxPyPzM ---------:" << std::endl;
  std::cout << "Transform: ";
  for (auto x_i :outPxPyPzM1) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "InvTransform: "; 
  for (auto x_i :outPxPyPzM2) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "Derivative: " << std::endl;
  for (auto row : PxPyPzM().Derivative(test)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << "InvDerivative: " << std::endl;
  for (auto row : PxPyPzM().InvDerivative(test)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }

  //Testing for class PxPyPzE
  Coordinates outPxPyPzE1 = PxPyPzE().Transform(test);
  Coordinates outPxPyPzE2 = PxPyPzE().InvTransform(test);
  std::cout << "--------- PxPyPzE class ---------:" << std::endl;
  std::cout << "Transform: ";
  for (auto x_i :outPxPyPzE1) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "InvTransform: ";
  for (auto x_i :outPxPyPzE2) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "Derivative: " << std::endl;
  for (auto row : PxPyPzE().Derivative(test)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << "InvDerivative: " << std::endl;
  for (auto row : PxPyPzE().InvDerivative(test)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }  

  //Testing for class PtEtaPhiM
  Coordinates outPtEtaPhiM1 = PtEtaPhiM().Transform(test);
  Coordinates outPtEtaPhiM2 = PtEtaPhiM().InvTransform(test);
  std::cout << "--------- PtEtaPhiM class ---------:" << std::endl;
  std::cout << "Transform: ";
  for (auto x_i :outPtEtaPhiM1) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "InvTransform: ";
  for (auto x_i :outPtEtaPhiM2) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "Derivative: " << std::endl;
  for (auto row : PtEtaPhiM().Derivative(test)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << "InvDerivative: " << std::endl;
  for (auto row : PtEtaPhiM().InvDerivative(test)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }

  std::cout << "---- Testing transform one PxPyPzE to PtEtaPhiM via internal representation----: " << std::endl;

  PxPyPzE repPxPyPzE;
  PtEtaPhiM repPtEtaPhiM;
  Coordinates out2transforms = ((CoordRepr*)&repPxPyPzE)->Transform(test,(CoordRepr*)&repPtEtaPhiM);
  for (auto x_i :out2transforms) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "MatrixAlgebra" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  //Coordinates for test vectors
  Coordinates V1{1.0,2.0,3.0,4.0};
  Coordinates V2{2.0,3.0,4.0,5.0};
  std::vector<Coordinates> M1{
    { 1.0, 2.0, 3.0, 4.0 },
    { 2.0, 4.0, 6.0, 8.0 },
    { 3.0, 6.0, 9.0, 12.0 },
    { 4.0, 8.0, 12.0, 16.0 },
  };

  std::vector<Coordinates> M2{
    { 1.0, 2.0, 3.0, 4.0 },
    { 0.0, 1.0, 2.0, 3.0 },
    { 0.0, 0.0, 1.0, 2.0 },
    { 0.0, 0.0, 0.0, 1.0 },
  };
  std::cout << "--------- Vector sum & diff ---------:" << std::endl;
  for (auto x_i :V1) std::cout << x_i << ", ";
  std::cout << " + ";
  for (auto x_i :V2) std::cout << x_i << ", ";
  std::cout << " = ";
  for (auto x_i : MatrixAlgebra::VectorSum(V1,V2,1.0)) std::cout << x_i << ", ";
  std::cout << std::endl;

  for (auto x_i :V1) std::cout << x_i << ", ";
  std::cout << " - ";
  for (auto x_i :V2) std::cout << x_i << ", ";
  std::cout << " = ";
  for (auto x_i : MatrixAlgebra::VectorSum(V1,V2,-1.0)) std::cout << x_i << ", ";
  std::cout << std::endl;

  std::cout << "--------- Matrix sum & diff ---------:" << std::endl;
  for (auto row : M1){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " + "<< std::endl;
  for (auto row : M2){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " = " << std::endl;
  for (auto row :  MatrixAlgebra::MatrixSum(M1,M2,1.0)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  for (auto row : M1){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " - "<< std::endl;
  for (auto row : M2){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " = " << std::endl;
  for (auto row :  MatrixAlgebra::MatrixSum(M1,M2,-1.0)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }

  std::cout << "--------- TwoVector Multiplication ---------:" << std::endl;
  for (auto x_i :V1) std::cout << x_i << ", ";
  std::cout << " * ";
  for (auto x_i :V2) std::cout << x_i << ", ";
  std::cout << " = " << MatrixAlgebra::VectorMultiplication(V1,V2) << std::endl;
  V2.erase(V2.begin());
  for (auto x_i :V1) std::cout << x_i << ", ";
  std::cout << " * ";
  for (auto x_i :V2) std::cout << x_i << ", ";
  std::cout << " = " << MatrixAlgebra::VectorMultiplication(V1,V2) << std::endl;

  std::cout << "--------- MatrixVector Multiplication ---------:" << std::endl;
  for (auto row : M2){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " * ";
  for (auto x_i :V1) std::cout << x_i << ", ";
  std::cout << " = ";
  for (auto x_i : MatrixAlgebra::MatrixVectorMultiplication(M2,V1)) std::cout << x_i << ", ";
  std::cout << std::endl;

  for (auto row : M1){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " * ";
  for (auto x_i :V2) std::cout << x_i << ", ";
  std::cout << " = ";
  for (auto x_i : MatrixAlgebra::MatrixVectorMultiplication(M1,V2)) std::cout << x_i << ", ";
  std::cout << std::endl;

  std::cout << "--------- TwoMatrix Multiplication ---------:" << std::endl;
  for (auto row : M1){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " * "<< std::endl;
  for (auto row : M2){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " = " << std::endl;
  for (auto row :  MatrixAlgebra::MatrixMultiplication(M1,M2)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::vector<Coordinates> M3{
    { 1.0, 2.0, 3.0 },
      { 4.0, 5.0, 6.0 },
	{ 7.0, 8.0, 9.0 },
	  { 10.0, 11.0, 12.0 },
	    };
  M2.erase(M2.begin());
  for (auto row : M1){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " * "<< std::endl;
  for (auto row : M2){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << " = "<< std::endl;
  for (auto row :  MatrixAlgebra::MatrixMultiplication(M1,M2)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }

  std::cout << "--------- Inverse Matrix ---------:" << std::endl;
  std::vector<Coordinates> M4{
    { 1.0, 3.0, 2.0 },
      { 4.0, 1.0, 3.0 },
        { 2.0, 5.0, 2.0 }
            };

  std::cout << " Inverse of "<< std::endl;
  for (auto row : M4){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << "with derminant = " << MatrixAlgebra::Determinant(&M4) << " is " << std::endl;
  for (auto row :  MatrixAlgebra::InverseMatrix(&M4)){
    for (auto col : row)
      std::cout << col << ", ";
    std::cout << std::endl;
  }
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "ProbDistFunc" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "ParticleObject" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  ParticleObject POtest(test, pxpypze, ptetaphim);
  Coordinates POcoord=POtest.GetCoordinates();
  for (auto x_i : POcoord) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << POtest.GetParametrisation()->name() << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "constraint" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::vector<Coordinates> List = {
    {1.0,2.0,3.0,4.0},
    {1.0,-1.0,2.0,4.0},
    {0.0,2.0,0.0,3.0}
  };

  std::vector<Coordinates> abcfitList = {
    {2.67377    ,
     13.5278    ,
     -14.7699   ,
     21.5475    },
    {-36.0753   ,
     -33.4628   ,
     32.5065    ,
     64.3082}    ,
    {35.4921    ,
     -19.6482   ,
     -32.833    ,
     52.1897}    ,
    {-2.09058   ,
     39.5834    ,
     15.0964    ,
     42.416}};
  /*C PARTICLE FIVE -- YES! WHY NOT? WHEN ABCFIT CAN DO IT!                                                                                                                                                                                        
C THIS COULD A THREE-JET LIKE W DECAY....                                                                                                                                                                                                      
P_REC(1,5)=-2.13801       ! PX                                                                                                                                                                                                           
P_REC(2,5)=18.2512        ! PY                                                                                                                                                                                                           
P_REC(3,5)=7.74521        ! PZ                                                                                                                                                                                                           
P_REC(4,5)=22.207         ! E   
  */

  std::vector<Coordinates> TopList = {
    {14.3186,
     -25.5041,
     -7.69136,
     30.243},
    {21.6725,
     23.0272,
     51.9512,
     60.8184},
    {8.1591,
     57.7045,
     -69.235,
     91.7773},
    {-30.4919,
     -51.0642,
     -12.7182,
     61.4842},
    {12.2401,
     -28.7724,
     27.6045,
     46.4874},
    {-25.8985,
     24.609,
     10.0887,
     37.7875}
  };
  std::vector<ParticleObject*> ParticleList;
  for (int i=0; i<List.size(); i++) {
    ParticleObject* po = new ParticleObject(List[i], pxpypze, abcd); //saves ParticleObject (object) to the heap
    ParticleList.push_back(po); //pushes bak pointer to find ParticleObjects again on the heap
  }
  std::vector<ParticleObject*> abcfitParticleList;
  for (int i=0; i<abcfitList.size(); i++) {
    ParticleObject* po = new ParticleObject(abcfitList[i], pxpypze, abcd); //saves ParticleObject (object) to the heap
    abcfitParticleList.push_back(po); //pushes bak pointer to find ParticleObjects again on the heap
    po->FixParameter(3);
  }
  CoordRepr* abcdlepton = new ABCD();
  CoordRepr* abcdbjet = new ABCD();
  CoordRepr* abcdlfjet = new ABCD();
  abcdlepton->SetParametrisationCovMatrixFunction(&LeptonParameterisationCovMatrix);
  abcdbjet->SetParametrisationCovMatrixFunction(&BjetParameterisationCovMatrix);
  abcdlfjet->SetParametrisationCovMatrixFunction(&LFjetParameterisationCovMatrix);
  abcdbjet->SetParametrisationFunction(&JetParameterisation);
  abcdlfjet->SetParametrisationFunction(&JetParameterisation);
  std::vector<ParticleObject*> TopParticleList;
  TopParticleList.push_back(new ParticleObject(TopList[0], pxpypze, abcdlepton));
  TopParticleList.push_back(new ParticleObject(TopList[1], pxpypze));
  TopParticleList.push_back(new ParticleObject(TopList[2], pxpypze, abcdbjet));
  TopParticleList.push_back(new ParticleObject(TopList[3], pxpypze, abcdbjet));
  TopParticleList.push_back(new ParticleObject(TopList[4], pxpypze, abcdlfjet));
  TopParticleList.push_back(new ParticleObject(TopList[5], pxpypze, abcdlfjet));
  
  for (int i=0; i<TopList.size(); i++) {
    TopParticleList[i]->FixParameter(3);    
  }
  
  Coordinates ABCvec = abcd->InvTransform(ParticleList[0]->GetCoordinates()); 
  std::cout << "ABCvec = ";
  MatrixAlgebra::print(ABCvec);

  std::cout << pxpypze << " " << pxpypzm << " " << abcd << std::endl;

  SumPxConstraint Px1(ListOfParticleObjects({ParticleList[0]}),0.0);
  SumPxConstraint Px2(ListOfParticleObjects({ParticleList[1]}),0.0);
  SumPxConstraint Px3(ListOfParticleObjects({ParticleList[2]}),0.0);
  CompositeConstraint SumPx({{&Px1,1.0},{&Px2,1.0},{&Px3,-2.0}});
  double SumPxFunc = SumPx.ConstraintFunction();
  Coordinates SumPxDeriv = SumPx.ConstraintFunctionDerivative();
  std::cout << "Composite constraint : " << "Sum = " << SumPxFunc << " Deriv = ";
  MatrixAlgebra::print(SumPxDeriv);

  InvMassConstraint InvMass1Particle(ListOfParticleObjects({ParticleList[0],ParticleList[1]}),0.0);  
  InvMassConstraint InvMass2Particle(ListOfParticleObjects({ParticleList[1],ParticleList[2]}),0.0);
  std::cout << "-------------------------------------------" << std::endl;
  MatrixAlgebra::print(InvMass1Particle.GetListOfParticleObjects());
  MatrixAlgebra::print(InvMass2Particle.GetListOfParticleObjects());
  std::cout << "-------------------------------------------" <<std::endl;
  CompositeConstraint EqualMass({{&InvMass1Particle,1.0}, {&InvMass2Particle,-1.0}});
  std::cout << "blah : " << std::endl;
  double SingleMFunc = InvMass1Particle.ConstraintFunction();
  std::cout << "blah : " << std::endl;
  Coordinates SingleMDeriv = InvMass1Particle.ConstraintFunctionDerivative();
  std::cout << "1 particle mass constraint : " << "Sum = " << SingleMFunc << " Deriv = ";
  MatrixAlgebra::print(SingleMDeriv);

  double EqualMassFunc = EqualMass.ConstraintFunction();
  Coordinates EqualMassDeriv = EqualMass.ConstraintFunctionDerivative();
  std::cout << "Equal mass constraint : " << "Sum = " << EqualMassFunc << " Deriv = ";
  MatrixAlgebra::print(EqualMassDeriv);


  SumPxConstraint PxConstraint(ParticleList,0.0); 
  SumPyConstraint PyConstraint(ParticleList,0.0); 
  SumPzConstraint PzConstraint(ParticleList,0.0); 
  SumEConstraint  EConstraint(ParticleList,0.0); 
  InvMassConstraint InvMConstraint(ParticleList,0.0);
  double PxFunc = PxConstraint.ConstraintFunction();
  double PyFunc = PyConstraint.ConstraintFunction();
  double PzFunc = PzConstraint.ConstraintFunction();
  double EFunc  = EConstraint.ConstraintFunction();
  double InvMFunc = InvMConstraint.ConstraintFunction();
  Coordinates PxDeriv = PxConstraint.ConstraintFunctionDerivative();
  Coordinates PyDeriv = PyConstraint.ConstraintFunctionDerivative();
  Coordinates PzDeriv = PzConstraint.ConstraintFunctionDerivative();
  Coordinates EDeriv  = EConstraint.ConstraintFunctionDerivative();
  Coordinates InvMDeriv = InvMConstraint.ConstraintFunctionDerivative();
  Coordinates InternalDeriv = InvMConstraint.ConstraintInternalDerivative();
  /*
  int Nparticles = InvMConstraint.GetParticleCoordinates().size();
  Matrix dpdy(4*Nparticles,Coordinates (4*Nparticles,0.0));
 
  for (unsigned int i=0; i<InvMConstraint.GetParticleCoordinates().size(); i++) {
    Coordinates p = InvMConstraint.GetParticleCoordinates()[i];
    //std::cout << "ladida " << p.GetCoordinates().size() << std::endl;
    Matrix matrix = p.GetParametrisation()->Derivative(p.GetParametrisation()->InvTransform(p));
    for (unsigned int j=0; j<4; j++) std::copy(matrix[j].begin(),matrix[j].end(),dpdy[j+i*4].begin()+4*i);
  }
  std::cout << "dpdy = " << std::endl;
  MatrixAlgebra::print(dpdy);
  std::cout << "Bmatrix = " << std::endl;
  MatrixAlgebra::print(MatrixAlgebra::VectorMatrixMultiplication(InternalDeriv,dpdy)); //dfdy=dfdp*dpdy
  */

  std::cout << "With list of particles: " << std::endl;
  for (auto X :  List){
    for (auto x_i : X)
      std::cout << x_i << ", ";
    std::cout << std::endl;
  }  

  std::cout << "Constraint Functions are :" << " px = " << PxFunc << " py = " << PyFunc << " pz = " << PzFunc << " E = " << EFunc << " InvMass = " << InvMFunc << std::endl;
  std::cout << "Constraint Functions derivatives are :" << std::endl;
  std::cout << "px = ";
  for (auto x_i: PxDeriv) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "py = ";
  for (auto x_i: PyDeriv) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "pz = ";
  for (auto x_i: PzDeriv) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "E = ";
  for (auto x_i: EDeriv) std::cout << x_i << ", ";
  std::cout << std::endl;
  std::cout << "InvMass = ";
  MatrixAlgebra::print(InvMDeriv);

  std::cout << "Internal E = ";
  for (auto x_i: InternalDeriv) std::cout << x_i << ", ";
  std::cout << std::endl;

  std:: cout << "Test parametrisation : " ;
  MatrixAlgebra::print(abcd->GetExpectation(test));
  MatrixAlgebra::print(abcd->GetExpectationCovMatrix(test));
  //abcd->SetParametrisationCovMatrixFunction(&TestParameterisationCovMatrix);
  MatrixAlgebra::print(abcd->GetExpectationCovMatrix(test));


  FitKernel myFit(std::vector<CompositeConstraint> {SumPx,EqualMass},10);
  FitResult FitTest = myFit.Fit();
  FitTest.print();

  std:: cout << "TestABCTest : " << std::endl;

  SumPxConstraint abcfitPxConstraint(abcfitParticleList,0.0); 
  SumPyConstraint abcfitPyConstraint(abcfitParticleList,0.0); 
  SumPzConstraint abcfitPzConstraint(abcfitParticleList,0.0); 
  SumEConstraint  abcfitEConstraint(abcfitParticleList,183.0); 
  InvMassConstraint abcfitInvMass1({abcfitParticleList[0],abcfitParticleList[1]});
  InvMassConstraint abcfitInvMass2({abcfitParticleList[2],abcfitParticleList[3]});
  SumEConstraint  abcfitEConstraint1({abcfitParticleList[0],abcfitParticleList[1]});
  SumEConstraint  abcfitEConstraint2({abcfitParticleList[2],abcfitParticleList[3]});

  //std::cout << "ladida 64738 " << abcfitEConstraint.GetConstraintValue() << std::endl;
  //CompositeConstraint EqualMass({{&InvMass1Particle,1.0}, {&InvMass2Particle,-1.0}})
  /*FitKernel abcfitfit(std::vector<CompositeConstraint> {CompositeConstraint ({{&abcfitPxConstraint,1.0}}), CompositeConstraint ({{&abcfitPyConstraint,1.0}}), CompositeConstraint ({{&abcfitPzConstraint,1.0}}), CompositeConstraint ({{&abcfitEConstraint,1.0}},183.0)},10);
  FitResult abcfitTest = abcfitfit.Fit();
  abcfitTest.print(true);
  
  std:: cout << "Equal mass : " << std::endl;

  FitKernel abcfitfitMASS(std::vector<CompositeConstraint> {CompositeConstraint ({{&abcfitPxConstraint,1.0}}), CompositeConstraint ({{&abcfitPyConstraint,1.0}}), CompositeConstraint ({{&abcfitPzConstraint,1.0}}), CompositeConstraint ({{&abcfitEConstraint,1.0}},183.0), CompositeConstraint ({{&abcfitInvMass1,1.0},{&abcfitInvMass2,-1.0}})},10);
  FitResult abcfitTestMASS = abcfitfitMASS.Fit();
  abcfitTestMASS.print(true);
  */
  std:: cout << "W mass : " << std::endl;
  GaussianPDF WmassDist(80.5,2.7);

  FitKernel abcfitfitWMASS(std::vector<CompositeConstraint> {CompositeConstraint ({{&abcfitPxConstraint,1.0}}), CompositeConstraint ({{&abcfitPyConstraint,1.0}}), CompositeConstraint ({{&abcfitPzConstraint,1.0}}), CompositeConstraint ({{&abcfitEConstraint,1.0}},183.0), CompositeConstraint ({{&abcfitInvMass1,1.0}},80.5), CompositeConstraint ({{&abcfitInvMass2,1.0}},80.5)},10);
  //FitKernel abcfitfitWMASS(std::vector<CompositeConstraint> {CompositeConstraint ({{&abcfitPxConstraint,1.0}}), CompositeConstraint ({{&abcfitPyConstraint,1.0}}), CompositeConstraint ({{&abcfitPzConstraint,1.0}}), CompositeConstraint ({{&abcfitEConstraint,1.0}},183.0), CompositeConstraint ({{&abcfitInvMass1,1.0}},&WmassDist), CompositeConstraint ({{&abcfitInvMass2,1.0}},&WmassDist)},10);
  FitResult abcfitTestWMASS = abcfitfitWMASS.Fit();
  abcfitTestWMASS.print(true);
    
  std:: cout << "Equal energy : " << std::endl;

  FitKernel abcfitfitENERGY(std::vector<CompositeConstraint> {CompositeConstraint ({{&abcfitPxConstraint,1.0}}), CompositeConstraint ({{&abcfitPyConstraint,1.0}}), CompositeConstraint ({{&abcfitPzConstraint,1.0}}), CompositeConstraint ({{&abcfitEConstraint,1.0}},183.0), CompositeConstraint ({{&abcfitEConstraint1,1.0},{&abcfitEConstraint2,-1.0}})},10);
  FitResult abcfitTestENERGY = abcfitfitENERGY.Fit();
  abcfitTestENERGY.print(true);

  std::cout << "-------------- TOP -----------" << std::endl;
  SumPxConstraint TopPxConstraint(TopParticleList,0.0); 
  SumPyConstraint TopPyConstraint(TopParticleList,0.0); 
  SumPzConstraint TopPzConstraint(TopParticleList,0.0); 
  SumEConstraint  TopEConstraint(TopParticleList,183.0); 
  InvMassConstraint TopInvMass1({TopParticleList[0],TopParticleList[1],TopParticleList[2]});
  InvMassConstraint TopInvMass2({TopParticleList[3],TopParticleList[4],TopParticleList[5]});
  InvMassConstraint WInvMass1({TopParticleList[0],TopParticleList[1]});
  InvMassConstraint WInvMass2({TopParticleList[4],TopParticleList[5]});
  SumEConstraint  TopEConstraint1({TopParticleList[0],TopParticleList[1],TopParticleList[2]});
  SumEConstraint  TopEConstraint2({TopParticleList[3],TopParticleList[4],TopParticleList[5]});

  GaussianPDF TopmassLepDist(173.5,10.0);
  GaussianPDF TopmassHadDist(173.5,20.0);
  std::cout << "Ladida 48.5" << std::endl;
  FitKernel TopfitMASS(std::vector<CompositeConstraint> {CompositeConstraint ({{&TopPxConstraint,1.0}}), 
							 CompositeConstraint ({{&TopPyConstraint,1.0}}),
							 CompositeConstraint ({{&TopPzConstraint,1.0}}),
							 CompositeConstraint ({{&TopEConstraint,1.0}},365.0),
							 CompositeConstraint ({{&TopInvMass1,1.0}},&TopmassLepDist),
							 CompositeConstraint ({{&TopInvMass2,1.0}},&TopmassHadDist),
							 CompositeConstraint ({{&WInvMass1,1.0}},&WmassDist),
                                                         CompositeConstraint ({{&WInvMass2,1.0}},&WmassDist)
	},10); 
  std::cout << "Ladida 49" << std::endl;
  FitResult TopTestMASS = TopfitMASS.Fit();
  TopTestMASS.print(true);
  //std::cout << " ABCD cov matrix = "; MatrixAlgebra::print(TopfitMASS.GetParameterCovMatrix(),true);
  //std::cout << " Momenta cov matrix = "; MatrixAlgebra::print(TopfitMASS.GetParticleCovMatrix(),true);


  // clean up memory
  for(auto ptr : ParticleList) delete ptr;
  return 0;
}
