/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.6 $
// $Date: 2007-07-27 19:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/BeamColumnJoint2dThermal.cpp,v $
                                                                        
// Written: NM (nmitra@u.washington.edu)
// Created: April 2002
// Last Revised: January 2007
//
// Description: This file contains the class implementation for beam-column joint.
// This element is a 4 noded 12 dof (3 dof at each node) finite area super-element.
// The element takes in 13 different material types in order to simulate
// the inelastic action observed in a reinforced beam column joint.
// Details about this element can be found in 
// Mitra and Lowes, J. Struct. Eng. ASCE, Volume 133, Number 1 (January 2007), pp. 105-120
// Updates: Several concerning Joint formulation (presently a revised formulation for joints)


#include <BeamColumnJoint2dThermal.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <MatrixUtil.h>

#include <UniaxialMaterial.h>
#include <string.h>
#include <math.h>
#include <ElementResponse.h>
#include <elementAPI.h>
#include <string>

//static Matrix CoupledZeroLengthM6(6, 6);   // class wide matrix for 6*6
static Matrix TwoNodeM6;   // class wide matrix for 6*6
static Vector TwoNodeV6;   // class wide Vector for size 6

void* OPS_BeamColumnJoint2dThermal()
{
    if (OPS_GetNumRemainingInputArgs() < 6) {
	opserr << "WARNING insufficient arguments\n";
	  opserr << "Want: element beamColumnJoint eleTag? node1? node2? \n";
	  opserr << "matTag1? matTag2? matTag3?\n";
	  //opserr << "<ElementHeightFactor? ElementWidthFactor?>\n";
	  return 0;
    }

    int idata[6];
    int numdata = 6;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    double data[2] = {1.0, 1.0};
    numdata = 2;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
	    opserr<<"WARNING: invalid double inputs\n";
	    return 0;
	}
    }

    UniaxialMaterial* mats[3];
    for (int i = 0; i < 3; i++) {
	mats[i] = OPS_getUniaxialMaterial(idata[3+i]);
	if (mats[i] == 0) {
	    opserr<<"WARNING: material "<<idata[3+i]<<" is not defined\n";
	    return 0;
	}
    }

    return new BeamColumnJoint2dThermal(idata[0],idata[1],idata[2],
				 *mats[0],*mats[1],*mats[2],data[0],data[1]);

}

// full constructors:
BeamColumnJoint2dThermal::BeamColumnJoint2dThermal(int tag,int Nd1, int Nd2,
				     UniaxialMaterial& theMat1,
				     UniaxialMaterial& theMat2,
				     UniaxialMaterial& theMat3):
  Element(tag,ELE_TAG_BeamColumnJoint2d), connectedExternalNodes(2),
  nodeDbTag(0), dofDbTag(0), theMatrix(0),ub(0), ubdot(0), qb(0), ul(0), Tgl(0, 0), Tlb(0, 0), theVector(0), theLoad(0)
{
	// ensure the connectedExternalNode ID is of correct size & set values
 
	if (connectedExternalNodes.Size() != 2)
      opserr << "ERROR : BeamColumnJoint::BeamColumnJoint " << tag << "failed to create an ID of size 4" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;


	MaterialPtr = new UniaxialMaterial*[3];

	for (int x = 0; x <3; x++)
	{	MaterialPtr[x] = 0; }



	nodePtr[0] = 0;
	nodePtr[1] = 0;

// get a copy of the material and check we obtained a valid copy
  MaterialPtr[0] = theMat1.getCopy();
  if (!MaterialPtr[0]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 1" << endln;}
  MaterialPtr[1] = theMat2.getCopy();
  if (!MaterialPtr[1]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 2"<< endln;}
  MaterialPtr[2] = theMat3.getCopy();
  if (!MaterialPtr[2]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 3"<< endln;}
  

}

// full constructors:
BeamColumnJoint2dThermal::BeamColumnJoint2dThermal(int tag,int Nd1, int Nd2, 
				     UniaxialMaterial& theMat1,
				     UniaxialMaterial& theMat2,
				     UniaxialMaterial& theMat3,
					 double elHgtFac, double elWdtFac):
  Element(tag,ELE_TAG_BeamColumnJoint2d), connectedExternalNodes(4),
  nodeDbTag(0), dofDbTag(0), ub(0), ubdot(0), qb(0), ul(0),
	Tgl(0, 0), Tlb(0, 0), theMatrix(0), theVector(0), theLoad(0)
{
	// ensure the connectedExternalNode ID is of correct size & set values
 
	if (connectedExternalNodes.Size() != 2)
      opserr << "ERROR : BeamColumnJoint::BeamColumnJoint " << tag << "failed to create an ID of size 4" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;
    

	MaterialPtr = new UniaxialMaterial*[3];

	for (int x = 0; x <3; x++)
	{	MaterialPtr[x] = 0; }


	// added
	nodePtr[0] = 0;
	nodePtr[1] = 0;


	// get a copy of the material and check we obtained a valid copy
  MaterialPtr[0] = theMat1.getCopy();
  if (!MaterialPtr[0]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 1" << endln;}
  MaterialPtr[1] = theMat2.getCopy();
  if (!MaterialPtr[1]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 2"<< endln;}
  MaterialPtr[2] = theMat3.getCopy();
  if (!MaterialPtr[2]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 3"<< endln;}

}

// default constructor:
BeamColumnJoint2dThermal::BeamColumnJoint2dThermal():
  Element(0,ELE_TAG_BeamColumnJoint2d), connectedExternalNodes(2),
  nodeDbTag(0), dofDbTag(0), ub(0), ubdot(0), qb(0), ul(0),
	Tgl(0, 0), Tlb(0, 0), theMatrix(0), theVector(0), theLoad(0)
{
	nodePtr[0] = 0;
	nodePtr[1] = 0;
	for (int x = 0; x <3; x++)
	{	MaterialPtr[x] = 0; }
   // does nothing (invoked by FEM_ObjectBroker)
}

//  destructor:
BeamColumnJoint2dThermal::~BeamColumnJoint2dThermal()
{
	for (int i =0; i<3; i++)
	{
	    if (MaterialPtr[i] != 0)
		delete MaterialPtr[i];
	}

	if (MaterialPtr)
		 delete [] MaterialPtr;
}

// public methods
int
BeamColumnJoint2dThermal::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
BeamColumnJoint2dThermal::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **BeamColumnJoint2dThermal::getNodePtrs(void)
{
	return nodePtr;
}

int
BeamColumnJoint2dThermal::getNumDOF(void) 
{
    return 6;
}

void
BeamColumnJoint2dThermal::setDomain(Domain *theDomain)
{
    if (theDomain == 0)
	opserr << "ERROR : BeamColumnJoint::setDomain -- Domain is null" << endln;
	
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nodePtr[0] = 0;
	nodePtr[1] = 0;
    }

	//node pointers set
    for (int i = 0; i < 2; i++ ) {
		nodePtr[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
		if (nodePtr[i] == 0)
		{
			opserr << "ERROR : BeamColumnJoint::setDomain -- node pointer is null"<< endln;
			exit(-1); // donot go any further - otherwise segmentation fault
		}
	}

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

	// ensure connected nodes have correct dof's
	int dofNd1 = nodePtr[0]->getNumberDOF();
	int dofNd2 = nodePtr[1]->getNumberDOF();


	if ((dofNd1 != 3) || (dofNd2 != 3) ) 
	{
			opserr << "ERROR : BeamColumnJoint::setDomain -- number of DOF associated with the node incorrect"<< endln;
			exit(-1); // donot go any further - otherwise segmentation fault
	}


    // obtain the nodal coordinates    
	const Vector &end1Crd = nodePtr[0]->getCrds();
	const Vector &end2Crd = nodePtr[1]->getCrds();


	Vector Node1(end1Crd);
	Vector Node2(end2Crd);

	theMatrix = &TwoNodeM6;
	theVector = &TwoNodeV6;

}   

int
BeamColumnJoint2dThermal::commitState(void)
{

	int errCode = 0;

	// commit material models
	for (int i = 0; i < 3; i++)
		errCode += MaterialPtr[i]->commitState();

	// commit the base class
	errCode += this->Element::commitState();

	return errCode;
}

int
BeamColumnJoint2dThermal::revertToLastCommit(void)
{
	int errCode = 0;
	// revert material models
	for (int i = 0; i < 3; i++)
		errCode += MaterialPtr[i]->revertToLastCommit();

	return errCode;
}

int
BeamColumnJoint2dThermal::revertToStart(void)
{
	int mcs = 0;
	for (int j=0; j<3; j++)
	{
		if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->revertToStart();
		if (mcs != 0) break;
	}
	
	return mcs;
}

int
BeamColumnJoint2dThermal::update(void)
{	
	int errCode = 0;

	// get global trial response
	const Vector& dsp1 = nodePtr[0]->getTrialDisp();
	const Vector& dsp2 = nodePtr[1]->getTrialDisp();
	const Vector& vel1 = nodePtr[0]->getTrialVel();
	const Vector& vel2 = nodePtr[1]->getTrialVel();

	int numDOF = 6;
	int numDOF2 = 3;
	Vector ug(numDOF), ugdot(numDOF), uldot(numDOF);
	for (int i = 0; i < numDOF2; i++) {
		ug(i) = dsp1(i);  ugdot(i) = vel1(i);
		ug(i + numDOF2) = dsp2(i);  ugdot(i + numDOF2) = vel2(i);
	}

	// transform response from the global to the local system
	ul.addMatrixVector(0.0, Tgl, ug, 1.0);
	uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);

	// transform response from the local to the basic system
	ub.addMatrixVector(0.0, Tlb, ul, 1.0);
	ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);
	//ub = (Tlb*Tgl)*ug;
	//ubdot = (Tlb*Tgl)*ugdot;

	// set trial response for material models
	for (int i = 0; i < 3; i++)
		errCode += MaterialPtr[i]->setTrialStrain(ub(i), ubdot(i));

	return errCode;
}

const Matrix &
BeamColumnJoint2dThermal::getTangentStiff(void)
{
	// zero the matrix
	theMatrix->Zero();

	return *theMatrix;
}

const Matrix &
BeamColumnJoint2dThermal::getInitialStiff(void)
{
	return getTangentStiff();
}

const Vector &
BeamColumnJoint2dThermal::getResistingForce(void)
{
	// zero the residual
	theVector->Zero();

	// get resisting forces
	for (int i = 0; i < 3; i++)
		qb(i) = MaterialPtr[i]->getStress();

	// determine resisting forces in local system
	//Vector ql(6);
	//ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);

	// add P-Delta effects to local forces

	// determine resisting forces in global system    
	//theVector->addMatrixTransposeVector(0.0, Tgl, ql, 1.0);

	return *theVector;
}


    
const Matrix &
BeamColumnJoint2dThermal::getDamp(void)
{
	//not applicable (stiffness being returned)
	return *theMatrix;
}

const Matrix &
BeamColumnJoint2dThermal::getMass(void)
{ 
	//not applicable  (stiffness being returned)
	return *theMatrix;
}

void 
BeamColumnJoint2dThermal::zeroLoad(void)
{
	//not applicable  
	return;
}

int 
BeamColumnJoint2dThermal::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	//not applicable
	return 0;
}

int 
BeamColumnJoint2dThermal::addInertiaLoadToUnbalance(const Vector &accel)
{
	//not applicable
	return 0;
}


const Vector &
BeamColumnJoint2dThermal::getResistingForceIncInertia()
{	
  //not applicable (residual being returned)
	return *theVector;
}

int
BeamColumnJoint2dThermal::sendSelf(int commitTag, Channel &theChannel)
{
	// yet to do.
	return -1;
}

int
BeamColumnJoint2dThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// yet to do.
	return -1;
}

int
BeamColumnJoint2dThermal::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	const Vector &node1Crd = nodePtr[0]->getCrds();
	const Vector &node2Crd = nodePtr[1]->getCrds();	


	const Vector &node1Disp = nodePtr[0]->getDisp();
	const Vector &node2Disp = nodePtr[1]->getDisp();    
 

	static Vector v1(3);
	static Vector v2(3);

	
	// calculate the current coordinates of four external nodes
	for (int i=0; i<2; i++) 
    {
		v1(i) = node1Crd(i)+node1Disp(i)*fact;
		v2(i) = node2Crd(i)+node2Disp(i)*fact;

	}
	

	return 0;

}

void
BeamColumnJoint2dThermal::Print(OPS_Stream &s, int flag)
{
	s << "Element: " << this->getTag() << " Type: Beam Column Joint " << endln;
	for (int i = 0; i<2; i++)
	{
		s << "Node :" << connectedExternalNodes(i);
		s << "DOF :" << nodePtr[i]->getNumberDOF();
	}
	return;
}

Response*
BeamColumnJoint2dThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  // we will compare argv[0] to determine the type of response required
  
  if (strcmp(argv[0],"node1BarSlipL") == 0 || strcmp(argv[0],"node1BarslipL") ==0 || strcmp(argv[0],"Node1BarSlipL") == 0)
    return MaterialPtr[0]->setResponse(&argv[1], argc-1, output);
  
  else if (strcmp(argv[0],"node1BarSlipR") == 0 || strcmp(argv[0],"node1BarslipR") ==0 || strcmp(argv[0],"Node1BarSlipR") == 0)
    return MaterialPtr[1]->setResponse(&argv[1], argc-1, output);
  
  else if (strcmp(argv[0],"node1InterfaceShear") == 0 || strcmp(argv[0],"node1Interfaceshear") ==0 || strcmp(argv[0],"Node1InterfaceShear") ==0 )
    return MaterialPtr[2]->setResponse(&argv[1], argc-1, output);
  
  else if (strcmp(argv[0],"node2BarSlipB") == 0 || strcmp(argv[0],"node2BarslipB") ==0 || strcmp(argv[0],"Node2BarSlipB") == 0)
    return MaterialPtr[3]->setResponse(&argv[1], argc-1, output);
  
  else if (strcmp(argv[0],"node2BarSlipT") == 0 || strcmp(argv[0],"node2BarslipT") ==0 || strcmp(argv[0],"Node2BarSlipT") == 0)
    return MaterialPtr[4]->setResponse(&argv[1], argc-1, output);
  
  else if (strcmp(argv[0],"node2InterfaceShear") == 0 || strcmp(argv[0],"node2Interfaceshear") ==0 || strcmp(argv[0],"Node2InterfaceShear") ==0 )
    return MaterialPtr[5]->setResponse(&argv[1], argc-1, output);

  else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"Deformation") == 0)
    return new ElementResponse(this,3,Vector(4));
  
  else
    return 0;
}

int 
BeamColumnJoint2dThermal::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1:  // global forces
		return eleInfo.setVector(this->getResistingForce());

	default:
		return -1;
	}
}

int
BeamColumnJoint2dThermal::setParameter (char **argv, int argc, Information &info)
{
  return -1;
}
    
int
BeamColumnJoint2dThermal::updateParameter (int parameterID, Information &info)
{
  return -1;
}
