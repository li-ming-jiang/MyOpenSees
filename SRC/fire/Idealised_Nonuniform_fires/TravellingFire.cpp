/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Fire & Heat Transfer modules developed by:                         **
**   Yaqiang Jiang (y.jiang@ed.ac.uk)                                 **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */

//
// Written by Yaqiang Jiang (y.jiang@ed.ac.uk)
//

#include <TravellingFire.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <PrescribedSurfFlux.h>
#include <ID.h>
#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
#include <PathTimeSeriesThermal.h>


TravellingFire::TravellingFire(int tag, PathTimeSeriesThermal* fireLocPath, double D,
	double Q, double H, int lineTag)
	:FireModel(tag, 7), FireLocPath(fireLocPath), fireLocs(3), d(D), q(Q), h(H), centerLine(lineTag)
{
    // check the direction of central line of a Hasemi fire
    // 1 indicates it is parrallel to x1 axis, 2 indicates
    // parallelt to x2 axis, 3 indicates parallel to x3 axis.
    if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "TravellingFire::TravellingFire - invalid line tag provided for Hasemi fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
		}
	if ((d > 10.0) || (q > 5e7)) {
		opserr << "TravellingFire::TravellingFire - error in specifying the fire, diameter "
			<< " shoudn't be greater than 10m, fire size shoudn't be greater than 50MW.\n";
		}
}


TravellingFire::~TravellingFire()
{

}



int
TravellingFire::setFirePars(double time) {
	
		Vector firePars = FireLocPath->getFactors(time);

		if (firePars.Size() == 3) {
			fireLocs(0) = firePars(0);
			fireLocs(1) = firePars(1);
			fireLocs(2) = firePars(2);
		}
		else if (firePars.Size() == 4) {
			fireLocs(0) = firePars(0);
			fireLocs(1) = firePars(1);
			fireLocs(2) = firePars(2);
			q = firePars(3);
		}
		else if (firePars.Size() == 5) {
			fireLocs(0) = firePars(0);
			fireLocs(1) = firePars(1);
			fireLocs(2) = firePars(2);
			q = firePars(3);
			d = firePars(4);
		}
		else {
			opserr << "WARNING! TravellingFire::getFlux failed to get the location of fire origin" << endln;
		}

		return 0;
}

double
TravellingFire::getFirePars(int ParTag) {
	if (ParTag == 1)
		return fireLocs(0);
	else if (ParTag == 2)
		return fireLocs(1);
	else if (ParTag == 3)
		return fireLocs(2);
	else if (ParTag == 4)
		return q;
	else if (ParTag == 5)
		return d;
	else {
		opserr << "WARNING! invalid tag for TravellingFire::getFirePars " << ParTag << endln;
		return -1;
	}
		

}

double 
TravellingFire::getFlux(HeatTransferNode* node, double time)
{
	if (FireLocPath != 0) {
		this->setFirePars(time);
	}
	double q_dot=0;
	
	// first calculate flame length
    double Lf = 0.0148 * pow(q,0.4) - 1.02 * d;
	if (Lf < h) {
		//opserr << "TravellingFire::getFlux() - flame is not impinging ceiling, method has not implemented.\n";
		q_dot = 0.0;
		}
	double constant = 1.11 * 1e6 * pow(d,2.5);
	double Qd_ast = q / constant;
	double z_acute;

	if (Qd_ast < 1.0) {
		double a = 0.66666666666666666666666666666667;
		double term = pow(Qd_ast,0.4) - pow(Qd_ast,a);
		z_acute = 2.4 * d * term;
		} else {
			double term = 1.0 - pow(Qd_ast, 0.4);
			z_acute = 2.4 * d * term;
		}
	double term = 1.11 * 1e6 * pow(h,2.5);
	double Qh_ast = q / term;

	// now calculate H plus Lh
	double Lt = 2.9 * h * pow(Qh_ast,0.33);

	// now calculate r	
	const Vector& coords = node->getCrds();
	int size = coords.Size();

	double deltaX1, deltaX2, sum, r;
   
	
	
	if (centerLine == 1) {
		deltaX1 = fireLocs(1) - coords(1);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = fireLocs(2) - coords(2);
			} else {
				// then heat transfer coordinate is 2D
				deltaX2 = fireLocs(2);
			}
	} else if (centerLine == 2) {
		deltaX1 = fireLocs(0) - coords(0);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = fireLocs(2) - coords(2);
		} else {
			// then heat transfer coordinate is 2D
			deltaX2 = fireLocs(2);
		}
	} else if (centerLine == 3) {
			deltaX1 = fireLocs(0) - coords(0);
			deltaX2 = fireLocs(1) - coords(1);
			}

    sum = deltaX1 * deltaX1 + deltaX2 * deltaX2;
    r = sqrt(sum);

	// now calculate y
	double y = (r + h + z_acute) / (Lt + z_acute);

	// now determine the flux
	
	if (y <= 0.3) {
		q_dot = 100000;
	} else if ((y > 0.3) && (y < 1.0)) {
		q_dot = 136300 - 121000 * y;
	} else if (y >= 1.0) {
		q_dot = 15000 * pow(y,-3.7);
	}
	int tag = node->getTag();
#ifdef _DEBUG
	//opserr<<" Travelling fire: "<<fireLocs(1)<<","<<x3<<", q"<<q<<" r:  "<<r<<"____ "<< " q:  "<<q_dot<<"____ ";
#endif
	return q_dot;
}


void
TravellingFire::applyFluxBC(HeatFluxBC* theFlux, double time)
{
    int flux_type = theFlux->getTypeTag();
    if (flux_type == 3) 
		{
		PrescribedSurfFlux* pflux = (PrescribedSurfFlux*) theFlux;

		//int flux_type = pflux->getTypeTag();
		int eleTag = pflux->getElementTag();
		int fTag = pflux->getFaceTag();
		HeatTransferDomain* theDomain = pflux->getDomain();
		if (theDomain == 0) {
			opserr << "TravellingFire::applyFluxBC() - HeatFluxBC has not been associated with a domain";
			exit(-1);
			}

		HeatTransferElement* theEle = theDomain->getElement(eleTag);
		if (theEle == 0) {
			opserr << "TravellingFire::applyFluxBC() - no element with tag " << eleTag << " exists in the domain";
			exit(-1);
			}

		const ID& faceNodes = theEle->getNodesOnFace(fTag);
		int size = faceNodes.Size();
		Vector nodalFlux(size);

		for (int i = 0; i < size; i++) {
			int nodTag = faceNodes(i);
			HeatTransferNode* theNode = theDomain->getNode(nodTag);
			if (theNode == 0) {
				opserr << "TravellingFire::applyFluxBC() - no node with tag " << nodTag << " exists in the domain";
				exit(-1);
				}
			nodalFlux(i) = this->getFlux(theNode,time); 
			//opserr << "Flux at node " << nodTag << " is " << nodalFlux(i) << endln;
			}

		pflux->setData(nodalFlux);
		pflux->applyFluxBC();
		} else {
			opserr << "TravellingFire::applyFluxBC() - incorrect flux type "
				<< flux_type << " provided\n";
			exit(-1);
		}
}

