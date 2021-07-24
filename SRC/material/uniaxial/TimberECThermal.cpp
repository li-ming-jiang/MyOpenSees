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

// Added by Liming Jiang (UoE)
// Created: 06/13
//
// Description: This file contains the class definition for 
// TimberECThermal. TimberECThermal is modified from Concrete02Thermal
// TimberECThermal is dedicated to provide a concrete material which 
// strictly satisfy Eurocode regarding the temperature dependent properties.

// Concrete02 is written by FMK in the year of 2006 and based on Concr2.f



#include <stdlib.h>
#include <TimberECThermal.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

void*
OPS_TimberECThermal(void)
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial* theMaterial = 0;

    int    iData[1];
    double dData[6];
    int numData = 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial TimberECThermal tag" << endln;
        return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();

    if (numData != 6) {
        opserr << "Invalid #args, want: uniaxialMaterial TimberECThermal " << iData[0] << "fc? epsc0? fcu? Epos? Eneg? ft?\n";
        return 0;
    }

    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid #args, want: uniaxialMaterial TimberECThermal " << iData[0] << "fc? epsc0? fcu? Epos? Eneg? ft?\n";
        return 0;
    }


    // Parsing was successful, allocate the material
    theMaterial = new TimberECThermal(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type TimberECThermal Material\n";
        return 0;
    }

    return theMaterial;
}


TimberECThermal::TimberECThermal(int tag, double _fc, double _epsc0, double _fcu,
    double _Epos, double _Eneg, double _ft) :
    UniaxialMaterial(tag, MAT_TAG_ConcreteECThermal),
    //fc(_fc), epsc0(_epsc0), fcu(_fcu), Epos(_Epos), Eneg(_Eneg), ft(_ft))
    fcT(_fc), epsc0T(_epsc0), fcuT(_fcu), EposT(_Epos), EnegT(_Eneg), ftT(_ft),E(_Epos)
{

    //JZ 07/10 /////////////////////////////////////////////////////////////start
    fc = fcT;
    epsc0 = epsc0T;
    fcu = fcuT;
    Epos = EposT;
    Eneg = EnegT;
    ft = ftT;
    //JZ 07/10 /////////////////////////////////////////////////////////////end 

//    eposP = Epos;
//    enegP = Eneg;
    //eP = 1.5*fc/epsc0; //for the euro code, the 2.0 should be changed into 1.5
    epsP = 0.0;
    sigP = 0.0;
    eps = 0.0;
    sig = 0.0;
    EposP = Epos;
    EnegP = Eneg;
    //e = 1.5*fc/epsc0;//for the euro code, the 2.0 should be changed into 1.5
     
    ThermalElongation = 0; //initialize 
    TempP = 20.0; //Pk add previous temp


}

TimberECThermal::TimberECThermal(void) :
    UniaxialMaterial(0, MAT_TAG_ConcreteECThermal)
{

}

TimberECThermal::~TimberECThermal(void)
{
    // Does nothing
}




UniaxialMaterial*
TimberECThermal::getCopy(void)
{
    TimberECThermal* theCopy = new TimberECThermal(this->getTag(), fc, epsc0, fcu, Epos, Eneg, ft);

    return theCopy;
}

double
TimberECThermal::getInitialTangent(void)
{
    return E;
}

int
TimberECThermal::setTrialStrain(double trialStrain, double FiberTemperature, double strainRate)
{
	//---------------------------------------------The Main stress determination--------------------
    fcu = 0.85 * fc;
    epsc0 = 1.25 * (fc / Epos);

    double k1 = (fcu) / (3 * Epos * (pow(epsc0, 4)) * (1 - fcu / fc));
    double k2 = 1 / Epos;
    double k3 = 1 / fc - 4 / (3 * Epos * (epsc0));
    double k4 = 1 / (3 * Epos * (pow(epsc0, 4)) * (1 - fcu / fc));
    eps = trialStrain;

    if (eps < 0) {
        sig = (eps + k1 * pow(eps, 4)) / (k2 + eps * k3 + k4 * pow(eps, 4));
        EposT = ((1 + 4 * k1 * pow(eps, 3)) * (k2 + k3 * eps + k4 * pow(eps, 4)) - (eps + k1 * pow(eps, 4)) * (k3 + 4 * k4 * pow(eps, 3))) / pow((k2 + k3 * eps + k4 * pow(eps, 4)), 2);
        E = EposT;
    }
        
    else {
        sig = Eneg * eps;
        E = Eneg;
        if (eps > (ft / Eneg)) {
            sig = ft;
            E = 1e-3;
        }
      
    }
   
   // opserr<<"trialStrain: "<<eps << "  Stress: "<<sig<< "Modulus: "<<E<<endln;
    return 0;
}



double
TimberECThermal::getStrain(void)
{
    return eps;
}

double
TimberECThermal::getStress(void)
{
    return sig;
}

double
TimberECThermal::getTangent(void)
{
    return E;
}

double
TimberECThermal::getThermalElongation(void) //***JZ
{
    return ThermalElongation;
}

double
TimberECThermal::getElongTangent(double TempT, double& ET, double& Elong, double TempTmax) //PK add to include max temp
{
    //material properties with temperature
    Temp = TempT;  //make up the 20 degree which is minus in the class of thermalfield
//    Tempmax = TempTmax; //PK add max temp for cooling
    // The datas are from EN 1992 part 1-2-1 
    // Tensile strength at elevated temperature

     //if (Temp >= 1080) {
    //	  opserr << "temperature " << " " << Temp <<endln;
    //}
    if (Temp <= 20) {
        ft = ftT;
        Eneg = EnegT;
    }
    else if (Temp <= 100) {
        ft = (1.0 - 0.35 * (Temp - 20) / 80) * ftT;
        Eneg = (1.0 - 0.50 * (Temp - 20) / 80) * EnegT;
        //Ets = (1.0 - 1.0*(Temp -80)/500)*EtsT;
    }
    else if (Temp <= 300) {
        ft = (0.65 - 0.65 * (Temp - 100) / 200) * ftT;
        Eneg = (0.50 - 0.50 * (Temp - 100) / 200) * EnegT;
    }
    else {
        ft = 1.0e-10;
        Eneg = 1.0e-10;
        //ft = 0;
        //Ets = 0;
    }

    // compression strength, at elevated temperature
    //   strain at compression strength, at elevated temperature
    //   ultimate (crushing) strain, at elevated temperature
    if (Temp <= 20) {
        fc = fcT;
        fcu = fcuT;
        Epos = EposT;
        //Ets = EtsT;  jz what is there the statement?
    }
    else if (Temp <= 100) {
        fc = fcT * (1 - ((Temp - 20) / 80) * 0.75);
        fcu = fcuT * (1 - ((Temp - 20) / 80) * 0.75);
        Epos = EposT * (1.0 - 0.65 * (Temp - 20) / 80);
    }
    else if (Temp <= 300) {
        fc = fcT * (0.25 - ((Temp - 100) / 200) * 0.25);
        fcu = fcuT * (0.25 - ((Temp - 100) / 200) * 0.25);
        Epos = EposT * (0.35 - 0.35 * (Temp - 100) / 200);
    }
    else {
        fc= 1.0e-10;
        fcu = 1.0e-10;
        Epos = 1.0e-10;
    }

    return 0;
}






int
TimberECThermal::commitState(void)
{
//    ecminP = ecmin;
//    deptP = dept;

    
    EposP = Epos;
    EnegP = Eneg;
    sigP = sig;
    epsP = eps;

    TempP = Temp; //PK add set the previous temperature

    return 0;
}

int
TimberECThermal::revertToLastCommit(void)
{
//    ecmin = ecminP;;
//    dept = deptP;


    Epos = EposP;
    Eneg = EnegP;
    sig = sigP;
    eps = epsP;

    Temp = TempP; //PK add set the previous temperature

    // NA ELENXW MIPWS EDW XANETAI TO TEMP LOGW MIN CONVERGENCE

    return 0;
}

int
TimberECThermal::revertToStart(void)
{
//    ecminP = 0.0;
//    deptP = 0.0;

    EposP = Epos;
    EnegP = Eneg;

    sig = 0.0;
    eps = 0.0;
    sigP = 0.0;
    epsP = 0.0;


    return 0;
}

int
TimberECThermal::sendSelf(int commitTag, Channel& theChannel)
{
    static Vector data(9);
    data(0) = fc;
    data(1) = epsc0;
    data(2) = fcu;
//    data(3) = epscu;
//    data(4) = rat;
    data(3) = EposP;
//    data(6) = Ets;
//    data(7) = ecminP;
//    data(8) = deptP;
    data(4) = EnegP;
    data(5) = ft;
    data(6) = sigP;
    data(7) = epsP;
    data(8) = this->getTag();

    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "TimberECThermal::sendSelf() - failed to sendSelf\n";
        return -1;
    }
    return 0;
}

int
TimberECThermal::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{

    static Vector data(9);

    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "TimberECThermal::recvSelf() - failed to recvSelf\n";
        return -1;
    }

    fc = data(0);
    epsc0 = data(1);
      fcu = data(2);
//  epscu = data(3);
//    rat = data(4);
    EposP = data(3);
    EnegP = data(4);
       ft = data(5);
     sigP = data(6);
     epsP = data(7);
     this->setTag(data(8));
//    Ets = data(6);
//    ecminP = data(7);
//    deptP = data(8);

      eps = epsP;
      sig = sigP;
    Eneg = EnegP;
    Epos = EposP;




    return 0;
}

void
TimberECThermal::Print(OPS_Stream& s, int flag)
{
    s << "TimberECThermal:(strain, stress, tangent) " << eps << " " << sig << " " << (Epos, Eneg) << endln;
}


int
TimberECThermal::getVariable(const char* varName, Information& theInfo)
{
    if (strcmp(varName, "ec") == 0) {
        theInfo.theDouble = epsc0;
        return 0;
    }
    else if (strcmp(varName, "ElongTangent") == 0) {
        Vector* theVector = theInfo.theVector;
        if (theVector != 0) {
            double tempT, ET, Elong, TempTmax;
            tempT = (*theVector)(0);
            ET = (*theVector)(1);
            Elong = (*theVector)(2);
            TempTmax = (*theVector)(3);
            this->getElongTangent(tempT, ET, Elong, TempTmax);
            (*theVector)(0) = tempT;
            (*theVector)(1) = ET;
            (*theVector)(2) = Elong;
            (*theVector)(3) = TempTmax;
        }
        return 0;
    }
    else if (strcmp(varName, "ThermalElongation") == 0) {
        theInfo.theDouble = ThermalElongation;
        return 0;
    }
    else if (strcmp(varName, "TempAndElong") == 0) {
        Vector* theVector = theInfo.theVector;
        if (theVector != 0) {
            (*theVector)(0) = Temp;
            (*theVector)(1) = ThermalElongation;
        }
        else {
            opserr << "null Vector in EC" << endln;
        }

        return 0;
    }
    return -1;
}



//this function is no use, just for the definiation of pure virtual function.
int
TimberECThermal::setTrialStrain(double strain, double strainRate)
{
    opserr << "TimberECThermal::setTrialStrain(double strain, double strainRate) - should never be called\n";
    return -1;
}
