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
// Note: This class was adapted from Material   
#include <HeatTransferMaterial.h>

HeatTransferMaterial::HeatTransferMaterial(int tag)
:TaggedObject(tag), k(0)
{


}


HeatTransferMaterial::~HeatTransferMaterial()
{
  // does nothing


}



Response*
HeatTransferMaterial::setResponse(const char** argv, int argc,
    OPS_Stream& theOutput)
{
    Response* theResponse = 0;

    if ((strcmp(argv[0], "TempElong") == 0) {
            theOutput.tag("ResponseType", "temp11");
            theOutput.tag("ResponseType", "Elong11");
            theResponse = new MaterialResponse(this, 7, Vector(2));
        }

        theOutput.endTag();
    }

    return theResponse;

}

int
HeatTransferMaterial::getResponse(int responseID, Information& matInfo)
{
    static Vector stressStrain(2);
    static Vector stressStrainTangent(3);

    static Vector tempData(2);  //L.jiang [SIF]
    static Information infoData(tempData);  //L.jiang [SIF]

    // each subclass must implement its own stuff   

    // added for sensitivity recorder. Quan 2009
    if ((responseID > 10000) && (responseID < 20000)) {
        matInfo.setDouble(this->getStressSensitivity(responseID - 10000, false));
        return 0;
    }
    else if (responseID > 20000) {
        matInfo.setDouble(this->getStrainSensitivity(responseID - 20000));
        return 0;
    }

    double kInit;
    double stress;
    double strain;

    switch (responseID) {
    case 1:
        matInfo.setDouble(this->getStress());
        return 0;
    default:
        return -1;
    }
}