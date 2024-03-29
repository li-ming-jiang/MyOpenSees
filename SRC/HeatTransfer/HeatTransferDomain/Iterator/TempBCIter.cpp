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
// Adapted by Yaqiang Jiang (y.jiang@ed.ac.uk)
//

#include <TempBCIter.h>
#include <TemperatureBC.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


TempBCIter::TempBCIter(TaggedObjectStorage* theStorage)
:myIter(theStorage->getComponents())
{

}

TempBCIter::~TempBCIter()
{

}    


void
TempBCIter::reset(void)
{
    myIter.reset();
}


TemperatureBC*
TempBCIter::operator()(void)
{
    // check if we still have SP_Constraints in the model
    // if not return 0, indicating we are done
    TaggedObject* theComponent = myIter();
    if (theComponent == 0)
		return 0;
    else {
		TemperatureBC* result = (TemperatureBC*)theComponent;
		return result;
		}
}