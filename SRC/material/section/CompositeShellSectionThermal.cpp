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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-05-21 23:03:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/CompositeShellSectionThermal.cpp,v $

// To define shell section for composite floor slabs
//
// Composite Shell Section
//
// Added for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 


#include <CompositeShellSectionThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <elementAPI.h>

//static vector and matrices
Vector  CompositeShellSectionThermal::stressResultant(8) ;
Matrix  CompositeShellSectionThermal::tangent(8,8) ;
//ID      CompositeShellSectionThermal::array(8) ;
//const double  CompositeShellSectionThermal::root56 = sqrt(5.0/6.0) ; //shear correction
//const double  CompositeShellSectionThermal::root56 = 1 ; //shear correction
//null constructor
void* OPS_CompositeShellSectionThermal()
{
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "WARNING insufficient arguments" << endln;
        opserr << "Want: section CompositeShellThermal tag? SectionTag1? r1? .SectionTag2? r2? " << endln;
        return 0;
    }

    int count = 2;;
    int tag, sectag1, sectag2;
    double ratio1, ratio2;

    double ribAngle = 0.0;
    SectionForceDeformation* theSection = 0;
    SectionForceDeformation* theSec1 = 0;
    SectionForceDeformation* theSec2 = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid section CompositeShellThermal tag" << endln;
        return 0;
    }
    
    if (OPS_GetIntInput(&numdata, &sectag1) < 0) {
        opserr << "WARNING invalid section tag 1" << endln;
        opserr << "CompositeShellThermal section: " << tag << endln;
        return 0;
    }
    
    if (OPS_GetDoubleInput(&numdata, &ratio1) < 0) {
        opserr << "WARNING invalid ratio1\n";
        opserr << "CompositeShellThermal section:  " << tag << endln;
        return 0;
    }
    
    if (OPS_GetIntInput(&numdata, &sectag2) < 0) {
        opserr << "WARNING invalid section tag 2" << endln;
        opserr << "CompositeShellThermal section: " << tag << endln;
        return 0;
    }

    if (OPS_GetDoubleInput(&numdata, &ratio2) < 0) {
        opserr << "WARNING invalid ratio2\n";
        opserr << "CompositeShellThermal section:  " << tag << endln;
        return 0;
    }

    if (OPS_GetNumRemainingInputArgs() < 1) {
        ribAngle = 0;
    }
    else {
        if (OPS_GetDoubleInput(&numdata, &ribAngle) < 0) {
            opserr << "WARNING invalid ribAngle\n";
            opserr << "CompositeShellThermal section:  " << tag << endln;
            return 0;
        }
    }

    theSec1 = OPS_GetSectionForceDeformation(sectag1);
    theSec2 = OPS_GetSectionForceDeformation(sectag2);

    theSection = new CompositeShellSectionThermal(tag, theSec1, theSec2, ratio1, ratio2, ribAngle);

    return theSection;

}










CompositeShellSectionThermal::CompositeShellSectionThermal( ) : 
SectionForceDeformation( 0, SEC_TAG_LayeredShellFiberSectionThermal ), 
strainResultant(8), theSection1(0),theSection2(0),Ratio1(0.0), Ratio2(0.0), ThermalElongation(0), ribAng(0.0)
{
    stiffratio1 = 0;
    stiffratio2 = 0;
    stiffratio3 = 0;
}

//full constructor
CompositeShellSectionThermal::CompositeShellSectionThermal(    
                                   int tag, 
                                   SectionForceDeformation* theSec1, 
                                   SectionForceDeformation* theSec2,
                                   double ratio1, double ratio2, double ribAngle) :
SectionForceDeformation( tag, SEC_TAG_LayeredShellFiberSectionThermal ),
strainResultant(8), theSection1(theSec1), theSection2(theSec2), Ratio1(ratio1),Ratio2(ratio2), ThermalElongation(0), ribAng(ribAngle)
{
    sT = new Vector(2);
    sT->Zero();

    stiffratio1 = 0;
    stiffratio2 = 0;
    stiffratio3 = 0;
}

//destructor
CompositeShellSectionThermal::~CompositeShellSectionThermal( ) 
{ 

} 

//make a clone of this material
SectionForceDeformation  *CompositeShellSectionThermal::getCopy( ) 
{
  CompositeShellSectionThermal *clone = 0;   //new instance of this class
  
  if (theSection1 != 0&&theSection2!=0)
  {

    clone = new CompositeShellSectionThermal( this->getTag(),
					  theSection1,
					  theSection2,
					  Ratio1, Ratio2,ribAng) ; //make the copy
  }
  return clone ;
}



//send back order of strainResultant in vector form
int CompositeShellSectionThermal::getOrder( ) const
{
  return 8 ;
}


//send back order of strainResultant in vector form
const ID& CompositeShellSectionThermal::getType( ) 
{
    return theSection1->getType();
}



//swap history variables
int CompositeShellSectionThermal::commitState( ) 
{
  int success = 0 ;
  /*
  for (int i = 0; i < nLayers; i++ )
    success += theFibers[i]->commitState( ) ;
  */
  

  return success ;
}

//revert to last saved state
int CompositeShellSectionThermal::revertToLastCommit( )
{
  int success = 0 ;

 // for (int i = 0; i < nLayers; i++ )
 //   success += theFibers[i]->revertToLastCommit( ) ;

  return success ;
}

//revert to start
int CompositeShellSectionThermal::revertToStart( )
{
  int success = 0 ;

 // for (int i = 0; i < nLayers; i++ )
  //  success += theFibers[i]->revertToStart( ) ;

  return success ;
}

//mass per unit area
double
CompositeShellSectionThermal::getRho( )
{

  double weight ;

  double rhoH = 0.0 ;

 // for ( int i = 0; i < nLayers; i++ ) {
   
  //  weight = ( 0.5*h ) * wg[i] ;

    //rhoH += ( theFibers[i]->getRho() ) * weight ;

  //}

  return rhoH ;

}





//receive the strainResultant 
int CompositeShellSectionThermal ::
setTrialSectionDeformation( const Vector &strainResultant_from_element)
{
 
  strainResultant = strainResultant_from_element ;

  Vector strain1(8);
  Vector strain2(8);
  strain1 = strainResultant_from_element;
  strain2 = strainResultant_from_element;
  //stiffratio = kb/ka;
  // ka A = kb B;
  //w1A + w2B =C;
  // A = stiffratio/(w1*stiffratio+w2)*C; B = 1/(w1*ratio+w2)*C

  if (abs(ribAng) < 1e-5) {
      //rib along the local x direction// for Qiu
      strain1(1) = stiffratio1 / (Ratio1 * stiffratio1 + Ratio2) * strainResultant_from_element(1);
      strain2(1) = 1 / (Ratio1 * stiffratio1 + Ratio2) * strainResultant_from_element(1);

      strain1(4) = stiffratio2 / (Ratio1 * stiffratio2 + Ratio2) * strainResultant_from_element(4);
      strain2(4) = 1 / (Ratio1 * stiffratio2 + Ratio2) * strainResultant_from_element(4);
  }
  else if (abs(ribAng) - 90 < 1e-5) {
      //rib along the local y direction// for Qiu


  }


  //Same strain along ribs; Need to add auto-identification of rib direction
  //Shear strains are set as the same for current solution.
  int success = 0 ;

  int i ;

  double z ;

 #ifdef _SDEBUG
	  if (strainResultant(8) ==111)
          opserr<< "Sec strain  "<<strainResultant<<endln;
#endif

      strain1 = 

      success = theSection1->setTrialSectionDeformation(strain1);
      success= success+ theSection2->setTrialSectionDeformation(strain2);

/*
for ( i = 0; i < nLayers; i++ ) {

      z = ( 0.5*h ) * sg[i]-Offset;

      strain(0) =  strainResultant(0)  - z*strainResultant(3)-ThermalElongation[i] ;

      strain(1) =  strainResultant(1)  - z*strainResultant(4)-ThermalElongation[i] ;

      strain(2) =  strainResultant(2)  - z*strainResultant(5) ;

      strain(3) =  root56*strainResultant(6) ;

      strain(4) =  root56*strainResultant(7) ;

    #ifdef _SDEBUG
      if (strainResultant(8) ==111&&i==0){
      opserr<<"LayeredShellSection z : "<<z<< "strain  "<<strain<<endln;
      strain(5)=1110;
      }
      else
        strain(5)=1115;
    #endif

//#ifdef _DEBUG
//opserr<<" Layer"<< i <<endln;
//#endif
      success += theFibers[i]->setTrialStrain( strain ) ;
      #ifdef _SDEBUG
      if (strainResultant(8) ==111&&i==0){
      opserr<<"LayeredShellSection z : "<<z<<"stress  "<<theFibers[i]->getStress()<<endln;
      opserr<<" tangent "<<theFibers[i]->getTangent()<<endln;
          }
    #endif

  } //end for i
*/
  

  return success ;
}


//send back the strainResultant
const Vector& CompositeShellSectionThermal::getSectionDeformation( )
{
  return this->strainResultant ;
}


const Vector&
CompositeShellSectionThermal::getTemperatureStress(const Vector& dataMixed)
{
   //this->countnGauss = 0;             
  /*
  double* ThermalTangent = new double[nLayers];
  //double ThermalElongation[5];
  for (int i = 0; i < nLayers; i++) {
          ThermalTangent[i]=0;
          ThermalElongation[i]=0;
  }
  double FiberTemperature = 0 ;
  double tangent, elongation;
  double averageThermalForce =0.0;
  double averageThermalMoment =0.0;
  //double *thickness = new double[nLayers];

  for (int i = 0; i < nLayers; i++) {
	  
	double thickness = 0.5*h*wg[i];

    double yi = ( 0.5*h ) * sg[i] - Offset;

	double tangent, elongation;

	FiberTemperature = this->determineFiberTemperature( dataMixed, yi);

	theFibers[i]->getThermalTangentAndElongation(FiberTemperature, tangent, elongation);

	ThermalTangent[i]=tangent;
	ThermalElongation[i]=elongation;
	averageThermalForce += elongation*thickness*tangent;
	averageThermalMoment+= yi*thickness*elongation*tangent;
	}

      (*sT)(0) = averageThermalForce - AverageThermalForceP;

      (*sT)(1) = averageThermalMoment - AverageThermalMomentP;
	  AverageThermalForceP = averageThermalForce;
	  AverageThermalMomentP = averageThermalMoment;
     
  */
    Vector vec1(18);
    Vector vec2(18);
    if (dataMixed.Size() != 36)
        opserr << "Warning::Composite section recived incorrect thermal action" << endln;
    
    for (int i = 0; i < 18; i++) {
        vec1(i) = dataMixed(i);
        vec2(i) = dataMixed(i+18);
    }
   // opserr << "vec1::" << vec1 << endln;
   // opserr << "vec2::" << vec2 << endln;
    theSection1->getTemperatureStress(vec1);
    theSection2->getTemperatureStress(vec2);
  
    return *sT;
}


//send back the stressResultant 
const Vector&  CompositeShellSectionThermal::getStressResultant( )
{

  static Vector stress(5) ;

  int i ;

  double z, weight ;

  stressResultant.Zero( ) ;

  Vector stress1 = theSection1->getStressResultant();
  Vector stress2 = theSection2->getStressResultant();

  stressResultant = stress1 * Ratio1 + stress2 * Ratio2;
  //currently using ratios of sections

  /*
  
  for ( i = 0; i < nLayers; i++ ) {

      z = ( 0.5*h ) * sg[i]- Offset;  //added for offset

      weight = ( 0.5*h ) * wg[i] ;

      stress = theFibers[i]->getStress( ) ;
  
      //membrane
      stressResultant(0)  +=  stress(0)*weight ;

      stressResultant(1)  +=  stress(1)*weight ;

      stressResultant(2)  +=  stress(2)*weight ;

      //bending moments
      stressResultant(3)  +=  ( z*stress(0) ) * weight ;

      stressResultant(4)  +=  ( z*stress(1) ) * weight ;

      stressResultant(5)  +=  ( z*stress(2) ) * weight ;

      //shear
      stressResultant(6)  += stress(3)*weight ;

      stressResultant(7)  += stress(4)*weight ;
  
  } //end for i

     //modify shear 
   stressResultant(6) *= root56 ;  
   stressResultant(7) *= root56 ;
  */
  


   return this->stressResultant ;
}


//send back the tangent 
const Matrix& CompositeShellSectionThermal::getSectionTangent( )
{
  static Matrix dd(5,5) ;

//  static Matrix Aeps(5,8) ;

//  static Matrix Asig(8,5) ;

  int i ;

  double z, weight ;

  tangent.Zero( ) ;

  Matrix Tangent1 = theSection1->getSectionTangent();
  Matrix Tangent2 = theSection2->getSectionTangent();

  tangent = Tangent1 * Ratio1 + Tangent2 * Ratio2;

  //To determine stiff ratios
  stiffratio1 = Tangent2(1, 1) / Tangent1(1, 1);
  stiffratio2 = Tangent2(3, 3) / Tangent1(3, 3);

  /*
  for ( i = 0; i < nLayers; i++ ) {

      z = ( 0.5*h ) * sg[i] - Offset;

      weight = (0.5*h) * wg[i] ; ----*/

   /*      //compute Aeps

      Aeps.Zero( ) ;

      Aeps(0,0) = 1.0 ;
      Aeps(0,3) = -z ;

      Aeps(1,1) = 1.0 ;
      Aeps(1,4) = -z ;

      Aeps(2,2) = 1.0 ;
      Aeps(2,5) = -z ;

      Aeps(3,6) = root56 ;
      Aeps(4,7) = root56 ;

      //compute Asig

      Asig.Zero( ) ;

      Asig(0,0) = 1.0 ;
      Asig(3,0) = z ;

      Asig(1,1) = 1.0 ;
      Asig(4,1) = z ;

      Asig(2,2) = 1.0 ;
      Asig(5,2) = z ;

      Asig(6,3) = root56 ;
      Asig(7,4) = root56 ;
*/

//compute the tangent
/*
  dd = theFibers[i]->getTangent();

  dd *= weight;

  //tangent +=  ( Asig * dd * Aeps ) ;   

//from MATLAB : tangent = 
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14,    d15]
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24,    d25]
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34,    d35]
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14,  z*d15]
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24,  z*d25]
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34,  z*d35]
//[       d41,           d42,           d43,        -d41*z,        -d42*z,        -d43*z,    d44,    d45]
//[       d51,           d52,           d53,        -d51*z,        -d52*z,        -d53*z,    d54,    d55]

      //row 1
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14,    d15]
  tangent(0, 0) += dd(0, 0);
  tangent(0, 1) += dd(0, 1);
  tangent(0, 2) += dd(0, 2);
  tangent(0, 3) += -z * dd(0, 0);
  tangent(0, 4) += -z * dd(0, 1);
  tangent(0, 5) += -z * dd(0, 2);
  tangent(0, 6) += dd(0, 3);
  tangent(0, 7) += dd(0, 4);

  //row 2
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24,    d25]
  tangent(1, 0) += dd(1, 0);
  tangent(1, 1) += dd(1, 1);
  tangent(1, 2) += dd(1, 2);
  tangent(1, 3) += -z * dd(1, 0);
  tangent(1, 4) += -z * dd(1, 1);
  tangent(1, 5) += -z * dd(1, 2);
  tangent(1, 6) += dd(1, 3);
  tangent(1, 7) += dd(1, 4);

  //row 3
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34,    d35]
  tangent(2, 0) += dd(2, 0);
  tangent(2, 1) += dd(2, 1);
  tangent(2, 2) += dd(2, 2);
  tangent(2, 3) += -z * dd(2, 0);
  tangent(2, 4) += -z * dd(2, 1);
  tangent(2, 5) += -z * dd(2, 2);
  tangent(2, 6) += dd(2, 3);
  tangent(2, 7) += dd(2, 4);

  //row 4
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14,  z*d15]
  tangent(3, 0) += z * dd(0, 0);
  tangent(3, 1) += z * dd(0, 1);
  tangent(3, 2) += z * dd(0, 2);
  tangent(3, 3) += -z * z * dd(0, 0);
  tangent(3, 4) += -z * z * dd(0, 1);
  tangent(3, 5) += -z * z * dd(0, 2);
  tangent(3, 6) += z * dd(0, 3);
  tangent(3, 7) += z * dd(0, 4);

  //row 5
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24,  z*d25]
  tangent(4, 0) += z * dd(1, 0);
  tangent(4, 1) += z * dd(1, 1);
  tangent(4, 2) += z * dd(1, 2);
  tangent(4, 3) += -z * z * dd(1, 0);
  tangent(4, 4) += -z * z * dd(1, 1);
  tangent(4, 5) += -z * z * dd(1, 2);
  tangent(4, 6) += z * dd(1, 3);
  tangent(4, 7) += z * dd(1, 4);

  //row 6
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34,  z*d35]
  tangent(5, 0) += z * dd(2, 0);
  tangent(5, 1) += z * dd(2, 1);
  tangent(5, 2) += z * dd(2, 2);
  tangent(5, 3) += -z * z * dd(2, 0);
  tangent(5, 4) += -z * z * dd(2, 1);
  tangent(5, 5) += -z * z * dd(2, 2);
  tangent(5, 6) += z * dd(2, 3);
  tangent(5, 7) += z * dd(2, 4);

  //row 7
//[  d41,    d42,    d43, -d41*z, -d42*z, -d43*z,  d44,  d45]
  tangent(6, 0) += dd(3, 0);
  tangent(6, 1) += dd(3, 1);
  tangent(6, 2) += dd(3, 2);
  tangent(6, 3) += -z * dd(3, 0);
  tangent(6, 4) += -z * dd(3, 1);
  tangent(6, 5) += -z * dd(3, 2);
  tangent(6, 6) += dd(3, 3);
  tangent(6, 7) += dd(3, 4);

  //row 8 
//[  d51,    d52,    d53, -d51*z, -d52*z, -d53*z,  d54,  d55]
  tangent(7, 0) += dd(4, 0);
  tangent(7, 1) += dd(4, 1);
  tangent(7, 2) += dd(4, 2);
  tangent(7, 3) += -z * dd(4, 0);
  tangent(7, 4) += -z * dd(4, 1);
  tangent(7, 5) += -z * dd(4, 2);
  tangent(7, 6) += dd(4, 3);
  tangent(7, 7) += dd(4, 4);

} //end for i
 */
  

  return tangent ;
}


//print out data
void  CompositeShellSectionThermal::Print( OPS_Stream &s, int flag )
{
  s << "CompositeSectionThermal tag: " << this->getTag() << endln ; 
 // s << "Total thickness h = " << h << endln ;


}

int 
CompositeShellSectionThermal::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;
  /*
  int dataTag = this->getDbTag();

  static ID iData(3);
  iData(0) = this->getTag();
  iData(1) = nLayers;

  res += theChannel.sendID(dataTag, commitTag, iData);
  if (res < 0) {
    opserr << "WARNING CompositeShellSectionThermal::sendSelf() - " << this->getTag() << " failed to send data" << endln;
    return res;
  }
  
  if (nLayers > 0)
  {
    Vector vecData(2*nLayers+1);
    int i;
    for (i = 0; i < nLayers; i++) {
      vecData(i)         = sg[i];
      vecData(i+nLayers) = wg[i];
    }
    vecData(2*nLayers) = h;
    res += theChannel.sendVector(dataTag, commitTag, vecData);
    if (res < 0) {
      opserr << "WARNING CompositeShellSectionThermal::sendSelf() - " << this->getTag() << " failed to send data" << endln;
      return res;
    }
    
    // Send the ids of its materials
    
    int matDbTag;
    ID idData(nLayers*2);
    for (i = 0; i < nLayers; i++) {
      idData(i) = theFibers[i]->getClassTag();
      matDbTag = theFibers[i]->getDbTag();
      // ensure that the material has a database tag
      if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        if (matDbTag != 0) theFibers[i]->setDbTag(matDbTag);
      }
      idData(i+nLayers) = matDbTag;
    }
    
    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
      opserr << "WARNING CompositeShellSectionThermal::sendSelf() - " << this->getTag() << " failed to send ID" << endln;
      return res;
    }
    
    // Finally, quad asks its material objects to send themselves
    for (i = 0; i < nLayers; i++) {
      res += theFibers[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "WARNING CompositeShellSectionThermal::sendSelf() - " << this->getTag() << " failed to send its Material" << endln;
        return res;
      }
    }
  }

  */
  
  return res;
}


int 
CompositeShellSectionThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  /*
  int dataTag = this->getDbTag();

  static ID iData(3);
  res += theChannel.recvID(dataTag, commitTag, iData);

  if (res < 0) {
    opserr << "WARNING CompositeShellSectionThermal::recvSelf() - " << this->getTag() << " failed to receive data" << endln;
   return res;
  } 

  this->setTag(iData(0));
  
  int i;
  if (nLayers != iData(1))
  {
    nLayers = iData(1);
    if (sg != 0) delete sg;
    sg = new double[nLayers];
    if (wg != 0) delete sg;
    wg = new double[nLayers];
    if (theFibers !=0)
    {
      for ( i = 0; i < nLayers; i++ )
      {
        if (theFibers[i] != 0) delete theFibers[i] ;
      }
      delete [] theFibers;
    }
    theFibers = new NDMaterial*[nLayers];
  }

  if (nLayers > 0)
  {
    Vector vecData(2*nLayers+1);
    res += theChannel.recvVector(dataTag, commitTag, vecData);
    if (res < 0) {
    opserr << "WARNING CompositeShellSectionThermal::recvSelf() - " << this->getTag() << " failed to receive data" << endln;
    return res;
    }  
    for (i = 0; i < nLayers; i++) {
      sg[i] = vecData[i];
      wg[i] = vecData[i+nLayers];
    }
    h = vecData[2*nLayers];
    ID idData(nLayers*2);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
      opserr << "WARNING CompositeShellSectionThermal::recvSelf() - " << this->getTag() << " failed to receive ID" << endln;
      return res;
    }

    for (i = 0; i < nLayers; i++) {
      int matClassTag = idData(i);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theFibers[i]->getClassTag() != matClassTag) {
        if (theFibers[i] != 0) delete theFibers[i];
        theFibers[i] = theBroker.getNewNDMaterial(matClassTag);
        if (theFibers[i] == 0) {
          opserr << "CompositeShellSectionThermal::recvSelf() - " << 
            "Broker could not create NDMaterial of class type" << matClassTag << endln;
         return -1;
        }
      }
      theFibers[i]->setDbTag(idData(i+nLayers));
      // Receive the material
      res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "CompositeShellSectionThermal::recvSelf() - material " << 
          i << ", failed to recv itself" << endln;
        return res;
      }
    }
  }
  */
  
    
  return res;
}
 



double
CompositeShellSectionThermal::determineFiberTemperature(const Vector& DataMixed, double fiberLoc) 
{
		double FiberTemperature = 0;

		double dataTempe[18]; 
		for (int i = 0; i < 18; i++) { 
			dataTempe[i] = DataMixed(i);
		}

		return FiberTemperature;
}



Response*
CompositeShellSectionThermal::setResponse(const char** argv, int argc,
    OPS_Stream& output)
{
    const ID& type = this->getType();
    int typeSize = this->getOrder();

    Response* theResponse = 0;

    output.tag("SectionOutput");
    output.attr("secType", this->getClassType());
    output.attr("secTag", this->getTag());

    // deformations
    if (strcmp(argv[0], "deformations") == 0 || strcmp(argv[0], "deformation") == 0) {
        output.tag("ResponseType", "eps11");
        output.tag("ResponseType", "eps22");
        output.tag("ResponseType", "gamma12");
        output.tag("ResponseType", "theta11");
        output.tag("ResponseType", "theta22");
        output.tag("ResponseType", "theta33");
        output.tag("ResponseType", "gamma13");
        output.tag("ResponseType", "gamma23");
        theResponse = new MaterialResponse(this, 1, this->getSectionDeformation());
        // forces
    }
    else if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0) {
        output.tag("ResponseType", "p11");
        output.tag("ResponseType", "p22");
        output.tag("ResponseType", "p12");
        output.tag("ResponseType", "m11");
        output.tag("ResponseType", "m22");
        output.tag("ResponseType", "m12");
        output.tag("ResponseType", "q1");
        output.tag("ResponseType", "q2");
        theResponse = new MaterialResponse(this, 2, this->getStressResultant());

        // force and deformation
    }
    else if (strcmp(argv[0], "forceAndDeformation") == 0) {
        output.tag("ResponseType", "eps11");
        output.tag("ResponseType", "eps22");
        output.tag("ResponseType", "gamma12");
        output.tag("ResponseType", "theta11");
        output.tag("ResponseType", "theta22");
        output.tag("ResponseType", "theta33");
        output.tag("ResponseType", "gamma13");
        output.tag("ResponseType", "gamma23");
        output.tag("ResponseType", "p11");
        output.tag("ResponseType", "p22");
        output.tag("ResponseType", "p12");
        output.tag("ResponseType", "m11");
        output.tag("ResponseType", "m22");
        output.tag("ResponseType", "m12");
        output.tag("ResponseType", "q1");
        output.tag("ResponseType", "q2");
        theResponse = new MaterialResponse(this, 4, Vector(2 * this->getOrder()));
    }
    else if (strcmp(argv[0], "comp") == 0 || strcmp(argv[0], "Comp") == 0) {
        if (argc < 3) {
            opserr << "CompositeShellSectionThermal::setResponse() - need to specify more data\n";
            return 0;
        }
        int pointNum = atoi(argv[1]);

        if (pointNum <= 0) {
            pointNum = 1;
        }
        else if (pointNum > 2) {
            pointNum = 2;
        }
        output.tag("FiberOutput");
        output.attr("number", pointNum);
        if(pointNum==1)
            theResponse = theSection1->setResponse(&argv[2], argc - 2, output);
        else if(pointNum == 2)
            theResponse = theSection2->setResponse(&argv[2], argc - 2, output);
        
        output.endTag();
    }
    output.endTag(); // SectionOutput
    return theResponse;
}

int
CompositeShellSectionThermal::getResponse(int responseID, Information& secInfo)
{
    switch (responseID) {
    case 1:
        return secInfo.setVector(this->getSectionDeformation());

    case 2:
        return secInfo.setVector(this->getStressResultant());

    case 4: {
        Vector& theVec = *(secInfo.theVector);
        const Vector& e = this->getSectionDeformation();
        const Vector& s = this->getStressResultant();
        for (int i = 0; i < 8; i++) {
            theVec(i) = e(i);
            theVec(i + 8) = s(i);
        }

        return secInfo.setVector(theVec);
    }
    default:
        return -1;
    }
}