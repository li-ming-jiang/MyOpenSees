#include <SIFYWall.h>

SIFYWall::SIFYWall(int tag, int typeTag, int jt1, int jt2, int jt3, int jt4,const ID& theMemInfo):SIFMember(tag,6,typeTag, theMemInfo)
{
  //typeTag indicates composite beam or single beam
  //typeTag==1 : single beam;
  //typeTag==2 : composite beam
  ConnectedJoints->resize(4);
  (*ConnectedJoints)(0) = jt1;
  (*ConnectedJoints)(1) = jt2;
  (*ConnectedJoints)(2) = jt3;
  (*ConnectedJoints)(3) = jt4;
}




SIFYWall::~SIFYWall()
{
  //
}





SIFMember*
SIFYWall::getConnectedBeam(void)
{
  if (TypeTag==2) {
    //return connectedSlab;
  }
  else
    return 0;
}


