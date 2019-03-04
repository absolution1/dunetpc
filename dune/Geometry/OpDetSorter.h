#ifndef GEO_OPDETSORTER_H
#define GEO_OPDETSORTER_H


#include "larcorealg/Geometry/OpDetGeo.h"

namespace geo {

  //----------------------------------------------------------------------------
  inline bool sortorderOpDet(const OpDetGeo* t1, const OpDetGeo* t2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    if(xyz1[0]!=xyz2[0])
      return xyz1[0]>xyz2[0];
    else if(xyz1[2]!=xyz2[2])
      return xyz1[2]>xyz2[2];
    else
    return xyz1[1]>xyz2[1];
  }

}
#endif // GEO_OPDETSORTER_H
