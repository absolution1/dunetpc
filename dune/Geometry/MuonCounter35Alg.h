////////////////////////////////////////////////////////////////////////
/// \file  MuonCounter35Alg.h
/// \brief Interface to algorithm class for DUNE 35t muon counters.
///
/// \version $Id:  $
/// \author  Matthew Worcester (mworcester@bnl.gov)
////////////////////////////////////////////////////////////////////////
#ifndef GEO_MUONCOUNTER35ALG_H
#define GEO_MUONCOUNTER35ALG_H

#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class MuonCounter35Alg : public ChannelMapAlg{

  public:

    // constructor and destructor

    explicit MuonCounter35Alg(fhicl::ParameterSet const& p);

    virtual ~MuonCounter35Alg();
    
    // read in muon counter geometry

    static int loadMuonCounterGeometry(char* filename, std::vector< std::vector<double> > &geometry);

    // test if a track intersects any counter

    static int testTrackInAllCounters(int trackID, TVector3 trackpoint, TVector3 trackvector, 
			       std::vector< std::vector<double> > &geometry,
			       std::vector< std::vector<double> > &hitcounters);

    // test if a track intersects a single counter

    static int testTrackInCounter(TVector3 trackpoint, TVector3 trackvector,
			   std::vector<double> &singlecountergeometry,
			   TVector3 &intersectionpoint);

    // function to test if a point is inside a 2D polygon

    static int testPointInPolygon(int nvert, double *vertx, double *verty, double testx, double testy);

  private:
    


  };

}
#endif // GEO_MUONCOUNTER35ALG_H

