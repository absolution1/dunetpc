////////////////////////////////////////////////////////////////////////
/// \file  MuonCounter35Alg.cxx
/// \brief Interface to algorithm class for DUNE 35t muon counters.
///
/// \version $Id:  $
/// \author  Matthew Worcester (mworcester@bnl.gov)
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/MuonCounter35Alg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include <fstream>

namespace geo{

  //----------------------------------------------------------------------------

  MuonCounter35Alg::MuonCounter35Alg(fhicl::ParameterSet const& p)
  {

  }

  //----------------------------------------------------------------------------

  MuonCounter35Alg::~MuonCounter35Alg(){}

  //----------------------------------------------------------------------------

  int MuonCounter35Alg::loadMuonCounterGeometry(char* filename, 
						std::vector< std::vector<double> > &geometry){

    // This function loads the muon counter geometry from an external text file
    // into a vector of vectors. Each sub-vector holds the geometry for one counter
    // with the format: counter ID, condition flag (0=off, non-zero=on), list of 
    // coordinates for each corner of the counter (x0,y0,z0,x1,y1,z1,...).
    // It returns 1 if some geometry data loads with the correct size, and 0 if 
    // it fails. The text file should contain one counter per line with exactly
    // the same format as the geometry vectors.

    int loaded = 1;

    // open the file
    std::ifstream counterfile(filename);

    // check that the file opened
    if(counterfile.is_open()){

      mf::LogInfo("MuonCounter35t") << "Loading 35t muon counter geometry from " << filename;
    
      // load the counter geometry

      std::string line;
      int ncounter = 0;

      // loop over the lines in the text file, each line represents a counter
      while(std::getline(counterfile,line)){

	// skip lines that begin with #
	if(line.at(0) != '#'){

	  // this is a counter, so add an empty vector to fill with the geometry
	  geometry.push_back(std::vector<double>());
	  ncounter = geometry.size()-1;

	  // format the line to easily extract the geometry
	  std::istringstream issline(line);
	  double value = 0.;

	  while(issline >> value){

	    // add the geometry to the vector
	    geometry[ncounter].push_back(value);

	  }

	}

      } // done looping over lines in the text file

      counterfile.close();

    }
    else{

      // the text file did not open

      mf::LogError("MuonCounter35t") << "Could not open " << filename;

      loaded = 0;

    }

    // check that some geometry data made it into the vectors

    if(geometry.size() != 0){

      // check that the amount of geometry data in each counter is correct

      for(unsigned int i=0; i<geometry.size(); i++){

	// there must be at least 3 verticies, so 9 values + 2 for the counter ID and flag

	if(geometry[i].size() < 11){

	  //std::cout << "ERROR: counter " << geometry[i][0] << " has size " << geometry[i].size() 
	  //	  << ", when it must have size > 11" << std::endl;

	  loaded = 0;

	}

	// there must be x3 + 2 values

	if((geometry[i].size()-2) % 3 != 0){

	  //std::cout << "ERROR: counter " << geometry[i][0] << " has size " << geometry[i].size() 
	  //	  << ", when it must have size = x3 + 2" << std::endl;

	  loaded = 0;

	}

      } // done looping over counter vector

    }
    else{

      // no geometry data made it into the vectors

      loaded = 0;

    }

    return loaded;

  }
  
  //----------------------------------------------------------------------------

  int MuonCounter35Alg::testTrackInAllCounters(int trackID, 
					       TVector3 trackpoint, TVector3 trackdirection, 
					       std::vector< std::vector<double> > &geometry, 
					       std::vector< std::vector<double> > &hitcounters){

    // This function tests whether a track, defined by a point and a direction, 
    // intersects with any of the muon counters, defined in geometry. The 
    // geometry defines each counter with an ID number, a condition flag (0=off), 
    // and an area defined by one coordinate set for each corner of the counter 
    // (this is to handle trapezoidal counters and other oddly-shaped bits).
    // It fills a vector of vectors, where each sub-vector contains the ID and 
    // condition flag of a hit counter, and the track ID and intersection point with
    // that counter. It returns the number of hit counters.

    int counterID = -1;
    int counterFlag = -1;

    // loop over the counters
    for(unsigned int igeo=0; igeo<geometry.size(); igeo++){

      // first two indicies of each counter geometry are ID and flag
      counterID = geometry[igeo][0];
      counterFlag = geometry[igeo][1];

      // check that the counter is on
      if(counterFlag != 0){

	// check if the track intersects the counter and get the intersection point
	
	TVector3 insectionPoint(0.,0.,0.);

	int inside = testTrackInCounter(trackpoint, trackdirection, geometry[igeo], insectionPoint);

	// if track intersection point is inside the counter area, fill the 
	// vector of hit counter data

	if(inside){

	  // add an empty vector to fill with the hit data
	  hitcounters.push_back(std::vector<double>());
	  int ncounter = hitcounters.size()-1;

	  // add the data to the vector
	  hitcounters[ncounter].push_back(counterID);
	  hitcounters[ncounter].push_back(counterFlag);
	  hitcounters[ncounter].push_back(trackID);
	  hitcounters[ncounter].push_back(insectionPoint.X());
	  hitcounters[ncounter].push_back(insectionPoint.Y());
	  hitcounters[ncounter].push_back(insectionPoint.Z());

	} // check if a track/counter plane intersection point was found

      } // check that the counter condition flag is non-zero

    } // loop over the counters

    return hitcounters.size();

  }

  //----------------------------------------------------------------------------

  int MuonCounter35Alg::testTrackInCounter(TVector3 trackpoint, TVector3 trackdirection, 
					   std::vector<double> &singlecountergeometry,
					   TVector3 &intersectionpoint){

    // This function tests whether a track, defined by a point and a direction, 
    // intersects with a muon counter, defined in singlecountergeometry. The 
    // singlecountergeometry contains an ID number, a condition flag (0=off), 
    // and an area defined by one coordinate set for each corner of the counter 
    // (this is to handle trapezoidal counters and other oddly-shaped bits).
    // It fills a TVector3 with the coordinates of the intersection point
    // between the track and counter. It returns 1 if the track intersects the 
    // counter and 0 otherwise.

    int intersect = 0;

    // sort the counter verticies into TVector3 vectors
    std::vector<TVector3> points;

    // after the first two indicies, each next three values are coordinates

    for(unsigned int icoord=2; icoord<singlecountergeometry.size(); icoord=icoord+3){

      TVector3 temppoint(singlecountergeometry[icoord],singlecountergeometry[icoord+1],
			 singlecountergeometry[icoord+2]);

      points.push_back(temppoint);

    }

    //std::cout << "Npoints " << points.size() << " first value " << points[0].X() 
    //	    << " last value " << points[points.size()-1].Z() << std::endl;

    // calculate the coefficients of the plane defined by the first 3 of those points

    double a = points[0].Y()*(points[1].Z()-points[2].Z())
      + points[1].Y()*(points[2].Z()-points[0].Z()) 
      + points[2].Y()*(points[0].Z()-points[1].Z());

    double b = points[0].Z()*(points[1].X()-points[2].X())
      + points[1].Z()*(points[2].X()-points[0].X())
      + points[2].Z()*(points[0].X()-points[1].X());

    double c = points[0].X()*(points[1].Y()-points[2].Y())
      + points[1].X()*(points[2].Y()-points[0].Y())
      + points[2].X()*(points[0].Y()-points[1].Y());

    double d = - points[0].X()*(points[1].Y()*points[2].Z() - points[2].Y()*points[1].Z())
      - points[1].X()*(points[2].Y()*points[0].Z() - points[0].Y()*points[2].Z())
      - points[2].X()*(points[0].Y()*points[1].Z() - points[1].Y()*points[0].Z());

    //std::cout << "cooefficients: a " << a << " b " << b << " c " << c << " d " << d << std::endl;

    // for completeness, check that the rest of the points are also in the plane 

    for(unsigned int p=3; p<points.size(); p++){

      double plane = a*points[p].X() + b*points[p].Y() + c*points[p].Z() + d;

      if(plane > 0.0001){

	//std::cout << "ERROR: the corners of muon counter " << singlecountergeometry[0]
	//	<< " do not lie in a plane." << std::endl;

      }

    }

    // find the normal to the plane of the counter

    double mag = sqrt(a*a + b*b + c*c);

    if(mag != 0){

      a = a/mag;
      b = b/mag;
      c = c/mag;

    }

    // find the intersection, if any, between the given track and the plane of the 
    // counter: start with any point in the plane and the normal to the plane

    TVector3 planepoint = points[2];
    TVector3 planenorm(a, b, c);

    // distance between a point in the plane of the counter and a point on the track

    TVector3 diff = planepoint - trackpoint;
    //std::cout << "diff " << diff[0] << ", " << diff[1] << ", " << diff[2] << std::endl;

    // pick out the distance between the track and plane in the direction normal 
    // to the plane; if num = 0, then the track point itself is in the plane of the counter

    double num = diff * planenorm;
    //std::cout << "num " << num << std::endl;

    // pick out the vector of the track normal to the plane of the counter
	
    double den = trackdirection * planenorm;
    //std::cout << "den " << den << std::endl;

    // if den = 0, then the track direction is perpendicular to the plane of the counter
    // and there are no intersections; except in the special case that num = 0, in which
    // case the entire track is in the plane of the counter, and the intersection
    // point is chosen to be the initial track point

    double r = 0.;
    int found_plane_intersection = 0;

    // calculate the intersection point between the track and plane of the counter

    TVector3 planeinterpoint(0.,0.,0.);

    if(fabs(den) > 0.0001){

      // case 1: the track direction is not perpendicular to the plane of the counter
      // so a single intersection point should be found; if num=0 then the track 
      // point itself is in the plane of the counter

      r = num/den;
      //std::cout << "r " << r << std::endl;
	      
      planeinterpoint = r*trackdirection + trackpoint;
      //std::cout << "planeinterpoint " << planeinterpoint.X() << ", " << planeinterpoint.Y()
      //	      << ", " << planeinterpoint.Z() << std::endl;

      found_plane_intersection = 1;

    }
    else{

      if(num == 0){

	// case 2: the entire track is in the plane of the counter; so we should 
	// check if the track crosses the counter area

      }

    }
	    
    // check if the plane intersection point is inside the counter

    if(found_plane_intersection){

      // find the coordinate plane that the plane of the counter is closest to 

      double xdiff = 0.; double xdiffmax = 0;
      double ydiff = 0.; double ydiffmax = 0;
      double zdiff = 0.; double zdiffmax = 0;

      for(unsigned int p=1; p<points.size(); p++){

	xdiff = points[p].X() - points[p-1].X();
	if(xdiff > xdiffmax) xdiffmax = xdiff;

	ydiff = points[p].Y() - points[p-1].Y();
	if(ydiff > ydiffmax) ydiffmax = ydiff;

	zdiff = points[p].Z() - points[p-1].Z();
	if(zdiff > zdiffmax) zdiffmax = zdiff;

      }

      //std::cout << xdiffmax << "  " << ydiffmax << "  " << zdiffmax << std::endl;

      // Project both the coordinates of the counter and the intersection point
      // down into the coordinate that the plane of the counter is closest to.
      // This gives us a 2D polygon and a point which we can test to see if 
      // it is inside the polygon.

      double *vert1 = NULL;
      vert1 = (double*) realloc(vert1,points.size()*sizeof(double));

      double *vert2 = NULL;
      vert2 = (double*) realloc(vert2,points.size()*sizeof(double));

      if(vert1 != NULL && vert2 != NULL){

	for(unsigned int i=0; i<points.size(); i++){

	  vert1[i] = 0;
	  vert2[i] = 0;

	}

      }

      double test1 = 0.;
      double test2 = 0.;

      if(xdiffmax < ydiffmax && xdiffmax <= zdiffmax){
		    
	// plane of the counter is closest to the y-z plane

	for(unsigned int i=0; i<points.size(); i++){

	  vert1[i] = points[i].Y();
	  vert2[i] = points[i].Z();

	}

	test1 = planeinterpoint.Y();
	test2 = planeinterpoint.Z();

      }
      else if(ydiffmax < xdiffmax && ydiffmax <= zdiffmax){

	// plane of the counter is closest to the x-z plane

	for(unsigned int i=0; i<points.size(); i++){

	  vert1[i] = points[i].X();
	  vert2[i] = points[i].Z();

	}

	test1 = planeinterpoint.X();
	test2 = planeinterpoint.Z();

      }
      else if (zdiffmax < xdiffmax && zdiffmax < ydiffmax){

	// plane of the counter is closest to the x-y plane

	for(unsigned int i=0; i<points.size(); i++){

	  vert1[i] = points[i].X();
	  vert2[i] = points[i].Y();

	}

	test1 = planeinterpoint.X();
	test2 = planeinterpoint.Y();

      }
      else{

	//std::cout << "ERROR: case not considered for projecting counter. " 
	//	<< "Testing the intersection point will fail." << std::endl;

      }

      // test for the intersection point inside the counter polygon
      // 1 = inside, 0 = outside

      intersect = testPointInPolygon(points.size(), vert1, vert2, test1, test2);

      free(vert1);
      free(vert2);

      if(intersect) intersectionpoint = planeinterpoint;

      //std::cout << "intersection? " << intersect << std::endl;

    } // found an intersection between track and plane of counter

    return intersect;

  }
 
  //----------------------------------------------------------------------------

  int MuonCounter35Alg::testPointInPolygon(int nvert, double *vertx, double *verty, 
					   double testx, double testy){

    // This function checks whether a given point (testx, testy) is inside a polygon
    // defined by 2D coordinate arrays (vertx and verty) of the polygon's
    // verticies. If the point is inside the polygon, the function will return 1; if
    // outside, it returns 0; if the point is exactly on one edge, 
    // (e.g. testx = vertx[1] = vertx[2]) it is treated as outside and returns 0.
    // Original author: W. Randolph Franklin. For the license to use, see the 
    // source: www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html 

    int i, j, c = 0;

    for(i = 0, j = nvert - 1; i < nvert; j = i++){

      if( ((verty[i] > testy) != (verty[j] > testy)) &&
	  (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ){

	c = !c;

      }

    }

    return c;

  }

  //----------------------------------------------------------------------------

} // namespace
