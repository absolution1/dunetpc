//
// BaselineDetrend.cc - Baseline detrend tool
#include "BaselineDetrend.h"

#include "TGraphSmooth.h"
#include "TGraph.h"

#include <iostream>
#include <string>

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;

// ctor
BaselineDetrend::BaselineDetrend(fhicl::ParameterSet const& ps ) :
  m_LogLevel( ps.get<int>("LogLevel") ), 
  m_Thresh( ps.get<float>("Thresh") ), 
  m_Pad( ps.get<unsigned>("Pad") ), 
  m_MinFrac( ps.get<float>("MinFrac") ), 
  m_Span( ps.get<float>("Span") )
{
  const string myname = "BaselineDetrend::ctor: ";
  if( m_LogLevel >= 1 ){
    cout<<myname<<"Threshold to drop samples from average "<<m_Thresh<<endl;
    cout<<myname<<"Number of ticks to pad around signals above threshold "<<m_Pad<<endl;
    cout<<myname<<"Min fraction of pedestal samples for detrending "<<m_MinFrac<<endl;
    cout<<myname<<"LOWESS span parameter "<<m_Span<<endl;
  }
  
  m_GS  = new TGraphSmooth( "smoother" );
}

// dtor
BaselineDetrend::~BaselineDetrend(){
  m_GS->Delete();
}

//
DataMap BaselineDetrend::update(AdcChannelData& acd) const {
  DataMap ret;
  AdcSignalVector& data = acd.samples;

  // init trend baseline
  //m_Trend.resize( data.size(), 0. );
  
  // computed baseline trend
  auto trend = Smoother( data );

  // detrending
  for(size_t i=0;i<data.size();++i){
    data[i] -= trend[i];
  }

  //m_Trend.clear();
  // finished
  return ret;
}


//
//
// smoother
AdcSignalVector BaselineDetrend::Smoother( const AdcSignalVector &data ) const
{
  const string myname = "BaselineDetrend::Smoother: ";
  AdcSignalVector trend( data.size(), 0 );
  
  vector<unsigned> idx;
  for( unsigned i = 0; i < data.size(); ++i ){
    if( data[i] > m_Thresh ) {
      unsigned count = 1;
      // remove last elements from the vector
      while( !idx.empty() ){
	if( count >= m_Pad ) break;
	//
	if( (int)idx.back() != ((int)i - (int)count) ) break;
	idx.pop_back();
	count++;
      }

      // move forward until we are below threshold + num of samples to skip after
      count = 0;
      for( unsigned j = i+1; j < data.size(); ++j ){
	if( data[j] > m_Thresh ) continue;
	if( count >= m_Pad ) { 
	  i = j;
	  break;
	}
	count++;
      }

      continue;
    }
   
    //
    idx.push_back( i );
  }
  
  // do nothing if we have only a limited number of ped samples
  if( ((float)idx.size())/data.size() < m_MinFrac ){
    if( m_LogLevel >= 2 ){
      cout<<myname<<"number of baseline samples is too low\n";
    }
    return trend;
  }
    
  //
  TGraph gin( idx.size() );
  for( unsigned i=0; i<idx.size(); ++i ){
    gin.SetPoint(i, idx[i], data[idx[i]] );
  }
  
  // do smoothing
  auto gout  = m_GS->SmoothLowess(&gin,"", m_Span);
  
  // our graph data are sorted in ticks
  gout->SetBit(TGraph::kIsSortedX);
  
  // compute detranding function
  unsigned nidx = idx.size();
  for( unsigned i = 0; i<nidx; ++i ){
    unsigned ii = idx[i];
    double x, y;
    gout->GetPoint(i, x, y);
    trend[ii]  = y;

    //
    if( i == 0 ){
      if( ii != 0 ){
	for(unsigned j=0;j<ii;++j)
	  trend[j] = gout->Eval( j );
      }
    }
    else if( i == nidx - 1 ) {
      unsigned iii = data.size() - 1;
      if( ii != iii){
	for(unsigned j=ii;j<=iii;++j)
	  trend[j] = gout->Eval( j );
      }
    }
    else {
      unsigned iii = idx[i-1];
      if( (ii - iii ) > 1 ){
	for(unsigned j=iii;j<ii;++j)
	  trend[j] = gout->Eval( j );
      }
    }
    
  }
  
  return trend;
}
