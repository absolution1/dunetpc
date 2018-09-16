////////////////////////////////////////////////////////////////////////
//
// File: CalibTreeRecord.h
// Author: Jason Stock (jason.stock@mines.sdsmt.edu)
// This file is a set of classes for writeout of data relative to calibration sources.
// The main purpose of this class is to allow easy fast analysis of calibration
// simulations.
//
////////////////////////////////////////////////////////////////////////

//includes
#ifndef DUNE_DUNEOBJ_CALIBTREERECORD_H
#define DUNE_DUNEOBJ_CALIBTREERECORD_H
#include <vector>
#include <TObject.h>
#include <iostream>

namespace CalibTreeRecord {
  typedef int64_t eve_number_t ;
  typedef int64_t particle_number_t ;
  //class OpHit : public TObject{

  class HCRec{
    /*ClassDef(CalibTreeRecord::HCRec,10);*/
    public:
    int eve_index;
    int part_index;
    int hit_index;
    HCRec(int e, int p, int h)
      :eve_index(e),part_index(p),hit_index(h)
    {}
    HCRec()
      :eve_index(0), part_index(0), hit_index(0)
    {}
  };

  class HitContributor {
    /*ClassDef(CalibTreeRecord::HitContributor,10);*/
    public:
      int index; //The internal hit index of the hit this record is describing.
      std::vector<HCRec> locations; //eve index, part index, hits index (vector position).
      HitContributor(int i)
        :index(i), locations()
      {}
      HitContributor()
        :index(), locations()
      {}
  };

  class PartialPDOpHit {
    public:
      Double_t pes;
      Double_t num_photons;
      Double_t time;
      Double_t width;
      Double_t energy;
      Double_t split;//fraction of hit energy in this partial.
      UInt_t opdet;
      int64_t index;
  };

  class PartialOpHit : public PartialPDOpHit {
    /*ClassDef(CalibTreeRecord::PartialOpHit,10);*/
    public:
      UInt_t opchan;
  };

  //class Hit : public TObject{
  class PartialHit  {
    /*ClassDef(CalibTreeRecord::PartialHit,10);*/
    public:
      Double_t charge;
      Double_t num_electrons;
      Double_t energy;
      Double_t time;
      Double_t width;
      Double_t split;//fraction of track energy in this hit.
      UInt_t wire;
      int64_t index;
      Bool_t is_collection_wire;
  };

/*  class TrackRecord {
      bool isEve;
      Int_t trackId;
      Double_t x_pos;
      Double_t y_pos;
      Double_t z_pos;
      std::vector<PartialHit> partial_hits;
      std::vector<PartialOpHit> partial_ophits;
  }*/ //without useful rescattering information this is not possible.

  //class ParticleRecord :public TObject {
  class ParticleRecord  {
    /*ClassDef(CalibTreeRecord::ParticleRecord,10);*/
    public:
      bool isEve;
      //particle_number_t number;
      Int_t trackId;
      UInt_t pdgid;
      Double_t dP;
      Double_t dE;
      Double_t x_pos;
      Double_t y_pos;
      Double_t z_pos;
      Double_t t_pos;
      std::vector<PartialHit> partial_hits;
      std::vector<PartialOpHit> partial_ophits;
  };

  class EveRecord{
    /*ClassDef(CalibTreeRecord::EveRecord,10);*/
    public:
      //eve_number_t number;
      Int_t trackId;
      UInt_t pdgid;
      Double_t x_pos;
      Double_t y_pos;
      Double_t z_pos;
      Double_t t_pos;
      std::vector<ParticleRecord> particles;
      std::string generator="";

  };

  //class CalibTreeRecord :public TObject {
  class CalibTreeRecord {
    /*ClassDef(CalibTreeRecord::CalibTreeRecord,10);*/
    public:
      bool bunch_hits;
      UInt_t run;
      UInt_t subrun;
      UInt_t event;
      std::vector<EveRecord> eves;
      std::vector<HitContributor> hits;
      std::vector<HitContributor> ophits;

      void Clear();
      void stdout_dump();
  };



}

#endif
