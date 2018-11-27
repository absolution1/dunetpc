//File: CRTID.h
//Brief: Tools for sorting CRT channels by detector groupings.  
//       Organizing CRT channels this way should help users 
//       identify overlapping channels for reconstruction.  
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRTID_H
#define CRTID_H

//c++ includes
#include <cstdint> //For uint8_t and friends
#include <cstddef> //For size_t
#include <map> //For CRT::map definition

namespace CRT
{
  //Programming model:
  //
  //I want to express the CRT geometry as a series of abstract layers 
  //that can be better defined based on what I know about the geometry 
  //in a particular application for particular run conditions.  To keep 
  //things as simple as possible, I'm defining "key" classes to uniquely 
  //identify layers of the geometry while letting users decide what 
  //they want to do with those layers.  I'm planning to use these 
  //identifiers with std::map, but you can probably find other ways to be 
  //even more creative.  
  //
  //The top-level geometry concept isa frame.  Inside of each frame are 
  //2 planes of modules in different orientations.  So, to specify a plane, 
  //one also has to specify what frame that plane is a part of.  There are 5 
  //layers like this in the CRT geometry: Frame, Plane, Module, Layer, and Strip.  
  //I have provided a typedef for each layer below so that you can "effortlessly" 
  //use these concepts in your code without ever knowing that I'm using class 
  //template.  

  //Implementation details:
  //
  //The original design of this scheme, where I wrote a separate struct 
  //for each abstract layer and hoped I didn't type any mistakes, failed 
  //miserably because it was very hard to maintain when I came across my 
  //first bug.  So, I made the compiler do (almost) all of the work for me! 
  //
  //Essentially, I want to create a tree at compile-time to sort 
  //based on multiple layers of CRT geometry abstraction.  To reduce 
  //the opportunity for typing mistakes and the number of points 
  //of maintenence, I'm going to write a tree node once as an ugly 
  //class template, then I'll typedef the appropriate class template 
  //instantiations to make things more user-friendly. 
  namespace detail
  {
    //Abandon all hope ye who enter here.  Somewhat-complicated class 
    //templates ahead.  If you want to read this, you should be comfortable 
    //with class templates, operator overloading, and partial template 
    //speicializations.  
    template <class BASE, class ID> 
    //BASE is the Node above this one.  I plan 
    //to write a partial template specialization 
    //to allow for a root node shortly. 
    //
    //ID is a comparable type that uniquely identifies 
    //elements within a layer of this compile-time octree.
    struct Node: public BASE
    {
      explicit Node(const BASE& base, const ID& id): BASE(base), fID(id) {}
      virtual ~Node() = default;

      //Compare Nodes that represent the same layer based on 
      //their identifier within this layer.  If Nodes are from 
      //a different layer, they should not compare equal. 
      //
      //Comparison operators are intentionally NOT virtual so 
      //that users can choose to sort strips only into modules 
      //for example.   
      bool operator ==(const Node<BASE, ID>& other) const
      {
        if(!(this->BASE::operator==(other))) return false;
        return fID == other.fID;
      }

      bool operator < (const Node<BASE, ID>& other) const
      {
        //if(this->BASE::operator<(other)) return true;
        if(!(this->BASE::operator==(other))) return this->BASE::operator<(other); 
        return fID < other.fID;
      }

      protected:
        const ID fID; //Unique identifier for an element in this layer of a compile-time octree.
    };
    
    //Specialization for BASE=void that defines a root node
    template <class ID>
    struct Node<void, ID>
    {
      explicit Node(const ID& id): fID(id) {}
      virtual ~Node() = default;

      bool operator ==(const Node<void, ID>& other) const
      {
        return fID == other.fID;
      }

      bool operator < (const Node<void, ID>& other) const
      {
        return fID < other.fID; 
      }

      protected:
        const ID fID; //Unique identifier for an element in this layer of a compile-time octree.
    };  
  }

  //A CRT Frame consists of 4 modules in 2 layers.  
  //One layer is horizontal and gives y positions, and 
  //the other layer is vertical and gives x positions.  
  //In the ProtoDUNE-SP Cosmic Ray Tagger, there are 
  //4 frames upstream of the TPC and 4 frames downstream.
  //
  //In terms of a tree, FrameID is the root node.  
  typedef detail::Node<void, uint8_t> FrameID;

  //Within each frame, there are two planes of CRT modules to give 
  //a 2D position readout.  A plane belongs to a frame, so a PlaneID 
  //is a FrameID.  
  typedef detail::Node<FrameID, bool> PlaneID;

  //Each CRT module has 2 layers of 32 scintillator strips 
  //each.  A CRT module belongs to a plane within a frame.  
  //So, a CRT::ModuleID is a PlaneID.  In the ProtoDUNE-SP 
  //CRT, there were two modules in each orientation in 
  //a frame.
  typedef detail::Node<PlaneID, bool> ModuleID;

  //A CRT module has 2 layers of 32 scintillator strips each. 
  //A layer is inside a module, so a LayerID is a ModuleID.  
  //Since modules may have been rotated during mounting, define 
  //layers based on which layer is closer to the HV connector on 
  //a CRT module's board.  The PMT can only be mounted in one 
  //direction.  
  typedef detail::Node<ModuleID, bool> LayerID;

  //There are 32 strips in each layer of a CRT module. 
  //A strip belongs to a layer, so a StripID is a LayerID.  
  //typedef detail::Node<LayerID, uint8_t> StripID;
  struct StripID: public  detail::Node<LayerID, uint8_t>
  {
    StripID(const LayerID& layer, const uint8_t local): detail::Node<LayerID, uint8_t>(layer, local) {}
    virtual ~StripID() = default;

    bool Overlaps(const StripID& other) const
    {
      if(!(this->ModuleID::operator==(other))) return false; //Strips from different modules can't overlap
      if(this->LayerID::operator==(other)) return false; //Strips from the same layer can't overlap

      const auto diff = fID - other.fID;

      //diff == 0 overlaps for either the top or the bottom.  
      if(diff == 0) return true;

      const auto layer = this->LayerID::fID;
      //If this is the bottom layer
      if(layer) return diff == 1;

      //If we survived this far, this is the top layer
      return diff == -1;
    }
  };

  //Provide a clone of std::map that has a less cumbersome interface.  We really just want to repeatedly 
  //cast a StripID to different types for different nodes' comparisons. 
  //
  //VALUE is any class that can be the second template argument to std::map.
  //Each of KEY and KEYS is an object that can be used as the first template argument to std::map.  
  template <class VALUE, class KEY, class ...KEYS>
  class map
  {
    private:
      std::map<KEY, CRT::map<VALUE, KEYS...>> fMap; //Implement this class with std::map and just give a signature 
                                                    //for operator [] that better matches what I want to do.

    public:
      map(): fMap() {}
      virtual ~map() = default;

      //Find the most-derived key type for the argument to operator[]
      using most_derived_key_type = typename CRT::map<VALUE, KEYS...>::most_derived_key_type; 

      //Public interface from std::map but with operator [] better suited to my needs.  
      //TODO: Every Map's operator[] needs to take the most-derived argument, but CRT::StripIDs need to be sorted 
      //      into less-derived ID objects.  
      inline VALUE& operator [](const most_derived_key_type& key) { return fMap[key][key]; }

      using iterator = decltype(fMap.begin());
      inline iterator begin() { return fMap.begin(); }
      inline iterator end() { return fMap.end(); }

      using const_iterator = decltype(fMap.cbegin());
      inline const_iterator begin() const { return fMap.cbegin(); }
      inline const_iterator end() const { return fMap.cend(); }

      inline size_t size() const { return fMap.size(); }
  };

  //Specialization for CRT::map with just one type as key.  Basically, expose std::map<KEY, VALUE>'s interface
  template <class VALUE, class KEY>
  class map<VALUE, KEY>
  {
    private:
      std::map<KEY, VALUE> fMap; //Implement this class with std::map and just give a signature
                                 //for operator [] that better matches what I want to do.

    public:
      map(): fMap() {}
      virtual ~map() = default;

      using most_derived_key_type = KEY;

      //Public interface from std::map but with operator [] better suited to my needs.
      inline VALUE& operator [](const KEY& key) { return fMap[key]; }

      using iterator = decltype(fMap.begin());
      inline iterator begin() { return fMap.begin(); }
      inline iterator end() { return fMap.end(); }

      using const_iterator = decltype(fMap.cbegin());
      inline const_iterator begin() const { return fMap.cbegin(); }
      inline const_iterator end() const { return fMap.cend(); }

      inline size_t size() const { return fMap.size(); }
  };

  //TODO: I could probably just use map instead of geoMap if I wrote some implementation class that splits the parameter pack above
  template <class VALUE>
  using geoMap = CRT::map<VALUE, CRT::FrameID, CRT::PlaneID, CRT::ModuleID, CRT::LayerID, CRT::StripID>;
}

#endif //CRTID_H 
