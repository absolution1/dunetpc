// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TickModTreeData_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TickModTreeData.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TickModTreeData(void *p = 0);
   static void *newArray_TickModTreeData(Long_t size, void *p);
   static void delete_TickModTreeData(void *p);
   static void deleteArray_TickModTreeData(void *p);
   static void destruct_TickModTreeData(void *p);
   static void streamer_TickModTreeData(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TickModTreeData*)
   {
      ::TickModTreeData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TickModTreeData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TickModTreeData", ::TickModTreeData::Class_Version(), "TickModTreeData.h", 10,
                  typeid(::TickModTreeData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TickModTreeData::Dictionary, isa_proxy, 16,
                  sizeof(::TickModTreeData) );
      instance.SetNew(&new_TickModTreeData);
      instance.SetNewArray(&newArray_TickModTreeData);
      instance.SetDelete(&delete_TickModTreeData);
      instance.SetDeleteArray(&deleteArray_TickModTreeData);
      instance.SetDestructor(&destruct_TickModTreeData);
      instance.SetStreamerFunc(&streamer_TickModTreeData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TickModTreeData*)
   {
      return GenerateInitInstanceLocal((::TickModTreeData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TickModTreeData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TickModTreeData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TickModTreeData::Class_Name()
{
   return "TickModTreeData";
}

//______________________________________________________________________________
const char *TickModTreeData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TickModTreeData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TickModTreeData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TickModTreeData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TickModTreeData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TickModTreeData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TickModTreeData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TickModTreeData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TickModTreeData::Streamer(TBuffer &R__b)
{
   // Stream an object of class TickModTreeData.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> run;
      R__b >> chan;
      R__b >> femb;
      R__b >> fembChan;
      R__b >> itkm;
      R__b >> nsample;
      R__b >> maxAdc;
      R__b >> maxAdc2;
      R__b >> meanAdc;
      R__b >> meanAdc2;
      R__b >> maxFraction;
      R__b >> zeroFraction;
      R__b >> oneFraction;
      R__b >> highFraction;
      R__b >> fitMean;
      R__b >> fitSigma;
      R__b >> fitExcess;
      R__b.CheckByteCount(R__s, R__c, TickModTreeData::IsA());
   } else {
      R__c = R__b.WriteVersion(TickModTreeData::IsA(), kTRUE);
      R__b << run;
      R__b << chan;
      R__b << femb;
      R__b << fembChan;
      R__b << itkm;
      R__b << nsample;
      R__b << maxAdc;
      R__b << maxAdc2;
      R__b << meanAdc;
      R__b << meanAdc2;
      R__b << maxFraction;
      R__b << zeroFraction;
      R__b << oneFraction;
      R__b << highFraction;
      R__b << fitMean;
      R__b << fitSigma;
      R__b << fitExcess;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TickModTreeData(void *p) {
      return  p ? new(p) ::TickModTreeData : new ::TickModTreeData;
   }
   static void *newArray_TickModTreeData(Long_t nElements, void *p) {
      return p ? new(p) ::TickModTreeData[nElements] : new ::TickModTreeData[nElements];
   }
   // Wrapper around operator delete
   static void delete_TickModTreeData(void *p) {
      delete ((::TickModTreeData*)p);
   }
   static void deleteArray_TickModTreeData(void *p) {
      delete [] ((::TickModTreeData*)p);
   }
   static void destruct_TickModTreeData(void *p) {
      typedef ::TickModTreeData current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TickModTreeData(TBuffer &buf, void *obj) {
      ((::TickModTreeData*)obj)->::TickModTreeData::Streamer(buf);
   }
} // end of namespace ROOT for class ::TickModTreeData

namespace {
  void TriggerDictionaryInitialization_TickModTreeData_Dict_Impl() {
    static const char* headers[] = {
"TickModTreeData.h",
0
    };
    static const char* includePaths[] = {
"/home/dladams/ups/root/v6_12_06a/Linux64bit+2.6-2.12-e15-prof/include",
"/home/dladams/dudev/dudev01/workdir/srcs/dunetpc/dune/DataPrep/Utility/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TickModTreeData_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TickModTreeData.h")))  TickModTreeData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TickModTreeData_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TickModTreeData.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TickModTreeData", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TickModTreeData_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TickModTreeData_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TickModTreeData_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TickModTreeData_Dict() {
  TriggerDictionaryInitialization_TickModTreeData_Dict_Impl();
}
