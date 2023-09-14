// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME vectorMatrix_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "RootTrajectorySummaryWriter.hpp"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *vectorlETMatrixTlEdoublegRsPgR_Dictionary();
   static void vectorlETMatrixTlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlETMatrixTlEdoublegRsPgR(void *p = nullptr);
   static void *newArray_vectorlETMatrixTlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlETMatrixTlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlETMatrixTlEdoublegRsPgR(void *p);
   static void destruct_vectorlETMatrixTlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TMatrixT<double> >*)
   {
      vector<TMatrixT<double> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TMatrixT<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TMatrixT<double> >", -2, "vector", 389,
                  typeid(vector<TMatrixT<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETMatrixTlEdoublegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TMatrixT<double> >) );
      instance.SetNew(&new_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlETMatrixTlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TMatrixT<double> > >()));

      ::ROOT::AddClassAlternate("vector<TMatrixT<double> >","std::vector<TMatrixT<double>, std::allocator<TMatrixT<double> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TMatrixT<double> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETMatrixTlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TMatrixT<double> >*>(nullptr))->GetClass();
      vectorlETMatrixTlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETMatrixTlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TMatrixT<double> > : new vector<TMatrixT<double> >;
   }
   static void *newArray_vectorlETMatrixTlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TMatrixT<double> >[nElements] : new vector<TMatrixT<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      delete (static_cast<vector<TMatrixT<double> >*>(p));
   }
   static void deleteArray_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      delete [] (static_cast<vector<TMatrixT<double> >*>(p));
   }
   static void destruct_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      typedef vector<TMatrixT<double> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TMatrixT<double> >

namespace {
  void TriggerDictionaryInitialization_vectorMatrix_dict_Impl() {
    static const char* headers[] = {
"RootTrajectorySummaryWriter.hpp",
nullptr
    };
    static const char* includePaths[] = {
"/home/niko/git/acts/Core/include/",
"/home/niko/anaconda3/envs/acts/include/eigen3/",
"/home/niko/git/acts/Examples/Framework/include/",
"/home/niko/git/acts/Fatras/include/",
"/home/niko/git/acts/Plugins/FpeMonitoring/include/",
"/home/niko/anaconda3/envs/acts/include/",
"/home/niko/git/acts/Examples/Io/Root/include/ActsExamples/Io/Root/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "vectorMatrix_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
template <class Element> class __attribute__((annotate("$clingAutoload$TMatrixT.h")))  __attribute__((annotate("$clingAutoload$RootTrajectorySummaryWriter.hpp")))  TMatrixT;

namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "vectorMatrix_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "RootTrajectorySummaryWriter.hpp"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("vectorMatrix_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_vectorMatrix_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_vectorMatrix_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_vectorMatrix_dict() {
  TriggerDictionaryInitialization_vectorMatrix_dict_Impl();
}
