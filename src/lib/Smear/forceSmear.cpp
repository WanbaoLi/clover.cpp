/*!
        @file    forceSmear.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/forceSmear.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "lib/Smear/forceSmear_APE.h"
#include "lib/Smear/forceSmear_HYP.h"

template<>
bool ForceSmear::init_factory()
{
  bool result = true;

  result &= ForceSmear_APE::register_factory();
  result &= ForceSmear_HYP::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
