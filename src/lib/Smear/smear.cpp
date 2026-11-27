/*!
        @file    smear.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/smear.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "lib/Smear/smear_APE.h"
#include "lib/Smear/smear_APE_spatial.h"
// pruned: #include "lib/Smear/smear_HYP.h"
// pruned: #include "lib/Smear/smear_APE_SF.h"
// pruned: #include "lib/Smear/smear_HYP_SF.h"

template<>
bool Smear::init_factory()
{
  bool result = true;

  result &= Smear_APE::register_factory();
  result &= Smear_APE_spatial::register_factory();
// pruned:   result &= Smear_HYP::register_factory();
// pruned:   result &= Smear_APE_SF::register_factory();
// pruned:   result &= Smear_HYP_SF::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
