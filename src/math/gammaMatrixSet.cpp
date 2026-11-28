/*!
        @file    gammaMatrixSet.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#include "math/gammaMatrixSet.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "math/gammaMatrixSet_Dirac.h"
#include "math/gammaMatrixSet_Chiral.h"

bool GammaMatrixSet::init_factory()
{
  bool result = true;

  result &= GammaMatrixSet_Dirac::register_factory();
  result &= GammaMatrixSet_Chiral::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
