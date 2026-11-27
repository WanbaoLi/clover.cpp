/*!
        @file    forceSmear.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCESMEAR_INCLUDED
#define FORCESMEAR_INCLUDED

#include "lib/Smear/aforceSmear.h"
#include "lib/Field/field.h"
#include "lib/IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "lib/Tools/factory.h"
#endif


//! Base class for force calculation of smeared operators.

/*!
    Base class for force calculation of smeared operators.
                                     [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks.
                                     [21 Mar 2015 Y.Namekawa]
    typedef of template base class.
                                     [21 Mar 2023 H.Matsufuru]
 */

typedef AForceSmear<Field> ForceSmear;

#endif
