/*!
        @file    projection.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef PROJECTION_INCLUDED
#define PROJECTION_INCLUDED

#include "lib/Smear/aprojection.h"

#ifdef USE_FACTORY
#include "lib/Tools/factory.h"
#endif

//! base class for projection operator into gauge group.

/*!
  The base class for projection operator into gauge group.

 - original version               [07 Apr 2012 H.Matsufuru]
 - typedef of AProjection<Field>  [07 Feb 2023 H.Matsufuru]
 */

typedef AProjection<Field> Projection;


#endif
