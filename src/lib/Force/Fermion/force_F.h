/*!
        @file    force_F.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$

        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_INCLUDED
#define FORCE_F_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Field/field.h"

//! Base class of fermion force calculation.

/*!
    This class defines the interface of fermion force calculation.
    force_udiv() and force_udiv1() are used recursively
    to determine the smeared fermion force.       [28 Dec 2011 HM]

    - set_mode() is added. This is for the cases when the force
    - calculation is nonhermitian.                [18 Jan 2012 HM]
    - converted to an instance of template class. [17 Apr 2023 HM]
*/

typedef AForce_F<Field> Force;

#endif
