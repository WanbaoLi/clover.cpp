/*!
        @file    fopr_Wilson_eo.h

        @brief

        @author  UEDA, Satoru
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FOPR_WILSON_EO_INCLUDED
#define FOPR_WILSON_EO_INCLUDED

//! Even-odd Wilson fermion operator.

/*!
    This class is an even-odd version of Wilson fermion operator.
    At present this is rough implementation, while correctly
    works, and to be updated by supplying complete functionality.
    Only the functions needed for even-odd preconditioned solver
    is ready.
                                     [20 Jun 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    Selector is implemented.         [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Implementation is separated to Fopr_Wilson_eo_impl class.
                                     [06 Jul 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */

#include "Imp/fopr_Wilson_eo_impl.h"

using Fopr_Wilson_eo = Imp::Fopr_Wilson_eo;

#endif
