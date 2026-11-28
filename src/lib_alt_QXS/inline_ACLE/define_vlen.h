/*!
        @file    define_vlen.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

/*!
    This file provides definitions of VLEN for QXS branch
    of alternative code in Bridge++.
                                 [H.Matsufuru 23 Feb 2021]
*/

#ifndef QXS_DEFINE_VLEN_INCLUDED
#define QXS_DEFINE_VLEN_INCLUDED

#define VLEND_QXS    8
#define VLENS_QXS    16


#define VLENXD       4
#define VLENYD       2

#define VLENXS       4
#define VLENYS       4

#define VLENS    VLENS_QXS
#define VLEND    VLEND_QXS

#if (VLEND != VLENXD * VLENYD)
#error bad VLENXD * VLENYD
#endif

#if (VLENS != VLENXS * VLENYS)
#error bad VLENXD * VLENYD
#endif

#endif
