/*
        @file    run_test_MG.cpp

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (2023年03月20日 (月)) $

        @version $LastChangedRevision: 2499 $
*/

#include <stdlib.h>
#include "testlist_MG.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//====================================================================
int run_test_MG()
{
  //test_coarse::test();
  //  test_restrict::test();
  //test_prolong::test();
  //test_afopr_coarse::test_parallel();
  //test_sap_mult::test();
  test_MG_solver::test();

  return EXIT_SUCCESS;
}
