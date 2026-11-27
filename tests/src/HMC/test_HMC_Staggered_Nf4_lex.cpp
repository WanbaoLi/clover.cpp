/*!
        @file    test_HMC_Staggered_Nf4_lex.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "test.h"

#include "lib/Action/Fermion/action_F_Staggered_lex.h"
#include "lib/Action/Fermion/action_F_Rational.h"

#include "lib/Force/Fermion/force_F_Staggered.h"
#include "lib/Force/Fermion/force_F_Smeared.h"
#include "lib/Measurements/Fermion/fprop_Standard_lex.h"

#include "lib/Fopr/fopr_Smeared.h"

#include "lib/HMC/hmc_General.h"
#include "lib/HMC/builder_Integrator.h"

#include "lib/IO/gaugeConfig.h"

#include "lib/Tools/file_utils.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Tools/randomNumbers_Mseries.h"

#include "lib/Smear/smear.h"
#include "lib/Smear/projection.h"
#include "lib/Smear/director_Smear.h"

//====================================================================
//! Test of staggered fermion.

/*!
    This test HMC update for dynamical staggered fermions, that
    works in the same way as the even-odd version.
                                          [20 Jan 2026 H.Matsufuru]
*/

namespace Test_HMC_Staggered {
  const std::string test_name = "HMC.Staggered.Nf4_lex";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HMC_Staggered_Nf4.yaml";
  }

  //- prototype declaration
  int update_Nf4_lex(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      update_Nf4_lex
      );
#endif
  }
#endif

//====================================================================
int update_Nf4_lex(void)
  {
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test       = params_all.lookup("Test_HMC_Staggered");
    const Parameters params_action_G   = params_all.lookup("Action_G");
    const Parameters params_action_F   = params_all.lookup("Action_F_Staggered");
    const Parameters params_staggered  = params_all.lookup("Fopr_Staggered");
    const Parameters params_solver     = params_all.lookup("Solver_Nf4");
    const Parameters params_smear      = params_all.lookup("Smearing");
    const Parameters params_integrator = params_all.lookup("Builder_Integrator");
    const Parameters params_hmc        = params_all.lookup("HMC_General");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        writefile        = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    int                 i_conf           = params_test.get_int("trajectory_number");
    const int           Ntraj            = params_test.get_int("trajectory_number_step");
    const int           i_save_conf      = params_test.get_int("save_config_interval");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string           str_action_G_type = params_action_G.get_string("action_type");
    const int              Nlevels           = params_integrator.get_int("number_of_levels");
    const std::vector<int> level_action      = params_integrator.get_int_vector("level_of_actions");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write  = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  writefile    = %s\n", writefile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  i_conf       = %d\n", i_conf);
    vout.general(vl, "  Ntraj        = %d\n", Ntraj);
    vout.general(vl, "  i_save_conf  = %d\n", i_save_conf);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
    err += ParameterCheck::non_negative(i_conf);
    err += ParameterCheck::non_negative(Ntraj);
    err += ParameterCheck::non_negative(i_save_conf);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n",
                   test_name.c_str());
      exit(EXIT_FAILURE);
    }

    //- NB. Stout_SU3 is used in this test.
    if (CommonParameters::Nc() != 3) {
      vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n",
                   CommonParameters::Nc());
      return EXIT_SKIP;
    }


    RandomNumberManager::initialize(str_rand_type, seed);


    // #####  object setup  #####
    Field_G U(Nvol, Ndim);

    if (str_gconf_status == "Continue") {
      GaugeConfig(str_gconf_read).read(U, readfile);
    } else if (str_gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(U);
    } else if (str_gconf_status == "Hot_start") {
      GaugeConfig("Random").read(U);
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n",
                   test_name.c_str(), str_gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    GaugeConfig gconf_write(str_gconf_write);

    // ### Gauge action ###
    unique_ptr<Action> action_G(
        Action::New(str_action_G_type, params_action_G));

    // ### Fermion action ###
    unique_ptr<Fopr>  fopr_Nf4(
        new Fopr_Staggered(params_staggered));

    unique_ptr<Force> forceF_Nf4(
        new Force_F_Staggered(params_staggered));

    //-- smear part
    std::string proj_type  = params_smear.get_string("projection_type");
    unique_ptr<Projection> proj(
        Projection::New(proj_type, params_smear));

    std::string smear_type = params_smear.get_string("smear_type");
    unique_ptr<Smear> smear(
        Smear::New(smear_type, proj.get(), params_smear));

    unique_ptr<ForceSmear> force_smear(
        ForceSmear::New(smear_type, proj.get(), params_smear));

    unique_ptr<Director_Smear> dr_smear(
        new Director_Smear(smear.get(), params_smear));

    unique_ptr<Fopr>  fopr_smear_Nf4(
        new Fopr_Smeared(fopr_Nf4.get(), dr_smear.get()));

    unique_ptr<Force> forceF_smear_Nf4(
        new Force_F_Smeared(forceF_Nf4.get(), force_smear.get(),
                            dr_smear.get()));

    //-- Fermion solver
    std::string solver_Nf4_type = params_solver.get_string("solver_type");
    unique_ptr<Solver> solver_Nf4(
        Solver::New(solver_Nf4_type, fopr_smear_Nf4.get(),
                    params_solver));

    unique_ptr<Fprop>  fprop_Nf4(
        new Fprop_Standard_lex(solver_Nf4.get()));

    //-- Fermion action
    unique_ptr<Action> action_F4(
        new Action_F_Staggered_lex(fopr_smear_Nf4.get(),
                                   forceF_smear_Nf4.get(),
                                   fprop_Nf4.get(),
                                   fprop_Nf4.get() ));

    // ### Setup of HMC ###
    ActionList actions(Nlevels);
    actions.append(level_action[0], action_F4.get());
    actions.append(level_action[1], action_G.get());

    std::vector<Director *> directors(1);
    directors[0] = static_cast<Director *>(dr_smear.get());

    unique_ptr<Builder_Integrator> builder(
        new Builder_Integrator(actions, directors, params_integrator));

    Integrator *integrator = builder->build();

    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));

    HMC_General hmc(actions, integrator, rand.get(), params_hmc);

    Timer timer(test_name);


    // ####  Execution main part  ####
    timer.start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj);

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      vout.general(vl, "\n");
      vout.general(vl, "traj = %d\n", traj);

      result = hmc.update(U);

      if ((i_conf + traj + 1) % i_save_conf == 0) {
        std::string filename = FileUtils::generate_filename("%s-%06d", writefile.c_str(), (i_conf + traj + 1));
        gconf_write.write_file(U, filename);
      }
    }

    gconf_write.write_file(U, writefile);

    timer.report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }

} // namespace Test_HMC_Staggered

//============================================================END=====
