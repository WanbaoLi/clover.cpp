/*!
        @file    test_HMC_Domainwall_Nf2_PV.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2013-04-27 12:28:50 #$
        @version $LastChangedRevision: 2693 $
*/

#include "test.h"

#include "lib/Action/Fermion/action_F_Standard_lex.h"
#include "lib/Action/Fermion/action_F_Ratio_lex.h"

#include "lib/Fopr/fopr_Smeared.h"

#include "lib/Force/Fermion/force_F_Domainwall.h"
#include "lib/Force/Fermion/force_F_Smeared.h"

#include "lib/HMC/hmc_General.h"
#include "lib/HMC/builder_Integrator.h"

#include "lib/IO/gaugeConfig.h"

#include "lib/Measurements/Fermion/fprop_Standard_lex.h"

#include "lib/Tools/file_utils.h"
#include "lib/Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of HMC update for domain-wall fermion operator.

/*!
    This class tests HMC update for dynamicalMoebius domain-wall
    fermions.
                                          [01 Feb 2026 H.Matsufuru]
 */

namespace Test_HMC_Domainwall {
  const std::string test_name = "HMC.Domainwall.Moebius_Nf2_PV";

  //- test-private parameters
  namespace {
    const std::string filename_input
                       = "test_HMC_Domainwall_Moebius_Nf2_PV.yaml";
  }

  //- prototype declaration
  int update_Moebius_Nf2_PV(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
  static const bool is_registered = TestManager::RegisterTest(
      test_name, update_Moebius_Nf2_PV);
#endif
  }
#endif

//====================================================================
int update_Moebius_Nf2_PV(void)
{
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test       = params_all.lookup("Test_HMC_Domainwall");
    const Parameters params_action_G   = params_all.lookup("Action_G");
    const Parameters params_proj       = params_all.lookup("Projection");
    const Parameters params_smear      = params_all.lookup("Smear");
    const Parameters params_dr_smear   = params_all.lookup("Director_Smear");
    const Parameters params_dw         = params_all.lookup("Fopr_Domainwall");
    const Parameters params_dw_pv      = params_all.lookup("Fopr_Domainwall_PauliVillars");
    const Parameters params_solver_MD  = params_all.lookup("Solver_MD");
    const Parameters params_solver_H   = params_all.lookup("Solver_H");
    const Parameters params_integrator = params_all.lookup("Builder_Integrator");
    const Parameters params_hmc        = params_all.lookup("HMC_General");

    const string gconf_status = params_test.get_string("gauge_config_status");
    const string gconf_read_type
                              = params_test.get_string("gauge_config_type_input");
    const string readfile     = params_test.get_string("config_filename_input");
    const string gconf_write_type
                              = params_test.get_string("gauge_config_type_output");
    const string writefile    = params_test.get_string("config_filename_output");
    const string rand_type    = params_test.get_string("random_number_type");
    int          i_conf       = params_test.get_int("trajectory_number");
    const int    Ntraj        = params_test.get_int("trajectory_number_step");
    const int    i_save_conf  = params_test.get_int("save_config_interval");
    const string vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string action_G_type  = params_action_G.get_string("action_type") ;
    const string gmset_type     = params_dw.get_string("gamma_matrix_type");
    const string proj_type      = params_proj.get_string("projection_type");
    const string smear_type     = params_smear.get_string("smear_type");
    const string solver_MD_type = params_solver_MD.get_string("solver_type");
    const string solver_H_type  = params_solver_H.get_string("solver_type");
    const int    Nlevels        = params_integrator.get_int("number_of_levels");

    const std::vector<int> level_action
                    = params_integrator.get_int_vector("level_of_actions");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status     = %s\n", gconf_status.c_str());
    vout.general(vl, "  gconf_read_type  = %s\n", gconf_read_type.c_str());
    vout.general(vl, "  readfile         = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write_type = %s\n", gconf_write_type.c_str());
    vout.general(vl, "  writefile        = %s\n", writefile.c_str());
    vout.general(vl, "  rand_type        = %s\n", rand_type.c_str());
    vout.general(vl, "  i_conf           = %d\n", i_conf);
    vout.general(vl, "  Ntraj            = %d\n", Ntraj);
    vout.general(vl, "  i_save_conf      = %d\n", i_save_conf);
    vout.general(vl, "  vlevel           = %s\n", vlevel.c_str());
    vout.general(vl, "  gmset_type       = %s\n", gmset_type.c_str());
    vout.general(vl, "  proj_type        = %s\n", proj_type.c_str());
    vout.general(vl, "  smear_type       = %s\n", smear_type.c_str());
    vout.general(vl, "  solver_MD_type   = %s\n", solver_MD_type.c_str());
    vout.general(vl, "  solver_H_type    = %s\n", solver_H_type.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(gconf_status);
    err += ParameterCheck::non_negative(i_conf);
    err += ParameterCheck::non_negative(Ntraj);
    err += ParameterCheck::non_negative(i_save_conf);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n",
                   test_name.c_str());
      exit(EXIT_FAILURE);
    }

    if (proj_type == "Stout_SU3") {
      if (CommonParameters::Nc() != 3) {
        vout.crucial(vl, "check skipped: Nc = 3 is needed, but Nc = %d.\n\n",
                     CommonParameters::Nc());
        return EXIT_SKIP;
      }
    }

    // #####  object setup  #####
    Field_G U(Nvol, Ndim);

    if (gconf_status == "Continue") {
      GaugeConfig(gconf_read_type).read(U, readfile);
    } else if (gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(U);
    } else if (gconf_status == "Hot_start") {
      GaugeConfig("Random").read(U);
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n",
                   test_name.c_str(), gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    GaugeConfig gconf_write(gconf_write_type);

    //####  gauge part  ####
    unique_ptr<Action> action_G(
        Action::New(action_G_type, params_action_G));

    //####  smearing  ####
    unique_ptr<Projection> proj(
        Projection::New(proj_type, params_proj));
    unique_ptr<Smear>      smear(
        Smear::New(smear_type, proj.get(), params_smear));

    unique_ptr<ForceSmear> force_smear(
        ForceSmear::New(smear_type, proj.get(), params_smear));

    unique_ptr<Director_Smear> dr_smear(
        new Director_Smear(smear.get(), params_dr_smear));

    //####  dynamical part  ####
    unique_ptr<Fopr> fopr_dw(
        Fopr::New("Domainwall", params_dw));

    unique_ptr<Force> forceF_dw(
        new Force_F_Domainwall(params_dw));

    unique_ptr<Fopr> fopr_smear(
        Fopr::New("Smeared", fopr_dw.get(), dr_smear.get()));

    unique_ptr<Force> forceF_smear(
        new Force_F_Smeared(forceF_dw.get(), force_smear.get(),
                            dr_smear.get()));

    unique_ptr<Solver> solver_MD(
        Solver::New(solver_MD_type, fopr_smear.get(), params_solver_MD));

    unique_ptr<Fprop> fprop_MD(
        new Fprop_Standard_lex(solver_MD.get()));

    unique_ptr<Solver> solver_H(
        Solver::New(solver_H_type, fopr_smear.get(), params_solver_H));

    unique_ptr<Fprop> fprop_H(
        new Fprop_Standard_lex(solver_H.get()));

    //####  Pauli-Villars part  ####
    unique_ptr<Fopr> fopr_dw_pv(
        Fopr::New("Domainwall", params_dw_pv));

    unique_ptr<Force> forceF_dw_pv(
        new Force_F_Domainwall(params_dw_pv));

    unique_ptr<Fopr>  fopr_smear_pv(
        Fopr::New("Smeared", fopr_dw_pv.get(), dr_smear.get()));

    unique_ptr<Force> forceF_smear_pv(
        new Force_F_Smeared(forceF_dw_pv.get(), force_smear.get(),
                            dr_smear.get()));

    unique_ptr<Solver> solver_H_pv(
        Solver::New(solver_H_type, fopr_smear_pv.get(), params_solver_H));

    unique_ptr<Fprop> fprop_H_pv(
        new Fprop_Standard_lex(solver_H_pv.get()));

    //####  settting ratio action  ####
    unique_ptr<Action> action_F(
        new Action_F_Ratio_lex(fopr_smear_pv.get(), forceF_smear_pv.get(),
                               fopr_smear.get(), forceF_smear.get(),
                               fprop_H_pv.get(),
                               fprop_MD.get(), fprop_H.get()));

    //####  settting actions  ####
    ActionList actions(Nlevels);
    actions.append(level_action[0], action_F.get());
    actions.append(level_action[1], action_G.get());

    std::vector<Director *> directors(1);
    directors[0] = static_cast<Director *>(dr_smear.get());

    unique_ptr<Builder_Integrator> builder(
        new Builder_Integrator(actions, directors, params_integrator));

    Integrator *integrator = builder->build();

    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));

    HMC_General hmc(actions, directors, integrator, rand.get(), params_hmc);

    Timer timer(test_name);

    // ####  Execution main part  ####
    timer.start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj);

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      vout.general(vl, "\n");
      vout.general(vl, "---------------------------------------------------\n");
      vout.general(vl, "traj = %d\n", traj);

      result = hmc.update(U);

      if ((i_conf + traj + 1) % i_save_conf == 0) {
        std::string filename = FileUtils::generate_filename("%s-%06d", writefile.c_str(), (i_conf + traj + 1));
        gconf_write.write_file(U, filename);
      }
    }

    gconf_write.write_file(U, writefile);

    timer.report();

    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }

}

} // namespace Test_HMC_Domainwall
