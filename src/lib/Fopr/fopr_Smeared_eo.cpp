/*!
        @file    fopr_Smeared_eo.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$

        @version $LastChangedRevision: 2693 $
*/

//#include "fopr_Smeared_eo.h"
#include "lib/Fopr/fopr_Smeared_eo.h"
#include "lib/Fopr/afopr_Smeared_eo-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Smeared_eo::register_factory();
}
#endif

//const std::string Fopr_Smeared_eo::class_name = "Fopr_Smeared_eo";
template<>
const std::string AFopr_Smeared_eo<Field>::class_name = "Fopr_Smeared_eo";

template class AFopr_Smeared_eo<Field>;


/*
//====================================================================
void Fopr_Smeared_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Fopr_Smeared_eo::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_Smeared_eo::set_config(Field *U)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_dr_smear->set_config(U);

  const int Nsmear = m_dr_smear->get_Nsmear();
  Field     *Uptr  = m_dr_smear->getptr_smearedConfig(Nsmear);

  m_fopr_eo->set_config(Uptr);
}


//====================================================================
double Fopr_Smeared_eo::flop_count()
{
  //- Counting of floating point operations in giga unit.
  //  not implemented, yet.

  vout.general(m_vl, "Warning at %s: flop_count() has not been implemented.\n", class_name.c_str());

  const double gflop = 0.0;

  return gflop;
}


//====================================================================
*/
//============================================================END=====
