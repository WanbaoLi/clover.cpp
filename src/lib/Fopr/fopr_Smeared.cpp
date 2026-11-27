/*!
        @file    fopr_Smeared.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Fopr/fopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Smeared::register_factory();
}
#endif

//const std::string Fopr_Smeared::class_name = "Fopr_Smeared";
template<>
const std::string AFopr_Smeared<Field>::class_name = "Fopr_Smeared";

template class AFopr_Smeared<Field>;

/*
//====================================================================
void Fopr_Smeared::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Fopr_Smeared::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Fopr_Smeared::set_config(Field *U)
{
  m_dr_smear->set_config(U);

  const int Nsmear = m_dr_smear->get_Nsmear();
  Field     *Uptr  = m_dr_smear->getptr_smearedConfig(Nsmear);

  m_fopr->set_config(Uptr);
}


//====================================================================
double Fopr_Smeared::flop_count()
{
  //- Counting of floating point operations in giga unit.
  //  not implemented, yet.

  vout.general(m_vl, "Warning at %s: flop_count() does not include smearing opertions.\n",
               class_name.c_str());

  const double gflop = m_fopr->flop_count();

  return gflop;
}


//====================================================================
*/

//============================================================END=====
