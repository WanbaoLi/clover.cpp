/*!
        @file    director_Smear.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "gauge/director_Smear.h"
#include "runtime/threadManager.h"

const std::string Director_Smear::class_name = "Director_Smear";

//====================================================================
void Director_Smear::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_Nsmear       = 0;
  m_U            = 0;
  m_status_linkv = 0;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Director_Smear::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_Nsmear       = 0;
  m_U            = 0;
  m_status_linkv = 0;

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Director_Smear::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int Nsmear;

  int err = 0;
  err += params.fetch_int("number_of_smearing", Nsmear);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Nsmear);
}


//====================================================================
void Director_Smear::set_parameters(const int Nsmear)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_Nsmear = Nsmear;

    m_Usmear.resize(m_Nsmear);
    if (m_Nsmear > 0) {
      for (int ismear = 0; ismear < m_Nsmear; ++ismear) {
        m_Usmear[ismear].reset(Nvol, Ndim);
        vout.detailed(m_vl, " size of Usmear[ismear] was set.\n");
      }
    }
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  Nsmear = %d\n", Nsmear);

}


//====================================================================
void Director_Smear::get_parameters(Parameters& params) const
{
  params.set_int("number_of_smearing", m_Nsmear);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Director_Smear::set_config(Field *U)
{
  m_U = (Field_G *)U;

  if (m_status_linkv == 0) smear();
}


//====================================================================
Field *Director_Smear::getptr_smearedConfig(const int ismear)
{
  assert(m_U != 0);

  if (ismear == 0) {
    return m_U;
  } else {
    return &m_Usmear[ismear-1];
  }
}


//====================================================================
Field_G *Director_Smear::get_config()
{
  return m_Nsmear == 0 ? m_U : &m_Usmear[m_Nsmear - 1];
}


//====================================================================
Field_G *Director_Smear::get_config(const int ismear)
{
  if (ismear == 0) return m_U;

  return &m_Usmear[ismear-1];

}

//====================================================================
void Director_Smear::smear()
{
  if (m_Nsmear > 0) {
    for (int ismear = 0; ismear < m_Nsmear; ++ismear) {
      if (ismear == 0) {
        m_smear->smear(m_Usmear[ismear], *m_U);
      } else {
        m_smear->smear(m_Usmear[ismear], m_Usmear[ismear-1]);
      }

    }
  }

  ++m_status_linkv;
}

//============================================================END=====
