/*!
        @file    force_F_Smeared.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Force/Fermion/force_F_Smeared.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Force_F_Smeared::class_name = "Force_F_Smeared";

//====================================================================
void Force_F_Smeared::init(const Parameters& params)
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

  set_parameters(params);

  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_force2 = new Field_G(Nvol, Ndim);
  m_force1 = new Field_G(Nvol, Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Force_F_Smeared::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_force2 = new Field_G(Nvol, Ndim);
  m_force1 = new Field_G(Nvol, Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Force_F_Smeared::tidyup()
{
  delete m_force2;
  delete m_force1;
}


//====================================================================
void Force_F_Smeared::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Force_F_Smeared::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Smeared::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_director_smear->set_config(U);
  m_force->set_config(m_director_smear->get_config());
}

//====================================================================
void Force_F_Smeared::force_udiv(Field& force_, const Field& eta)
{
#pragma omp barrier

  const int Nsmear = m_director_smear->get_Nsmear();

  if (Nsmear == 0) {
    m_force->force_udiv(*m_force2, eta);
  } else {
    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);
    m_force->force_udiv(*m_force2, eta);

    mult_jacobian(*m_force2);
  }

  copy(force_, *m_force2);
#pragma omp barrier
}


//====================================================================
void Force_F_Smeared::force_udiv1(Field& force_,
                                  const Field& zeta, const Field& eta)
{
#pragma omp barrier

  const int Nsmear = m_director_smear->get_Nsmear();

  if (Nsmear == 0) {
    m_force->force_udiv1(*m_force2, zeta, eta);
  } else {
    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);
    m_force->force_udiv1(*m_force2, zeta, eta);

    mult_jacobian(*m_force2);
  }

  copy(force_, *m_force2);
#pragma omp barrier

}


//====================================================================
void Force_F_Smeared::mult_jacobian(Field_G& force)
{
#pragma omp barrier

  const int Nsmear = m_director_smear->get_Nsmear();

  copy(*m_force1, force);
#pragma omp barrier

  for (int ismear = Nsmear - 1; ismear >= 0; --ismear) {
    Field_G *Uptr = m_director_smear->get_config(ismear);

    m_force_smear->force_udiv(force, *m_force1, *Uptr);

    if (ismear > 0) copy(*m_force1, force);  // ftmp = force;
#pragma omp barrier
  }
}


//============================================================END=====
