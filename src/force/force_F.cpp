/*!
        @file    force_F.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$

        @version $LastChangedRevision: 2693 $
*/

#include "force/force_F.h"

#include "field/field_thread-inc.h"

template<>
const std::string Force::class_name = "Force_F";

//====================================================================
template<>
void Force::init()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_force_tmp = new Field_G(Nvol, Ndim);
  //m_force1    = new Field_G(Nvol, 1);
}

//====================================================================
template<>
void Force::tidyup()
{
  delete m_force_tmp;
  // delete m_force1;
}

//====================================================================
template<>
void Force::mult_generator(Field& force)
{
  vout.crucial("%s: mult_generator(Field&) must not be called\n",
	       class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<>
void Force::mult_generator(Field_G& force)
{
#pragma omp barrier

  const int Nvol = force.nvol();
  const int Ndim = force.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int site = is; site < ns; ++site) {
      Mat_SU_N u = m_U->mat(site, mu);
      u *= force.mat(site, mu);
      u.at();
      u *= -2.0;
      force.set_mat(site, mu, u);
    }
  }
#pragma omp barrier

  /*
  for (int mu = 0; mu < Ndim; ++mu) {
    mult_Field_Gnn(*m_force1, 0, *m_U, mu, force, mu);
    at_Field_G(*m_force1,0);
    scal(*m_force1, -2.0);
    copy(force, mu, *m_force1, 0);
#pragma omp barrier
  }
  */
}

//====================================================================
// default templates for core and core1
template<>
void Force::force_core(Field& force_, const Field& eta)
{
  force_udiv(*m_force_tmp, eta);

  mult_generator(*m_force_tmp);

  copy(force_, *m_force_tmp);
#pragma omp barrier

}


//====================================================================
template<>
void Force::force_core1(Field& force_, const Field& zeta, const Field& eta)
{
  force_udiv1(*m_force_tmp, zeta, eta);

  mult_generator(*m_force_tmp);

  copy(force_, *m_force_tmp);
#pragma omp barrier

}


//====================================================================

// explicit instanciation.
template class AForce_F<Field>;

//============================================================END=====
