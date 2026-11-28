/*!
        @file    fopr_CloverTerm_eo.h

        @brief

        @author  UEDA, Satoru (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_CLOVERTERM_EO_IMPL_IMP_INCLUDED
#define FOPR_CLOVERTERM_EO_IMPL_IMP_INCLUDED

#include "operator/fopr.h"

#include "field/index_eo.h"
#include "field/field_F.h"
#include "math/gammaMatrix.h"

#include "io/bridgeIO.h"
using Bridge::vout;

class Staple_eo;
class ShiftField_eo;
class Solver_CG;
class Fopr_Clover_eo;

//! Clover term operator.

/*!
    This class is an even-odd version of clover term operator.
    The field strength is calculate when the function
    set_config() is called, and the inverse of site-diagonal
    term is determined.
    The `mode' for setting fermion operator mode is now only
    defined to the case 'D'.
    Original code was implemented by Satoru Ueda in June 2012.

    Modify this code to work.      [03 Mar 2013 Y.Namekawa]
    Multi-threading was applied to D() and mult_csw_inv().
    Previous version explicitly implements the Dirac representaion
    of gamma-matrices and thus failed for other representation.
    In the cpp files in Imp/ and other improved performance versions,
    now chiral representation is available in addition to the Dirac.
    For these implementation, performance tuning was also applied.
                                     [31 Jul 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                   [21 Mar 2015 Y.Namekawa]

    Note: mult with mode 'even' or 'odd' multiplies

      1 - csw kappa sigma_{mu nu} F_{mu nu}.

    (this is different from that of fopr_CloverTerm with mode 'D',
     which multiplies csw kappa sigma_{mu nu} F_{mu nu}. )
                                        [22 Jan 2019 I.Kanamori]
    Multi-threading is applied to set_config.
                                       [02 Dec 2021 H.Matsufuru]
 */

namespace Imp {
  class Fopr_CloverTerm_eo : public Fopr
  {
    // This class returns D_ee = 1-f_ee or D_oo = 1-f_oo

   public:
    static const std::string class_name;

   private:
    // input parameters
    double m_kappa;
    double m_cSW;
    std::vector<int> m_boundary;
    std::string m_repr;
    Bridge::VerboseLevel m_vl;

    // internal data members
    std::string m_mode;

    int m_Nvol, m_Nvol2, m_Ndim;
    int m_Nc, m_Nd, m_Ndm2;
    int m_NinF;

    //! Gamma Matrix and Sigma_{mu,nu} = -i [Gamma_mu, Gamma_nu] /2
    std::vector<GammaMatrix> m_GM, m_SG;

    Index_eo m_idx;
    // Only the enclosing Schur operator can request inverse-only storage.
    friend class ::Fopr_Clover_eo;
    bool m_inverse_only = false;
    Fopr_CloverTerm_eo(const Parameters& params, bool inverse_only)
      : m_inverse_only(inverse_only) { init(params); }

    Staple_eo *m_staple = nullptr;
    ShiftField_eo *m_shift_eo = nullptr;
    Solver_CG *m_solver = nullptr;

    Field_G m_Ueo{0, 0};                //!< even-odd gauge configuration
    Field_G m_T{0, 0};                  //!< m_T = 1 - kappa c_SW sigma F / 2

    Field_F m_Fee_inv{0, 0}, m_Foo_inv{0, 0}; //!<inverse of site-diagonal part
    Field_F m_w1{0, 0}, m_w2{0, 0};           //!< working vectors
    Field_G m_Ft{0, 0};                 //!< working vectors
    Field_G m_Cup{0, 0}, m_Cdn{0, 0}, m_Umu{0, 0};  //!< working vectors
    Field_G m_ut1{0, 0}, m_ut2{0, 0}, m_ut3{0, 0};  //!< working vectors
    // try to reduce number of working vector in set_fieldstrength !

   public:

    DEPRECATED
    Fopr_CloverTerm_eo(std::string repr) { init(repr); }

    //! standard constructor.
    Fopr_CloverTerm_eo(const Parameters& params) { init(params); }

    ~Fopr_CloverTerm_eo() { tidyup(); }

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa, const double cSW,
                        const std::vector<int> bc);

    void get_parameters(Parameters& params) const;

    void set_config(Field *U);

    void set_mode(const std::string mode);

    std::string get_mode() const { return m_mode; }

    //! return D = D^dag = 1-f_ee or 1-f_oo
    void mult(Field& v, const Field& f);

    void mult_dag(Field& v, const Field& f);

    void mult_isigma(Field_F&, const Field_F&,
                     const int mu, const int nu);

    //! multiplies 1-csw kappa sigma_{mu nu} F_{mu nu}
    void D(Field& v, const Field& f, const int ieo);

    //! explicit implementation for Dirac representation (for Imp-version).
    void D_dirac(Field& v, const Field& f, const int ieo);

    //! explicit implementation for Chiral representation (for Imp-version).
    void D_chiral(Field& v, const Field& f, const int ieo);

    //! multiplies [ 1-csw kappa sigma_{mu nu} F_{mu nu} ]^{-1}
    void mult_csw_inv(Field&, const Field&, const int ieo);

    void mult_csw_inv_dirac(Field&, const Field&, const int ieo);

    void mult_csw_inv_chiral(Field&, const Field&, const int ieo);

    //   std::vector<double> csmatrix(const int&);

    int field_nvol() { return m_Nvol2; }

    int field_nin() { return 2 * m_Nc * m_Nd; }

    int field_nex() { return 1; }

    //! returns number of floating point operations.
    double flop_count();

   private:
    void init(const std::string repr);

    void init(const Parameters& params);

    void tidyup();

    void setup();
    void allocate_fieldstrength_workspace();
    void release_fieldstrength_workspace();
    void allocate_inverse_workspace();

    void setup_gamma_matrices();

    void set_config_omp(Field *U);

    void set_config_impl(Field *U);

    void solve_csw_inv();

    void set_csw();

    //! explicit implementation for Dirac representation (for Imp-version).
    void set_csw_dirac();

    //! explicit implementation for Chiral representation (for Imp-version).
    void set_csw_chiral();

    void mult_csw(Field_F&, const Field_F&, const int ieo);

    void set_fieldstrength(Field_G&, const int, const int);

    int sg_index(const int mu, const int nu) { return mu * m_Ndim + nu; }
  };
}

using Fopr_CloverTerm_eo = Imp::Fopr_CloverTerm_eo;

#endif /* FOPR_CLOVERTERM_EO_IMPL_IMP_INCLUDED */
