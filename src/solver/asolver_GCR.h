#ifndef ASOLVER_GCR_H
#define ASOLVER_GCR_H

#include "solver/asolver.h"
#include "solver/aprecond_MG.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// Restarted flexible GCR using Bridge++ fields and the existing MG cycle.
// The paired direction/image basis follows multigrid.cpp's native GCR.
template<class Field>
class ASolver_GCR : public ASolver<Field>
{
  using Real = typename Field::real_t;
  using Complex = typename Field::complex_t;
  using Preconditioner = APrecond_MG<Field, AField<float, Field::IMPL>>;

 public:
  ASolver_GCR(AFopr<Field>* op, Preconditioner* preconditioner)
    : m_op(op), m_preconditioner(preconditioner) {}

  void set_parameters(const Parameters& parameters) override
  {
    const long long budget = static_cast<long long>(parameters.get_int("maximum_number_of_iteration"))
                           * parameters.get_int("maximum_number_of_restart");
    const int restart = parameters.is_set("gcr_restart") ? parameters.get_int("gcr_restart") : 16;
    if (budget < 1 || budget > (std::numeric_limits<int>::max()-2)/2 || restart < 1)
      throw std::runtime_error("invalid GCR iteration budget or restart length");
    m_budget = static_cast<int>(budget);
    const int window = std::min(restart, m_budget);
    m_recursive = parameters.get_double("convergence_criterion_squared");
    m_verified = parameters.is_set("gcr_verified_criterion_squared")
      ? parameters.get_double("gcr_verified_criterion_squared") : m_recursive;
    m_directions.resize(window);
    m_images.resize(window);
    const int nin = m_op->field_nin(), nvol = m_op->field_nvol(), nex = m_op->field_nex();
    m_residual.reset(nin, nvol, nex);
    for (int i = 0; i < window; ++i) {
      m_directions[i].reset(nin, nvol, nex);
      m_images[i].reset(nin, nvol, nex);
    }
  }

  int iterations() const { return m_iterations; }

  double flop_count() override
  {
    return m_preconditioner->flop_count() + m_applications * m_op->flop_count()
           + m_vector_flops;
  }

  void solve(Field& x, const Field& b, int& nconv, Real& diff) override
  {
#pragma omp master
    { m_preconditioner->reset_flop_count(); }
#pragma omp barrier
    const double b2 = norm2(b);
    int iterations = 0, applications = 0, active = 0;
    double vector_cost = 0, residual2 = 0;
    bool refreshed = true;
    if (b2 > 0) {
      copy(x, b);
      residual2 = refresh(x, b, applications, vector_cost);
      while (iterations < m_budget && std::isfinite(residual2) &&
             (refreshed ? residual2 > b2*m_verified : residual2 > 0)) {
        Field& z = m_directions[active];
        Field& q = m_images[active];
        m_preconditioner->mult(z, m_residual);
        m_op->mult(q, z);
        ++applications;

        // Modified Gram-Schmidt applies the same transform to z and D*z.
        for (int j = 0; j < active; ++j) {
          const Complex beta = dotc(m_images[j], q);
          axpy(q, -beta, m_images[j]);
          axpy(z, -beta, m_directions[j]);
        }
        vector_cost += 24.0 * active;
        Real image2, rhs2;
        Complex projection;
        dotc_and_norm2(projection, image2, rhs2, q, m_residual);
        vector_cost += 16.0;
        if (!(image2 > std::numeric_limits<Real>::min()) || !std::isfinite(image2))
          break;

        const Real norm = std::sqrt(image2);
        const Complex alpha = projection / norm;
        q.scal(Real(1)/norm);
        z.scal(Real(1)/norm);
        axpy(x, alpha, z);
        axpy(m_residual, -alpha, q);
        residual2 = norm2(m_residual);
        vector_cost += 24.0;
        ++iterations;
        ++active;
        refreshed = false;
        if (residual2 <= b2*m_recursive || active == static_cast<int>(m_images.size())) {
          residual2 = refresh(x, b, applications, vector_cost);
          refreshed = true;
          active = 0;
        }
      }
      if (!refreshed)
        residual2 = refresh(x, b, applications, vector_cost);
    } else {
      x.set(Real(0));
    }

#pragma omp master
    {
      m_iterations = iterations;
      nconv = m_applications = applications;
      diff = b2 == 0 ? Real(0) : static_cast<Real>(std::sqrt(residual2/b2));
      m_vector_flops = vector_cost * (b.nin()/2) * b.nex() * CommonParameters::Lvol();
      vout.general(this->m_vl, "GCR: %d iterations, %d fine operator applications, relative residual %.8e\n",
                   iterations, applications, static_cast<double>(diff));
    }
#pragma omp barrier
  }

 private:
  double refresh(const Field& x, const Field& b, int& applications, double& vector_cost)
  {
    m_op->mult(m_residual, x);
    ++applications;
    vector_cost += 6.0;
    return m_residual.residual_norm2(b);
  }

  AFopr<Field>* m_op;
  Preconditioner* m_preconditioner;
  std::vector<Field> m_directions, m_images;
  Field m_residual;
  int m_budget = 0, m_iterations = 0, m_applications = 0;
  double m_recursive = 0, m_verified = 0, m_vector_flops = 0;
};

#endif
