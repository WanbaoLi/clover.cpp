#pragma once
#include "Force/Fermion/force_F_Clover_Nf2.h"
#include <stdexcept>

// One HMC force request supplies an immutable smeared gauge field. Bridge's
// smearing wrapper revisits set_config before each bilinear; the Clover term
// and its staples can be prepared once for this field and this mass.
class FixedGaugeCloverForce final : public Force_F_Clover_Nf2 {
  Field* const links_;
public:
  FixedGaugeCloverForce(Field_G* links, double kappa, double csw,
                        const std::vector<int>& boundary)
    : Force_F_Clover_Nf2("Dirac"), links_(links) {
    set_parameters(kappa,csw,boundary);
    Force_F_Clover_Nf2::set_config(links);
  }
  void set_config(Field* links) override {
    if(links!=links_) throw std::logic_error("HMC force used with a different gauge field");
  }
};
