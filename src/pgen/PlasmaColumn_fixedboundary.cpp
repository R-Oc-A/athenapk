//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
/*
   This is a problem generator for the application of a fixed velocity profile 
   applied on othe lower boundary
   of a fixed loop. 
   This Test relies on the plasma column advection test, it differs in that
   1) There'll be no advection
   2) The boundary conditions are not from the default vals, of course.alignas
   3) There might be a difference on the profile of the magnetic field.
*/
//========================================================================================
//! \file PlasmaColumn_fixedboundary.cpp
//! \brief Problem generator for advection of a plasma column test.
//!
//! Can only be run in 2D or 3D.  Input parameters are:
//!  -  problem/rad   = radius of plasma column
//!  -  problem/amp   = amplitude of vector potential (and therefore B)
//!  -  problem/vflow = flow velocity
//!  -  problem/drat  = density ratio in loop, to test density advection and conduction
//! The flow is automatically set to run along the diagonal.
//!
//! Various test cases are possible:
//!  - (iprob=1): plasma column in x1-x2 plane (cylinder in 3D)
//!  - (iprob=2): plasma column in x2-x3 plane (cylinder in 3D)
//!  - (iprob=3): plasma column in x3-x1 plane (cylinder in 3D)
//!  - (iprob=4): rotated cylindrical plasma column in 3D.
//!  - (iprob=5): spherical plasma column in rotated plane
//!
//! REFERENCE: T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
//! constrined transport", JCP, 205, 509 (2005)
//========================================================================================

// C++ headers
#include <algorithm> // min, max
#include <cmath>     // sqrt()
#include <cstdio>    // fopen(), fprintf(), freopen()
#include <iostream>  // endl
#include <sstream>   // stringstream
#include <stdexcept> // runtime_error
#include <string>    // c_str()

// Parthenon headers
#include "mesh/mesh.hpp"
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include "../../external/parthenon/src/bvals/boundary_conditions.hpp"

// Athena headers
#include "../main.hpp"
#include "outputs/outputs.hpp"

namespace PlasmaColumn_fixedboundary {
using namespace parthenon::package::prelude;

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief field loop advection problem generator for 2D/3D problems.
//========================================================================================

void ProblemGenerator(MeshBlock *pmb, ParameterInput *pin) {
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  Real gm1 = pin->GetReal("hydro", "gamma") - 1.0;

  // Read initial conditions, diffusion coefficients (if needed)
  Real rad = pin->GetReal("problem/PlasmaColumn_fixedboundary", "rad");
  Real amp = pin->GetReal("problem/PlasmaColumn_fixedboundary", "amp");
  Real B0_ = amp;
  Real vflow = pin->GetReal("problem/PlasmaColumn_fixedboundary", "vflow");
  Real drat = pin->GetOrAddReal("problem/PlasmaColumn_fixedboundary", "drat", 1.0);
  int iprob = pin->GetInteger("problem/PlasmaColumn_fixedboundary", "iprob");
  Real ang_2, cos_a2(0.0), sin_a2(0.0), lambda(0.0);

  Real x1size =
      pmb->pmy_mesh->mesh_size.xmax(X1DIR) - pmb->pmy_mesh->mesh_size.xmin(X1DIR);
  Real x2size =
      pmb->pmy_mesh->mesh_size.xmax(X2DIR) - pmb->pmy_mesh->mesh_size.xmin(X2DIR);
  Real x3size =
      pmb->pmy_mesh->mesh_size.xmax(X3DIR) - pmb->pmy_mesh->mesh_size.xmin(X3DIR);

  // Use vector potential to initialize field loop
  auto &coords = pmb->coords;

  //if ((SQR(coords.Xc<1>(i)) + SQR(coords.Xc<2>(j)) + SQR(coords.Xc<3>(k))) <
  //    rad * rad) {
  //}
  // Initialize density and momenta.  If drat != 1, then density and temperature will be
  // different inside loop than background values

  auto &mbd = pmb->meshblock_data.Get();
  auto &u_dev = mbd->Get("cons").data;
  // initializing on host
  auto u = u_dev.GetHostMirrorAndCopy();
  for (int k = kb.s; k <= kb.e; k++) {
    for (int j = jb.s; j <= jb.e; j++) {
      for (int i = ib.s; i <= ib.e; i++) {
        u(IDN, k, j, i) = 1.0;
        if ((SQR(coords.Xc<1>(i)) + SQR(coords.Xc<2>(j)) ) <
            rad * rad) {
          u(IDN, k, j, i) = drat;
        }
        u(IM1, k, j, i) = 0.0; //u(IDN, k, j, i) * vflow * x1size;
        u(IM2, k, j, i) = 0.0;//u(IDN, k, j, i) * vflow * x2size;
        u(IM3, k, j, i) = 0.0;//u(IDN, k, j, i) * vflow * x3size;
        u(IB1, k, j, i) = 0.0;
        u(IB2, k, j, i) = 0.0;
        u(IB3, k, j, i) = B0_;
        u(IEN, k, j, i) = 
            1.0 / gm1 +
            0.5 * (SQR(u(IB1, k, j, i)) + SQR(u(IB2, k, j, i)) + SQR(u(IB3, k, j, i))) +
            0.5 * (SQR(u(IM1, k, j, i)) + SQR(u(IM2, k, j, i)) + SQR(u(IM3, k, j, i))) /
                u(IDN, k, j, i);
      }
    }
  }
  // copy initialized vars to device
  u_dev.DeepCopy(u);
}

void CustomX3(std::shared_ptr<MeshBlockData<Real>> &mbd, bool coarse) {
  parthenon::BoundaryFunction::OutflowInnerX3(mbd, coarse);
  auto pmb = mbd->GetBlockPointer();
  auto cons = mbd->PackVariables(std::vector<std::string>{"cons"}, coarse);
  // TODO(pgrete) Add par_for_bndry to Parthenon without requiring nb
  const auto nb = IndexRange{0, 3};
  const auto rho_wind_ = 0.0;
  const auto mom_wind_ = 5.0e-2;
  const auto rhoe_wind_ = 0.0;
  const auto Bx_ = 0.0;
  const auto By_ = 0.0;
  const auto Bz_ = 0.0;
  const bool fine = false;
  pmb->par_for_bndry(
      "CustomX3", nb, IndexDomain::inner_x3, parthenon::TopologicalElement::CC,
      coarse, fine, KOKKOS_LAMBDA(const int &, const int &k, const int &j, const int &i) {
        cons(IM2, k, j, i) += cons(IDN,k,j,i)*mom_wind_;
      });
}

} // namespace PlasmaColumn_fixedboundary