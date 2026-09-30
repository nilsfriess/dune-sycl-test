#pragma once

#include <cassert>
#include <cmath>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <memory>
#include <vector>

template <class Scalar = double>
struct PoissonProblem {
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;
  using Vector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;

  /** @brief Assembles the Q1 stiffness matrix and load vector of
   *
   *    -div(a(x) grad u) = f(x) in Omega,  u = 0 on the Dirichlet part of the boundary.
   *
   *  @param gv            grid view to assemble on
   *  @param is_dirichlet  predicate on a global coordinate; true for the vertices carrying a
   *                       homogeneous Dirichlet condition. It sees global coordinates, so ranks
   *                       sharing a vertex classify it identically.
   *  @param a             scalar diffusion coefficient, evaluated at global coordinates
   *  @param f             source term, evaluated at global coordinates
   */
  template <class GridView, class IsDirichlet, class Coefficient, class Source>
  PoissonProblem(const GridView& gv, IsDirichlet is_dirichlet, Coefficient a, Source f)
  {
    using DF = typename GridView::ctype;
    constexpr int dim = GridView::dimension;

    auto& indexset = gv.indexSet();
    const int n = indexset.size(dim);

    A = std::make_shared<Matrix>();
    A->setBuildMode(Matrix::BuildMode::implicit);
    A->setImplicitBuildModeParameters(std::pow(3, dim), 0.05);
    A->setSize(n, n);

    // Create sparsity pattern, dropping couplings to vertices outside the patch
    for (const auto& e : elements(gv)) {
      auto ndofs = e.subEntities(dim);
      for (unsigned int i = 0; i < ndofs; ++i)
        for (unsigned int j = 0; j < ndofs; ++j) A->entry(indexset.subIndex(e, i, dim), indexset.subIndex(e, j, dim)) = 0;
    }
    A->compress();

    b.resize(n);
    b = 0.;

    dirichlet.assign(n, false);
    for (const auto& v : vertices(gv)) {
      const auto lidx = indexset.index(v);
      dirichlet[lidx] = is_dirichlet(v.geometry().corner(0));
    }

    // Assemble the matrix entries
    Dune::LagrangeLocalFiniteElementCache<DF, Scalar, dim, 1> fecache;

    std::vector<Dune::FieldVector<Scalar, 1>> phi;                      // shape function values
    std::vector<Dune::FieldMatrix<Scalar, 1, dim>> reference_gradients; // gradients on the reference element
    std::vector<Dune::FieldVector<Scalar, dim>> gradients;              // ... pushed forward to the element
    std::vector<std::size_t> indices;                                   // global index of each local dof

    for (const auto& e : elements(gv)) {
      const auto& fe = fecache.get(e.type());
      const auto& localbasis = fe.localBasis();
      const auto geo = e.geometry();
      const std::size_t ndofs = localbasis.size();

      phi.resize(ndofs);
      reference_gradients.resize(ndofs);
      gradients.resize(ndofs);

      indices.resize(ndofs);
      for (std::size_t i = 0; i < ndofs; ++i) {
        const auto& key = fe.localCoefficients().localKey(i);
        assert(key.codim() == dim);
        indices[i] = indexset.subIndex(e, key.subEntity(), dim);
      }

      const auto& rule = Dune::QuadratureRules<DF, dim>::rule(e.type(), 2 * localbasis.order());
      for (const auto& qp : rule) {
        const auto& pos = qp.position();
        const auto jit = geo.jacobianInverseTransposed(pos);
        const double weight = qp.weight() * geo.integrationElement(pos);

        localbasis.evaluateFunction(pos, phi);
        localbasis.evaluateJacobian(pos, reference_gradients);
        for (std::size_t i = 0; i < ndofs; ++i) jit.mv(reference_gradients[i][0], gradients[i]);

        const auto global = geo.global(pos);
        const Scalar a_x = a(global);
        const Scalar f_x = f(global);

        for (std::size_t i = 0; i < ndofs; ++i) {
          b[indices[i]][0] += f_x * phi[i][0] * weight;
          for (std::size_t j = 0; j < ndofs; ++j) (*A)[indices[i]][indices[j]][0][0] += a_x * (gradients[i] * gradients[j]) * weight;
        }
      }
    }

    // Homogeneous Dirichlet conditions: replace each constrained row by the identity row and zero
    // its load entry. The columns are left alone, so A is not symmetric.
    for (auto ri = A->begin(); ri != A->end(); ++ri) {
      if (not dirichlet[ri.index()]) continue;
      for (auto ci = ri->begin(); ci != ri->end(); ++ci) *ci = (ci.index() == ri.index()) ? 1. : 0.;
      b[ri.index()] = 0.;
    }
  }

  std::shared_ptr<Matrix> A;
  Vector b;

  /// Vertices carrying a homogeneous Dirichlet condition of the global problem.
  std::vector<bool> dirichlet;
};
