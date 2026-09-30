#include <algorithm>
#include <cmath>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solver.hh>
#include <dune/istl/solvers.hh>
#include <iostream>
#include <sycl/sycl.hpp>

#include "mat.hh"
#include "poisson_problem.hh"
#include "vec.hh"

int main()
{
  sycl::queue q{sycl::property::queue::in_order{}};

  const auto device = q.get_device();
  std::cout << "SYCL device: " << device.get_info<sycl::info::device::name>() << ", vendor: " << device.get_info<sycl::info::device::vendor>() << "\n";

  const int dim = 2;
  const int gridsize = 512;

  using Grid = Dune::YaspGrid<dim>;
  Grid grid({1., 1.}, {gridsize, gridsize});
  auto gv = grid.leafGridView();

  // Homogeneous Dirichlet conditions on the whole boundary of the unit square.
  auto is_dirichlet = [](const auto& x) {
    const double eps = 1e-10;
    return std::ranges::any_of(x, [=](double xi) { return xi < eps or xi > 1. - eps; });
  };
  auto coefficient = [](const auto&) { return 1.; };
  auto source = [](const auto&) { return 1.; };

  PoissonProblem p(gv, is_dirichlet, coefficient, source);
  auto A = SyclMat<double>::from_bcrs(q, *p.A);
  auto b = SyclVec<double>::from_host_vector(q, p.b);
  auto x = b;
  x = 0;

  Dune::MatrixAdapter<SyclMat<double>, SyclVec<double>, SyclVec<double>> op(A);
  Dune::Richardson<SyclVec<double>, SyclVec<double>> prec;

  double reduction = 1e-8;
  int maxit = 1000;
  int verbose = 1;
  Dune::CGSolver<SyclVec<double>> solver(op, prec, reduction, maxit, verbose);

  Dune::InverseOperatorResult res;
  solver.apply(x, b, res);

  return res.converged ? 0 : 1;
}
