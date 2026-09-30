# dune-sycl-test

Experiment with using SYCL to accelerate Dune solvers.

## Build

```sh
git clone --recursive <repo-url>
cd dune-sycl-test
mkdir build && cd build
cmake .. -DAdaptiveCpp_DIR=/path/to/AdaptiveCpp/lib/cmake/AdaptiveCpp
make
```
