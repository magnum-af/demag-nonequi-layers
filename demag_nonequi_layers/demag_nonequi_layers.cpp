#include <algorithm>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS // surpress warning
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cassert>
#include <iostream>
#include <thread>
#include <vector>

namespace magnumaf {

struct NonequiMesh {
  std::size_t nx, ny; //!< Number of cells in x, y
  double dx, dy;  //!< Distance between equidistant x, y cells in [m]
  std::vector<double>
      z_spacing; //!< Thickness of each layer along z-axis in [m]
  std::size_t nz() const { return z_spacing.size(); }
};

namespace util {
// inline int ij2k(const int i, const int j, const int n) {
inline long int ij2k(long int i, long int j, long int n) {
  // template <typename T> T ij2k(const T i, const T j, const T n) {
  return (n * (n + 1) / 2) - (n - i) * ((n - i) + 1) / 2 + j - i;
}
} // namespace util

namespace newell_nonequi {

inline std::size_t nx_expanded(std::size_t nx) { return 2 * nx; }
inline std::size_t ny_expanded(std::size_t ny) { return 2 * ny; }

double f(double x, double y, double z) {
  x = std::fabs(x);
  y = std::fabs(y);
  z = std::fabs(z);
  const double R = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
  const double xx = std::pow(x, 2);
  const double yy = std::pow(y, 2);
  const double zz = std::pow(z, 2);

  double result = 1.0 / 6.0 * (2.0 * xx - yy - zz) * R;
  if (xx + zz > 0) {
    result += y / 2.0 * (zz - xx) * asinh(y / (sqrt(xx + zz)));
  }
  if (xx + yy > 0) {
    result += z / 2.0 * (yy - xx) * asinh(z / (sqrt(xx + yy)));
  }
  if (x * R > 0) {
    result += -x * y * z * atan(y * z / (x * R));
  }
  return result;
}

double g(double x, double y, double z) {
  z = fabs(z);
  const double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  const double xx = pow(x, 2);
  const double yy = pow(y, 2);
  const double zz = pow(z, 2);

  double result = -x * y * R / 3.0;
  if (xx + yy > 0) {
    result += x * y * z * asinh(z / (sqrt(xx + yy)));
  }
  if (yy + zz > 0) {
    result += y / 6.0 * (3.0 * zz - yy) * asinh(x / (sqrt(yy + zz)));
  }
  if (xx + zz > 0) {
    result += x / 6.0 * (3.0 * zz - xx) * asinh(y / (sqrt(xx + zz)));
  }
  if (z * R > 0) {
    result += -pow(z, 3) / 6.0 * atan(x * y / (z * R));
  }
  if (y * R != 0) {
    result += -z * yy / 2.0 * atan(x * z / (y * R));
  }
  if (x * R != 0) {
    result += -z * xx / 2.0 * atan(y * z / (x * R));
  }
  return result;
}

// double F2(const double x, const double y, const double z){
//    return f(x, y, z);
//    //Last three terms cancel out//return f(x, y, z) - f(x, 0, z) - f(x, y,
//    0)
//    + f(x, 0, 0);
//}

double F1(const double x, const double y, const double z, const double dz,
          const double dZ) {
  return f(x, y, z + dZ) - f(x, y, z) - f(x, y, z - dz + dZ) + f(x, y, z - dz);
}

double F0(const double x, const double y, const double z, const double dy,
          const double dY, const double dz, const double dZ) {
  return F1(x, y + dY, z, dz, dZ) - F1(x, y, z, dz, dZ) -
         F1(x, y - dy + dY, z, dz, dZ) + F1(x, y - dy, z, dz, dZ);
}

double Nxx(const double x, const double y, const double z, const double dx,
           const double dy, const double dz, const double dX, const double dY,
           const double dZ) {
  // x, y, z is vector from source cuboid to target cuboid
  // dx, dy, dz are source cuboid dimensions
  // dX, dY, dZ are target cuboid dimensions
  // const double tau = dX * dY * dZ;// Defining dX, dY, dZ as target cuboid
  // (one could alternatively choose dx, dy, dz with implications on x, y, z)
  // return -1./(4.0 * M_PI * tau) * (
  return -1. / (4.0 * M_PI) *
         (F0(x, y, z, dy, dY, dz, dZ) - F0(x - dx, y, z, dy, dY, dz, dZ) -
          F0(x + dX, y, z, dy, dY, dz, dZ) +
          F0(x - dx + dX, y, z, dy, dY, dz, dZ));
}

// double G2(const double x, const double y, const double z){
//    return g(x, y, z);
//    //return g(x, y, z) - g(x, y, 0);
//    //return g(x, y, z) - g(x, 0, z) - g(x, y, 0) + g(x, 0, 0);
//}

double G1(const double x, const double y, const double z, const double dz,
          const double dZ) {
  return g(x, y, z + dZ) - g(x, y, z) - g(x, y, z - dz + dZ) + g(x, y, z - dz);
}

double G0(const double x, const double y, const double z, const double dy,
          const double dY, const double dz, const double dZ) {
  return G1(x, y + dY, z, dz, dZ) - G1(x, y, z, dz, dZ) -
         G1(x, y - dy + dY, z, dz, dZ) + G1(x, y - dy, z, dz, dZ);
}

double Nxy(const double x, const double y, const double z, const double dx,
           const double dy, const double dz, const double dX, const double dY,
           const double dZ) {
  // x, y, z is vector from source cuboid to target cuboid
  // dx, dy, dz are source cuboid dimensions
  // dX, dY, dZ are target cuboid dimensions
  // const double tau = dX * dY * dZ;// Defining dX, dY, dZ as target cuboid
  // (one could alternatively choose dx, dy, dz with implications on x, y, z)
  // return -1./(4.0 * M_PI * tau) * (
  return -1. / (4.0 * M_PI) *
         (G0(x, y, z, dy, dY, dz, dZ) - G0(x - dx, y, z, dy, dY, dz, dZ) -
          G0(x + dX, y, z, dy, dY, dz, dZ) +
          G0(x - dx + dX, y, z, dy, dY, dz, dZ));
}

double nonequi_index_distance(const std::vector<double> &spacings,
                              const std::size_t i, const std::size_t j,
                              const bool verbose) {
  // Calculates the signed distance beween elements by summing up i < j:
  // sum_(k=i)^(j-1)[spacings[k]] or i > j: sum_(k=j)^(i-1)[ - spacings[k]]
  // Note that spacings[spacings.size()-1] is never used
  if (verbose and (i == spacings.size() or j == spacings.size())) {
    std::cout
        << "Warning: in nonequi_index_distance: index == vector.size(), the "
           "distance includes the last element which is not wanted "
           "behaviour"
        << std::endl;
  }

  double result = 0;
  if (i > j) {
    for (std::size_t k = i; k > j; k--) {
      result -= spacings[k - 1];
    }
  } else {
    for (std::size_t k = i; k < j; k++) {
      result += spacings[k];
    }
  }
  return result;
}

void init_N(const NonequiMesh &nemesh, std::vector<double> &N,
            std::size_t ix_start, std::size_t ix_end) {
  for (std::size_t ix = ix_start; ix < ix_end; ix++) {
    const long int jx =
        (ix + nx_expanded(nemesh.nx) / 2) % nx_expanded(nemesh.nx) -
        nx_expanded(nemesh.nx) / 2;
    for (std::size_t iy = 0; iy < ny_expanded(nemesh.ny); iy++) {
      const long int jy =
          (iy + ny_expanded(nemesh.ny) / 2) % ny_expanded(nemesh.ny) -
          ny_expanded(nemesh.ny) / 2;
      for (std::size_t i_source = 0; i_source < nemesh.nz(); i_source++) {
        for (std::size_t i_target = 0; i_target < nemesh.nz(); i_target++) {

          if (i_source <= i_target) {
            const long int idx =
                6 * (util::ij2k(i_source, i_target, nemesh.nz()) +
                     ((nemesh.nz() * (nemesh.nz() + 1)) / 2) *
                         (iy + ny_expanded(nemesh.ny) * ix));
            // std::cout << "idx=" << idx << " of " <<
            // nx_expanded(nemesh.nx) * ny_expanded(nemesh.ny) *
            // (nemesh.nz() * (nemesh.nz() + 1))/2 * 6
            // << std::endl;
            const double x = nemesh.dx * static_cast<double>(jx);
            const double y = nemesh.dy * static_cast<double>(jy);
            const double z = nonequi_index_distance(nemesh.z_spacing, i_source,
                                                    i_target, true);

            N[idx + 0] = newell_nonequi::Nxx(
                x, y, z, nemesh.dx, nemesh.dy, nemesh.z_spacing[i_source],
                nemesh.dx, nemesh.dy, nemesh.z_spacing[i_target]);
            N[idx + 1] = newell_nonequi::Nxy(
                x, y, z, nemesh.dx, nemesh.dy, nemesh.z_spacing[i_source],
                nemesh.dx, nemesh.dy, nemesh.z_spacing[i_target]);
            N[idx + 2] = newell_nonequi::Nxy(
                x, z, y, nemesh.dx, nemesh.z_spacing[i_source], nemesh.dy,
                nemesh.dx, nemesh.z_spacing[i_target], nemesh.dy);
            N[idx + 3] = newell_nonequi::Nxx(
                y, z, x, nemesh.dy, nemesh.z_spacing[i_source], nemesh.dx,
                nemesh.dy, nemesh.z_spacing[i_target], nemesh.dx);
            N[idx + 4] = newell_nonequi::Nxy(
                y, z, x, nemesh.dy, nemesh.z_spacing[i_source], nemesh.dx,
                nemesh.dy, nemesh.z_spacing[i_target], nemesh.dx);
            N[idx + 5] = newell_nonequi::Nxx(
                z, x, y, nemesh.z_spacing[i_source], nemesh.dx, nemesh.dy,
                nemesh.z_spacing[i_target], nemesh.dx, nemesh.dy);
          }
        }
      }
    }
  }
}

std::vector<double> calculate_N(const NonequiMesh &nemesh,
                                std::size_t nthreads_in) {
  std::size_t nthreads =
      nthreads_in > 0 ? nthreads_in : std::thread::hardware_concurrency();
  assert(nthreads > 0);

  std::vector<double> N_values(nx_expanded(nemesh.nx) * ny_expanded(nemesh.ny) *
                               (nemesh.nz() * (nemesh.nz() + 1)) / 2 * 6);

  std::vector<std::thread> threads;
  for (std::size_t i = 0; i < nthreads; i++) {
    std::size_t ix_start =
        i * static_cast<double>(nx_expanded(nemesh.nx)) / nthreads;
    std::size_t ix_end =
        (i + 1) * static_cast<double>(nx_expanded(nemesh.nx)) / nthreads;
    threads.emplace_back(init_N, std::ref(nemesh), std::ref(N_values), ix_start,
                         ix_end);
  }

  for (auto &thread : threads) {
    thread.join();
  }

  return N_values;
}

  // TODO add option returning fft transform if arrayfire is installed

  // af::array Naf(6, (nemesh.nz() * (nemesh.nz() +
  // 1)) / 2, ny_expanded(nemesh.ny), nx_expanded(nemesh.nx),
  //               N_values.data());
  // Naf = af::reorder(Naf, 3, 2, 1, 0);
  // Naf = af::fftR2C<2>(Naf);
  // return Naf;

} // namespace newell_nonequi
} // namespace magnumaf

namespace bn = boost::python::numpy;
namespace bp = boost::python;

// const auto vec_from_pylist = []<typename T>(bp::list &values) {
template <typename T> std::vector<T> vec_from_pylist(bp::list &values) {
  std::vector<T> result{};
  for (bp::ssize_t i = 0; i < bp::len(values); ++i) {
    result.push_back(bp::extract<T>(values[i]));
  }
  return result;
}

bn::ndarray cp_vec_to_ndarray(std::vector<double> const &v,
                              std::vector<Py_intptr_t> const &shape) {
  auto result =
      bn::empty(shape.size(), shape.data(), bn::dtype::get_builtin<double>());
  std::copy(v.begin(), v.end(), reinterpret_cast<double *>(result.get_data()));
  return result;
}

template <typename From, typename To>
std::vector<To> statcast_each(std::initializer_list<From> &&from) {
  std::vector<To> to;
  for (auto elem : from) {
    to.emplace_back(static_cast<To>(elem));
  }
  return to;
}

bn::ndarray setup_demagtensor(size_t nx, size_t ny, double dx, double dy,
                              bp::list &dz_list) {

  auto dz_values = vec_from_pylist<double>(dz_list);
  // auto dz_values = vec_from_pylist.operator()<double>(dz_list);

  const auto nemesh = magnumaf::NonequiMesh{nx, ny, dx, dy, dz_values};
  auto N = magnumaf::newell_nonequi::calculate_N(
      nemesh, std::thread::hardware_concurrency());

  const auto nx_exp = magnumaf::newell_nonequi::nx_expanded(nemesh.nx);
  const auto ny_exp = magnumaf::newell_nonequi::ny_expanded(nemesh.ny);
  const auto nz = nemesh.nz();

  // std::vector<Py_intptr_t> shape = {
  //     nx_exp, ny_exp, static_cast<Py_intptr_t>(nz * (nz + 1)) / 2, 6};

  std::vector<Py_intptr_t> shape = statcast_each<std::size_t, Py_intptr_t>(
      {nx_exp, ny_exp, (nz * (nz + 1)) / 2, 6});

  const auto result = cp_vec_to_ndarray(N, shape);
  return result;
}

BOOST_PYTHON_MODULE(demag_nonequi_layers) {
  bn::initialize(); // NOTE: This is essential, otherwise returning bn::ndarray
                    // segfaults!
  bp::def("setup_demagtensor", setup_demagtensor,
          (bp::arg("nx"), bp::arg("ny"), bp::arg("dx"), bp::arg("dy"),
           bp::arg("dz_list")));
}
