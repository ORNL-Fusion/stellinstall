#include <vector>
#include <iostream>
#include <random>
#include <chrono>

#include "vmec_quantities.hpp"
#include "siesta_quantities.hpp"
#include "commandline_parser.hpp"

#define RUN_TEST_3D(vmec, siesta, method, name, min, max, tol)(   \
test([&vmec]   (const double s, const double u, const double v) { \
        return vmec.method(s,u,v);                                \
     },                                                           \
     [&siesta] (const double s, const double u, const double v) { \
        return siesta.method(s,u,v);                              \
     },                                                           \
     name, tol, min, max)                                         \
)

#define RUN_TEST_1D(vmec, siesta, method, name, min, max, tol)( \
test([&vmec]   (const double s) {return vmec.method(s);},       \
     [&siesta] (const double s) {return siesta.method(s);},     \
     name, tol, min, max)                                       \
)

//------------------------------------------------------------------------------
///  @brief Dump vmec and siesta profile quantities.
///
///  @param[in] vmec   Vmec quantities.
///  @param[in] siesta Siesta quantites.
///  @param[in] u      Poloidal angle.
///  @param[in] v      Toroidal angle.
//------------------------------------------------------------------------------
void dump(const vmec_quantities &vmec,
          const siesta_quantities &siesta,
          const double u,
          const double v) {
    std::cout << std::endl;
    for (double s = 0.0; s <= 1; s += 0.005) {
        std::cout << s << " ";
        std::cout << vmec.get_r(s, u, v) << " ";
        std::cout << siesta.get_r(s, u, v) << " ";
        std::cout << vmec.get_r_prime(s, u, 0.0) << " ";
        std::cout << siesta.get_r_prime(s, u, 0.0) << " ";
        std::cout << vmec.get_z(s, u, 0.0) << " ";
        std::cout << siesta.get_z(s, u, 0.0) << " ";
        std::cout << vmec.get_z_prime(s, u, 0.0) << " ";
        std::cout << siesta.get_z_prime(s, u, 0.0) << " ";
        std::cout << vmec.get_phipf(s) << " ";
        std::cout << siesta.get_phipf(s) << " ";
        std::cout << vmec.get_chipf(s) << " ";
        std::cout << siesta.get_chipf(s) << " ";
        std::cout << siesta.get_jbsups(s, u, 0.0) << " ";
        std::cout << siesta.get_jbsups(s, u, 0.0) << " ";
        std::cout << vmec.get_jbsupu(s, u, 0.0) << " ";
        std::cout << siesta.get_jbsupu(s, u, 0.0) << " ";
        std::cout << vmec.get_jbsupv(s, u, 0.0) << " ";
        std::cout << siesta.get_jbsupv(s, u, 0.0) << " ";
        std::cout << std::endl;
    }
}

//------------------------------------------------------------------------------
///  @brief Test comparitor.
///
///  Randomly samples s,u,v positions to check that the value fall with in the
///  tolerance value.
///
///  @param[in] vmec      Getter function for the vmec value.
///  @param[in] siesta    Getter function for the siesta value.
///  @param[in] name      Name of the test.
///  @param[in] tolarance Tolarnce threshold that values must fall with in.
///  @param[in] smin      Minimum radial extent.
///  @param[in] smax      Maximum radial extent.
//------------------------------------------------------------------------------
void test(std::function<double(const double,const double,const double)> vmec,
          std::function<double(const double,const double,const double)> siesta,
          const std::string &name,
          const double tolarance,
          const double smin,
          const double smax) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> min_to_max(smin, smax);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = min_to_max(engine);
        const double u = zero_to_twopi(engine);
        const double v = zero_to_twopi(engine);
        double vmec_value = vmec(s, u, v);
        double siesta_value = siesta(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        const bool pass = diff < tolarance;

        if (!pass) {
            std::cout << name << ": failed at s = " << s
                      <<                    " u = " << u
                      <<                    " v = " << v << std::endl;
            std::cout <<            "vmec value = " << vmec_value
                      <<         " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.00001) {
                vmec_value = vmec(s, u, v);
                siesta_value = siesta(s, u, v);
                diff = abs(vmec_value - siesta_value);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << diff << std::endl;
            }

            exit(1);
        }
    }

    std::cout << name << " : passed" << std::endl;
}

//------------------------------------------------------------------------------
///  @brief Test comparitor.
///
///  Randomly samples s positions to check that the value fall with in the
///  tolerance value.
///
///  @param[in] vmec      Getter function for the vmec value.
///  @param[in] siesta    Getter function for the siesta value.
///  @param[in] name      Name of the test.
///  @param[in] tolarance Tolarnce threshold that values must fall with in.
///  @param[in] smin      Minimum radial extent.
///  @param[in] smax      Maximum radial extent.
//------------------------------------------------------------------------------
void test(std::function<double(const double)> vmec,
          std::function<double(const double)> siesta,
          const std::string &name,
          const double tolarance,
          const double smin,
          const double smax) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> min_to_max(smin, smax);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = min_to_max(engine);
        double vmec_value = vmec(s);
        double siesta_value = siesta(s);
        double diff = abs(vmec_value - siesta_value);

        const bool pass = diff < tolarance;

        if (!pass) {
            std::cout << name << ": failed at s = " << std::endl;
            std::cout <<            "vmec value = " << vmec_value
                      <<         " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.00001) {
                vmec_value = vmec(s);
                siesta_value = siesta(s);
                diff = abs(vmec_value - siesta_value);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << diff << std::endl;
            }

            exit(1);
        }
    }

    std::cout << name << " : passed" << std::endl;
}

//------------------------------------------------------------------------------
///  @brief Main test program.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Array of arguments strings.
//------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    const commandline_parser args(argc, argv);

    const vmec_quantities vmec(args.get("-wout_file"));
    const siesta_quantities siesta(args.get("-restart_file"));

    RUN_TEST_3D(vmec, siesta, get_r,     "test_r1",     0.0,  0.01, 1.0E-2);
    RUN_TEST_3D(vmec, siesta, get_r,     "test_r2",     0.01, 0.04, 1.0E-3);
    RUN_TEST_3D(vmec, siesta, get_r,     "test_r3",     0.04, 0.2,  1.0E-4);
    RUN_TEST_3D(vmec, siesta, get_r,     "test_r4",     0.2,  1.0,  1.0E-5);

    RUN_TEST_3D(vmec, siesta, get_r_prime, "test_r_prime", 0.0, 1.0, 1.0E-5);

    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_du1", 0.0,  0.01, 1.0E-2);
    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_du2", 0.01, 0.04, 1.0E-3);
    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_du3", 0.04, 0.2,  1.0E-4);
    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_du4", 0.2,  1.0,  1.0E-5);

    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_dv1", 0.0,  0.01, 1.0E-2);
    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_dv2", 0.01, 0.04, 1.0E-3);
    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_dv3", 0.04, 0.2,  1.0E-4);
    RUN_TEST_3D(vmec, siesta, get_dr_du, "test_dr_dv4", 0.2,  1.0,  1.0E-5);

    RUN_TEST_3D(vmec, siesta, get_z,     "test_z1",     0.0,  0.01, 1.0E-2);
    RUN_TEST_3D(vmec, siesta, get_z,     "test_z2",     0.01, 0.04, 1.0E-3);
    RUN_TEST_3D(vmec, siesta, get_z,     "test_z3",     0.04, 0.2,  1.0E-4);
    RUN_TEST_3D(vmec, siesta, get_z,     "test_z4",     0.2,  1.0,  1.0E-5);

//    RUN_TEST_3D(vmec, siesta, get_z_prime, "test_z_prime", 0.0, 1.0, 1.0E-5);

    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_du1", 0.0,  0.01, 1.0E-2);
    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_du2", 0.01, 0.04, 1.0E-3);
    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_du3", 0.04, 0.2,  1.0E-4);
    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_du4", 0.2,  1.0,  1.0E-5);

    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_dv1", 0.0,  0.01, 1.0E-2);
    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_dv2", 0.01, 0.04, 1.0E-3);
    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_dv3", 0.04, 0.2,  1.0E-4);
    RUN_TEST_3D(vmec, siesta, get_dz_du, "test_dz_dv4", 0.2,  1.0,  1.0E-5);

    RUN_TEST_1D(vmec, siesta, get_phipf,   "test_phipf",   0.0, 1.0, 1.0E-17);
    RUN_TEST_1D(vmec, siesta, get_chipf,   "test_chipf",   0.0, 1.0, 2.0E-6);

//    RUN_TEST_3D(vmec, siesta, get_jbsups,  "test_jbsups",  0.0, 1.0, 1.0E-5);
    RUN_TEST_3D(vmec, siesta, get_jbsupu,  "test_jbsupu",  0.0, 1.0, 1.0E-4);
    RUN_TEST_3D(vmec, siesta, get_jbsupv,  "test_jbsupv",  0.0, 1.0, 1.0E-4);

//    dump(vmec, siesta, 0.0, 0.0);
//    dump(vmec, siesta, M_PI/2.0, 0.0);
}
