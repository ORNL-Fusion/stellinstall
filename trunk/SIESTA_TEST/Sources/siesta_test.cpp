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

#define RUN_TEST_0D(vmec, siesta, method, name, tol)( \
test([&vmec]   () {return vmec.method();},            \
     [&siesta] () {return siesta.method();},          \
     name, tol)                                       \
)

#define RUN_TEST_DIV(eq, method, name, min, max, tol)(        \
test([&eq] (const double s, const double u, const double v) { \
        return eq.method(s,u,v);                              \
     },                                                       \
     name, tol, min, max)                                     \
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
///  Randomly samples s,u,v positions to check that the value fall with in the
///  tolerance value.
///
///  @param[in] div       Getter function for the diveregnce value.
///  @param[in] name      Name of the test.
///  @param[in] tolarance Tolarnce threshold that values must fall with in.
///  @param[in] smin      Minimum radial extent.
///  @param[in] smax      Maximum radial extent.
//------------------------------------------------------------------------------
void test(std::function<double(const double,const double,const double)> div,
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
        double diff = abs(div(s, u, v));

        const bool pass = diff < tolarance;

        if (!pass) {
            std::cout << name << ": failed at s = " << s
                      <<                    " u = " << u
                      <<                    " v = " << v << std::endl;
            std::cout << " divergence = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.00001) {
                diff = abs(div(s, u, v));
                std::cout << s << " "
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
void test(std::function<double(void)> vmec,
          std::function<double(void)> siesta,
          const std::string &name,
          const double tolarance) {

    double vmec_value = vmec();
    double siesta_value = siesta();
    double diff = abs(vmec_value - siesta_value);

    const bool pass = diff < tolarance;

    if (!pass) {
        std::cout << name << ": failed " << std::endl;
        std::cout <<         "vmec value = " << vmec_value
                  <<         " siesta value = " << siesta_value
                  << " difference between vmec and siesta = "
                  << diff << std::endl;

        exit(1);
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

    const vmec_quantities vmec(args.get<std::string> ("-wout_file"));
    const siesta_quantities siesta(args.get<std::string> ("-restart_file"));

    const std::string test_name = args.get<std::string> ("-test_name");

    if (test_name == "r") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_r, test_name, min,  max, tol);
    } else if (test_name == "z") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_z, test_name, min,  max, tol);
    } else if (test_name == "drdu") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_dr_du, test_name, min,  max, tol);
    } else if (test_name == "drdv") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_dr_dv, test_name, min,  max, tol);
    } else if (test_name == "dzdu") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_dz_du, test_name, min,  max, tol);
    } else if (test_name == "dzdv") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_dz_dv, test_name, min,  max, tol);
    } else if (test_name == "phipf") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_1D(vmec, siesta, get_phipf, test_name, min,  max, tol);
    } else if (test_name == "chipf") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_1D(vmec, siesta, get_chipf, test_name, min,  max, tol);
    } else if (test_name == "jbsups") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_jbsups, test_name, min,  max, tol);
    } else if (test_name == "jbsupu") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_jbsupu, test_name, min,  max, tol);
    } else if (test_name == "jbsupv") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_jbsupv, test_name, min,  max, tol);
    } else if (test_name == "jacobian") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_3D(vmec, siesta, get_jacobian, test_name, min,  max, tol);
    } else if (test_name == "curtor") {
        const double tol = args.get<double> ("-tol");
        RUN_TEST_0D(vmec, siesta, get_curtor, test_name, tol);
    } else if (test_name == "vmec_divb") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_DIV(vmec, get_divb, test_name, min,  max, tol);
    } else if (test_name == "siesta_divb") {
        const double tol = args.get<double> ("-tol");
        const double min = args.get<double> ("-min");
        const double max = args.get<double> ("-max");
        RUN_TEST_DIV(siesta, get_divb, test_name, min,  max, tol);
    } else {
        std::cout << "Unknown test name " << test_name << "." << std::endl;
    }
}
