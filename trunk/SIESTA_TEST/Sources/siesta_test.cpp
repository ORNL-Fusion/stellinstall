#include <vector>
#include <iostream>
#include <random>
#include <chrono>

#include "vmec_quantities.hpp"
#include "siesta_quantities.hpp"
#include "commandline_parser.hpp"

void dump(const vmec_quantities &vmec,
          const siesta_quantities &siesta,
          const double u) {
    std::cout << std::endl;
    for (double s = 0.0; s <= 1; s += 0.005) {
        std::cout << s << " ";
        std::cout << vmec.r.get(s, u, 0.0) << " ";
        std::cout << siesta.r.get(s, u, 0.0) << " ";
        std::cout << vmec.r.get_prime(s, u, 0.0) << " ";
        std::cout << siesta.r.get_prime(s, u, 0.0) << " ";
        std::cout << vmec.signj*sqrt(s)*vmec.phipf.get(s)/M_PI << " ";
        std::cout << siesta.phipf.get(s) << " ";
        std::cout << vmec.signj*sqrt(s)*vmec.chipf.get(s)/M_PI << " ";
        std::cout << siesta.chipf.get(s) << " ";
        std::cout << siesta.jbsups.get(s, u, 0.0) << " ";
        std::cout << 2.0*sqrt(s)*vmec.j.get(s, u, 0.0)*vmec.bsupu.get(s, u, 0.0) << " ";
        std::cout << siesta.jbsupu.get(s, u, 0.0) << " ";
        std::cout << 2.0*sqrt(s)*vmec.j.get(s, u, 0.0)*vmec.bsupv.get(s, u, 0.0) << " ";
        std::cout << siesta.jbsupv.get(s, u, 0.0) << " ";
        std::cout << siesta.bsups.get(s, u, 0.0) << " ";
        std::cout << vmec.bsupu.get(s, u, 0.0) << " ";
        std::cout << siesta.bsupu.get(s, u, 0.0) << " ";
        std::cout << vmec.bsupv.get(s, u, 0.0) << " ";
        std::cout << siesta.bsupv.get(s, u, 0.0) << " ";
        std::cout << std::endl;
    }
}

void test_r(const vmec_quantities &vmec,
            const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.r.get(s, u, v);
        double siesta_value = siesta.r.get(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s == 0.0) {
            pass = diff < 1.0E-15;
        } else if (s > 0.1) {
            pass = diff < 2.0E-5;
        } else if (s > 0.01) {
            pass = diff < 2.0E-3;
        } else {
            pass = diff < 1.0E-1;
        }

        if (!pass) {
            std::cout << "test_r: failed at s = " << s
                      <<                  " u = " << u
                      <<                  " v = " << v << std::endl;
            std::cout <<    "vmec value = " << vmec_value
                      << " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.r.get(s, u, v);
                siesta_value = siesta.r.get(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_r : passed" << std::endl;
}

void test_r_prime(const vmec_quantities &vmec,
                  const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.r.get_prime(s, u, v);
        double siesta_value = siesta.r.get_prime(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s > 0.3 && s < 0.99) {
            pass = diff < 3.0E-5;
        } else if (s > 0.2 && s <= 0.3) {
            pass = diff < 7.0E-5;
        }

        if (!pass) {
            std::cout << "test_r_prime: failed at s = " << s
                      <<                        " u = " << u
                      <<                        " v = " << v << std::endl;
            std::cout <<                "vmec value = " << vmec_value
                      <<             " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.r.get_prime(s, u, v);
                siesta_value = siesta.r.get_prime(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_r_prime : passed" << std::endl;
}

void test_dr_du(const vmec_quantities &vmec,
            const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.r.get_du(s, u, v);
        double siesta_value = siesta.r.get_du(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s == 0.0) {
            pass = diff < 1.0E-15;
        } else if (s > 0.1) {
            pass = diff < 2.0E-5;
        } else if (s > 0.01) {
            pass = diff < 2.0E-3;
        } else {
            pass = diff < 1.0E-1;
        }

        if (!pass) {
            std::cout << "test_dr_du: failed at s = " << s
                      <<                      " u = " << u
                      <<                      " v = " << v << std::endl;
            std::cout <<              "vmec value = " << vmec_value
                      <<           " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.r.get_du(s, u, v);
                siesta_value = siesta.r.get_du(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_dr_du : passed" << std::endl;
}

void test_dr_dv(const vmec_quantities &vmec,
                const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.r.get_dv(s, u, v);
        double siesta_value = siesta.r.get_dv(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s == 0.0) {
            pass = diff < 1.0E-15;
        } else if (s > 0.1) {
            pass = diff < 2.0E-5;
        } else if (s > 0.01) {
            pass = diff < 2.0E-3;
        } else {
            pass = diff < 1.0E-1;
        }

        if (!pass) {
            std::cout << "test_dr_dv: failed at s = " << s
                      <<                      " u = " << u
                      <<                      " v = " << v << std::endl;
            std::cout <<              "vmec value = " << vmec_value
                      <<           " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.r.get_dv(s, u, v);
                siesta_value = siesta.r.get_dv(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_dr_dv : passed" << std::endl;
}

void test_z(const vmec_quantities &vmec,
            const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.z.get(s, u, v);
        double siesta_value = siesta.z.get(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s == 0.0) {
            pass = diff < 1.0E-15;
        } else if (s > 0.1) {
            pass = diff < 2.0E-5;
        } else if (s > 0.01) {
            pass = diff < 2.0E-3;
        } else {
            pass = diff < 1.0E-1;
        }

        if (!pass) {
            std::cout << "test_z: failed at s = " << s
                      <<                  " u = " << u
                      <<                  " v = " << v << std::endl;
            std::cout <<    "vmec value = " << vmec_value
                      << " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.z.get(s, u, v);
                siesta_value = siesta.z.get(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_z : passed" << std::endl;
}

void test_z_prime(const vmec_quantities &vmec,
                  const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.z.get_prime(s, u, v);
        double siesta_value = siesta.z.get_prime(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s > 0.3 && s < 0.99) {
            pass = diff < 3.0E-5;
        } else if (s > 0.2 && s <= 0.3) {
            pass = diff < 7.0E-5;
        }

        if (!pass) {
            std::cout << "test_z_prime: failed at s = " << s
                      <<                        " u = " << u
                      <<                        " v = " << v << std::endl;
            std::cout <<                "vmec value = " << vmec_value
                      <<             " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.z.get_prime(s, u, v);
                siesta_value = siesta.z.get_prime(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_z_prime : passed" << std::endl;
}

void test_dz_du(const vmec_quantities &vmec,
            const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.z.get_du(s, u, v);
        double siesta_value = siesta.z.get_du(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s == 0.0) {
            pass = diff < 1.0E-15;
        } else if (s > 0.1) {
            pass = diff < 2.0E-5;
        } else if (s > 0.01) {
            pass = diff < 2.0E-3;
        } else {
            pass = diff < 1.0E-1;
        }

        if (!pass) {
            std::cout << "test_dz_du: failed at s = " << s
                      <<                      " u = " << u
                      <<                      " v = " << v << std::endl;
            std::cout <<              "vmec value = " << vmec_value
                      <<           " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.z.get_du(s, u, v);
                siesta_value = siesta.z.get_du(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_dz_du : passed" << std::endl;
}

void test_dz_dv(const vmec_quantities &vmec,
                const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = vmec.z.get_dv(s, u, v);
        double siesta_value = siesta.z.get_dv(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = true;
        if (s == 0.0) {
            pass = diff < 1.0E-15;
        } else if (s > 0.1) {
            pass = diff < 2.0E-5;
        } else if (s > 0.01) {
            pass = diff < 2.0E-3;
        } else {
            pass = diff < 1.0E-1;
        }

        if (!pass) {
            std::cout << "test_dz_dv: failed at s = " << s
                      <<                      " u = " << u
                      <<                      " v = " << v << std::endl;
            std::cout <<              "vmec value = " << vmec_value
                      <<           " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.z.get_dv(s, u, v);
                siesta_value = siesta.z.get_dv(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_dz_dv : passed" << std::endl;
}

void test_phipf(const vmec_quantities &vmec,
                const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double vmec_value = -sqrt(s)*vmec.phipf.get(s)/M_PI;
        double siesta_value = siesta.phipf.get(s);
        double diff = abs(vmec_value - siesta_value);
        bool pass = diff < 1.0E-10;

        if (!pass) {
            std::cout << "test_phipf: failed at s = " << s << std::endl;
            std::cout <<              "vmec value = " << vmec_value
                      <<           " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.signj*vmec.signj*sqrt(s)*vmec.phipf.get(s)/M_PI;
                siesta_value = siesta.phipf.get(s);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_phipf : passed" << std::endl;
}

void test_chipf(const vmec_quantities &vmec,
                const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double vmec_value = vmec.signj*sqrt(s)*vmec.chipf.get(s)/M_PI;
        double siesta_value = siesta.chipf.get(s);
        double diff = abs(vmec_value - siesta_value);

        bool pass = diff < 1.0E-5;

        if (!pass) {
            std::cout << "test_chipf: failed at s = " << s << std::endl;
            std::cout <<              "vmec value = " << vmec_value
                      <<           " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = vmec.signj*sqrt(s)*vmec.phipf.get(s)/M_PI;
                siesta_value = siesta.phipf.get(s);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_chipf : passed" << std::endl;
}

void test_jbsups(const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double siesta_value = siesta.jbsups.get(s, u, v);
        double diff = abs(siesta_value);

        bool pass = siesta_value < 1.0E-15;

        if (!pass) {
            std::cout << "test_jbsups: failed at s = " << s << std::endl;
            std::cout <<            " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                siesta_value = siesta.jbsups.get(s, u, v);
                std::cout << s << " "
                          << siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_jbsups : passed" << std::endl;
}

void test_jbsupu(const vmec_quantities &vmec,
                 const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = 2.0*sqrt(s)*vmec.j.get(s, u, v)*vmec.bsupu.get(s, u, v);
        double siesta_value = siesta.jbsupu.get(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = diff < 2.0E-4;

        if (!pass) {
            std::cout << "test_jbsupu: failed at s = " << s
                      <<                       " u = " << u
                      <<                       " v = " << v << std::endl;
            std::cout <<               "vmec value = " << vmec_value
                      <<            " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = 2.0*sqrt(s)*vmec.j.get(s, u, v)*vmec.bsupu.get(s, u, v);
                siesta_value = siesta.jbsupu.get(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_jbsupu : passed" << std::endl;
}

void test_jbsupv(const vmec_quantities &vmec,
                 const siesta_quantities &siesta) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> zero_to_one(0.0, 1.0);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = zero_to_one(engine);
        double u = zero_to_twopi(engine);
        double v = zero_to_twopi(engine);
        double vmec_value = 2.0*sqrt(s)*vmec.j.get(s, u, v)*vmec.bsupv.get(s, u, v);
        double siesta_value = siesta.jbsupv.get(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        bool pass = diff < 2.0E-4;

        if (!pass) {
            std::cout << "test_jbsupv: failed at s = " << s
                      <<                       " u = " << u
                      <<                       " v = " << v << std::endl;
            std::cout <<               "vmec value = " << vmec_value
                      <<            " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.005) {
                vmec_value = 2.0*sqrt(s)*vmec.j.get(s, u, v)*vmec.bsupv.get(s, u, v);
                siesta_value = siesta.jbsupv.get(s, u, v);
                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << vmec_value - siesta_value << std::endl;
            }

            exit(1);
        }
    }

    std::cout << "test_jbsupv : passed" << std::endl;
}

int main(int argc, const char * argv[]) {
    const commandline_parser args(argc, argv);

    const vmec_quantities vmec(args.get("-wout_file"));
    const siesta_quantities siesta(args.get("-restart_file"));

    test_r(vmec, siesta);
    test_r_prime(vmec, siesta);
    test_dr_du(vmec, siesta);
    test_dr_dv(vmec, siesta);

    test_z(vmec, siesta);
    test_z_prime(vmec, siesta);
    test_dz_du(vmec, siesta);
    test_dz_dv(vmec, siesta);

    test_phipf(vmec, siesta);
    test_chipf(vmec, siesta);

    test_jbsups(siesta);
    test_jbsupu(vmec, siesta);
    test_jbsupv(vmec, siesta);

    dump(vmec, siesta, 0.0);
    dump(vmec, siesta, M_PI/2.0);
}
