//******************************************************************************
///  @file siesta_quantities.hpp
///  @brief Contains a context to store the siesta quantities.
//******************************************************************************

#ifndef siesta_quantities_hpp
#define siesta_quantities_hpp

#include <string>
#include <netcdf.h>

#include <iostream>

#include "grid_quantity.hpp"

//------------------------------------------------------------------------------
///  @brief Siesta quantities.
///
///  Class representing a quantities from siesta.
//------------------------------------------------------------------------------
class siesta_quantities {
public:
///  r
    const siesta_fourier<full_grid, cosine> r;
///  z
    const siesta_fourier<full_grid, sine> z;
///  Radial derivative of toroidal flux.
    const siesta_grid<full_grid> phipf;
///  Radial derivative of poloidal flux.
    const siesta_grid<full_grid> chipf;
///  JB^s
    const siesta_fourier<half_grid, sine> jbsups;
///  JB^u
    const siesta_fourier<half_grid, cosine> jbsupu;
///  JB^v
    const siesta_fourier<half_grid, cosine> jbsupv;
///  JB^s
    const siesta_fourier<half_grid, sine> bsups;
///  JB^u
    const siesta_fourier<half_grid, cosine> bsupu;
///  JB^v
    const siesta_fourier<half_grid, cosine> bsupv;

//------------------------------------------------------------------------------
///  @brief Siesta quantities.
///
///  For now this only supports symmetric cases only.
///
///  @param[in] restart_file File name of a restart file.
//------------------------------------------------------------------------------
    siesta_quantities(const std::string &restart_file) :
    r(load_fourier<full_grid, cosine> (restart_file, "rmnc_m_n_r_")),
    z(load_fourier<full_grid, sine> (restart_file, "zmns_m_n_r_")),
    phipf(load<full_grid> (restart_file, "phipf_r_")),
    chipf(load<full_grid> (restart_file, "chipf_r_")),
    jbsups(load_fourier_denorm<half_grid, sine> (restart_file, "JBsupssh_m_n_r_")),
    jbsupu(load_fourier_denorm<half_grid, cosine> (restart_file, "JBsupuch_m_n_r_")),
    jbsupv(load_fourier_denorm<half_grid, cosine> (restart_file, "JBsupvch_m_n_r_")),
    bsups(load_fourier<half_grid, sine> (restart_file, "bsupsmnsh_m_n_r_")),
    bsupu(load_fourier<half_grid, cosine> (restart_file, "bsupumnch_m_n_r_")),
    bsupv(load_fourier<half_grid, cosine> (restart_file, "bsupvmnch_m_n_r_")) {}

//------------------------------------------------------------------------------
///  @brief Factory method to load a siesta quantity.
///
///  @param[in] restart_file File name of a siesta restart file.
///  @param[in] name          Name of the quantity to load.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS>
    static siesta_grid<GIRD_CLASS> load(const std::string &restart_file,
                                        const std::string &name) {
        int ncid;
        nc_open(restart_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int nrad;
        nc_inq_varid(ncid, "nrad", &varid);
        nc_get_var(ncid, varid, &nrad);

        std::vector<double> temp(nrad);
        nc_inq_varid(ncid, name.c_str(), &varid);
        nc_get_var(ncid, varid, temp.data());

        nc_close(ncid);

        return siesta_grid<GIRD_CLASS> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a fourier siesta quantity.
///
///  @param[in] restart_file File name of a siesta restart file.
///  @param[in] name         Name of the quantity to load.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static siesta_fourier<GIRD_CLASS, PARITY> load_fourier(const std::string &restart_file,
                                                           const std::string &name) {
        int ncid;
        nc_open(restart_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int nrad;
        nc_inq_varid(ncid, "nrad", &varid);
        nc_get_var(ncid, varid, &nrad);

        int mpol;
        nc_inq_varid(ncid, "mpol", &varid);
        nc_get_var(ncid, varid, &mpol);

        int ntor;
        nc_inq_varid(ncid, "ntor", &varid);
        nc_get_var(ncid, varid, &ntor);

        int nfp;
        nc_inq_varid(ncid, "nfp", &varid);
        nc_get_var(ncid, varid, &nfp);

        siesta_quantity<GIRD_CLASS> quantity;
        for (size_t m = 0, em = mpol; m <= em; m++) {
            for (ssize_t n = -ntor, en = ntor; n <= en; n++) {
                const size_t ni = n + ntor;

                std::vector<double> temp(nrad);
                nc_inq_varid(ncid, name.c_str(), &varid);

                const std::array<size_t, 3> start = {0, ni, m};
                const std::array<size_t, 3> end = {static_cast<size_t> (nrad), 1, 1};
                nc_get_vara(ncid, varid,
                            start.data(), end.data(),
                            temp.data());

                if (n == -ntor) {
                    siesta_quantity_1d<GIRD_CLASS> new_grid;
                    new_grid.push_back(temp);
                    quantity.push_back(new_grid);
                } else {
                    quantity[m].push_back(siesta_grid<GIRD_CLASS> (temp));
                }
            }
        }

        nc_close(ncid);

        return siesta_fourier<GIRD_CLASS, PARITY> (quantity, mpol, ntor, nfp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a fourier denormalized siesta quantity.
///
///  @param[in] restart_file File name of a siesta restart file.
///  @param[in] name         Name of the quantity to load.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static siesta_fourier<GIRD_CLASS, PARITY> load_fourier_denorm(const std::string &restart_file,
                                                                  const std::string &name) {
        int ncid;
        nc_open(restart_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int nrad;
        nc_inq_varid(ncid, "nrad", &varid);
        nc_get_var(ncid, varid, &nrad);

        int mpol;
        nc_inq_varid(ncid, "mpol", &varid);
        nc_get_var(ncid, varid, &mpol);

        int ntor;
        nc_inq_varid(ncid, "ntor", &varid);
        nc_get_var(ncid, varid, &ntor);

        int nfp;
        nc_inq_varid(ncid, "nfp", &varid);
        nc_get_var(ncid, varid, &nfp);

        double bfactor;
        nc_inq_varid(ncid, "b_factor", &varid);
        nc_get_var(ncid, varid, &bfactor);

        siesta_quantity<GIRD_CLASS> quantity;
        for (size_t m = 0, em = mpol; m <= em; m++) {
            for (ssize_t n = -ntor, en = ntor; n <= en; n++) {
                const size_t ni = n + ntor;

                std::vector<double> temp(nrad);
                nc_inq_varid(ncid, name.c_str(), &varid);

                const std::array<size_t, 3> start = {0, ni, m};
                const std::array<size_t, 3> end = {static_cast<size_t> (nrad), 1, 1};
                nc_get_vara(ncid, varid,
                            start.data(), end.data(),
                            temp.data());

                for (size_t i = 0; i < nrad; i++) {
                    if (m == 0 && n == 0) {
                        temp[i] = temp[i]/bfactor;
                    } else {
                        temp[i] = sqrt(2.0)*temp[i]/bfactor;
                    }
                }

                if (n == -ntor) {
                    siesta_quantity_1d<GIRD_CLASS> new_grid;
                    new_grid.push_back(temp);
                    quantity.push_back(new_grid);
                } else {
                    quantity[m].push_back(siesta_grid<GIRD_CLASS> (temp));
                }
            }
        }

        nc_close(ncid);

        return siesta_fourier<GIRD_CLASS, PARITY> (quantity, mpol, ntor, nfp);
    }

//------------------------------------------------------------------------------
///  @brief Get R quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns R at the s,u,v position.
//------------------------------------------------------------------------------
    double get_r(const double s, const double u, const double v) const {
        return r.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get R' quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns R' at the s,u,v position.
//------------------------------------------------------------------------------
    double get_r_prime(const double s, const double u, const double v) const {
        return r.get_prime(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dRdu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dRdu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dr_du(const double s, const double u, const double v) const {
        return r.get_du(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dRdv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dRdv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dr_dv(const double s, const double u, const double v) const {
        return r.get_dv(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get Z quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns Z at the s,u,v position.
//------------------------------------------------------------------------------
    double get_z(const double s, const double u, const double v) const {
        return z.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get Z' quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns R' at the s,u,v position.
//------------------------------------------------------------------------------
    double get_z_prime(const double s, const double u, const double v) const {
        return z.get_prime(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dZdu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dZdu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dz_du(const double s, const double u, const double v) const {
        return z.get_du(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dZdv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dZdv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dz_dv(const double s, const double u, const double v) const {
        return z.get_dv(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get phipf quantity.
///
///  @param[in] s Radial position.
///  @returns phipf at the s position.
//------------------------------------------------------------------------------
    double get_phipf(const double s) const {
        return phipf.get(s);
    }

//------------------------------------------------------------------------------
///  @brief Get chipf quantity.
///
///  @param[in] s Radial position.
///  @returns chipf at the s position.
//------------------------------------------------------------------------------
    double get_chipf(const double s) const {
        return chipf.get(s);
    }

//------------------------------------------------------------------------------
///  @brief Get jbsups quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jbsups at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jbsups(const double s, const double u, const double v) const {
        return jbsups.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get jbsupu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jbsupu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jbsupu(const double s, const double u, const double v) const {
        return jbsupu.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get jbsupv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jbsupv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jbsupv(const double s, const double u, const double v) const {
        return jbsupv.get(s, u, v);
    }
};

#endif
