//******************************************************************************
///  @file vmec_quantities.hpp
///  @brief Contains a context to store the vmec quantities.
//******************************************************************************

#ifndef vmec_quantities_hpp
#define vmec_quantities_hpp

#include <string>
#include <netcdf.h>
#include <array>

#include "grid_quantity.hpp"

//------------------------------------------------------------------------------
///  @brief Vmec quantities.
///
///  Class representing a quantities from vmec.
//------------------------------------------------------------------------------
class vmec_quantities {
public:
///  r
    const vmec_fourier<full_grid, cosine> r;
///  z
    const vmec_fourier<full_grid, sine> z;
///  Radial derivative of toroidal flux.
    const vmec_grid<full_grid> phipf;
///  Radial derivative of poloidal flux.
    const vmec_grid<full_grid> chipf;
///  Sign of vmec jacobian.
    const double signj;
///  Jacobian
    const vmec_fourier<half_grid, cosine> j;
///  B^u
    const vmec_fourier<half_grid, cosine> bsupu;
///  B^v
    const vmec_fourier<half_grid, cosine> bsupv;

//------------------------------------------------------------------------------
///  @brief Vmec quantities.
///
///  For now this only supports symmetric cases only.
///
///  @param[in] wout_file File name of a vmec wout file.
//------------------------------------------------------------------------------
    vmec_quantities(const std::string &wout_file) :
    r(load_fourier<full_grid, cosine> (wout_file, "rmnc")),
    z(load_fourier<full_grid, sine> (wout_file, "zmns")),
    phipf(load<full_grid> (wout_file, "phipf")),
    chipf(load<full_grid> (wout_file, "chipf")),
    signj(load_scalar<int> (wout_file, "signgs")),
    j(load_fourier_nyq<half_grid, cosine> (wout_file, "gmnc")),
    bsupu(load_fourier_nyq<half_grid, cosine> (wout_file, "bsupumnc")),
    bsupv(load_fourier_nyq<half_grid, cosine> (wout_file, "bsupvmnc")) {}

//------------------------------------------------------------------------------
///  @brief Factory method to load a scalar vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the quantity to load.
//------------------------------------------------------------------------------
    template<typename TYPE>
    static double load_scalar(const std::string &wout_file,
                              const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        TYPE temp;
        nc_inq_varid(ncid, name.c_str(), &varid);
        nc_get_var(ncid, varid, &temp);

        nc_close(ncid);

        return static_cast<double> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the quantity to load.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS>
    static vmec_grid<GIRD_CLASS> load(const std::string &wout_file,
                                      const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int ns;
        nc_inq_varid(ncid, "ns", &varid);
        nc_get_var(ncid, varid, &ns);

        std::vector<double> temp(ns);
        nc_inq_varid(ncid, name.c_str(), &varid);
        nc_get_var(ncid, varid, temp.data());

        nc_close(ncid);

        return vmec_grid<GIRD_CLASS> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a fourier vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the quantity to load.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static vmec_fourier<GIRD_CLASS, PARITY> load_fourier(const std::string &wout_file,
                                                         const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int ns;
        nc_inq_varid(ncid, "ns", &varid);
        nc_get_var(ncid, varid, &ns);

        int mnmax;
        nc_inq_varid(ncid, "mnmax", &varid);
        nc_get_var(ncid, varid, &mnmax);

        std::vector<double> xm(mnmax);
        nc_inq_varid(ncid, "xm", &varid);
        nc_get_var(ncid, varid, xm.data());

        std::vector<double> xn(mnmax);
        nc_inq_varid(ncid, "xn", &varid);
        nc_get_var(ncid, varid, xn.data());

        vmec_quantity<GIRD_CLASS> quantity;
        for (size_t i = 0; i < mnmax; i++) {
            std::vector<double> temp(ns);
            nc_inq_varid(ncid, name.c_str(), &varid);

            const std::array<size_t, 2> start = {0, i};
            const std::array<size_t, 2> end = {static_cast<size_t> (ns), 1};
            nc_get_vara(ncid, varid,
                        start.data(), end.data(),
                        temp.data());

            quantity.push_back(vmec_grid<GIRD_CLASS> (temp));
        }

        nc_close(ncid);

        return vmec_fourier<GIRD_CLASS, PARITY> (quantity, xm, xn);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a nyquest fourier vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the nyquest quantity to load.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static vmec_fourier<GIRD_CLASS, PARITY> load_fourier_nyq(const std::string &wout_file,
                                                             const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int ns;
        nc_inq_varid(ncid, "ns", &varid);
        nc_get_var(ncid, varid, &ns);

        int mnmax;
        nc_inq_varid(ncid, "mnmax_nyq", &varid);
        nc_get_var(ncid, varid, &mnmax);

        std::vector<double> xm(mnmax);
        nc_inq_varid(ncid, "xm_nyq", &varid);
        nc_get_var(ncid, varid, xm.data());

        std::vector<double> xn(mnmax);
        nc_inq_varid(ncid, "xn_nyq", &varid);
        nc_get_var(ncid, varid, xn.data());

        vmec_quantity<GIRD_CLASS> quantity;
        for (size_t i = 0; i < mnmax; i++) {
            std::vector<double> temp(ns);
            nc_inq_varid(ncid, name.c_str(), &varid);

            const std::array<size_t, 2> start = {0, i};
            const std::array<size_t, 2> end = {static_cast<size_t> (ns), 1};
            nc_get_vara(ncid, varid,
                        start.data(), end.data(),
                        temp.data());

            quantity.push_back(vmec_grid<GIRD_CLASS> (temp));
        }

        nc_close(ncid);

        return vmec_fourier<GIRD_CLASS, PARITY> (quantity, xm, xn);
    }
};

#endif
