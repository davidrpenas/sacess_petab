#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_ode.h"
#include "amici/splinefunctions.h"

namespace amici {

class Solver;

namespace model_Oliveira_NatCommun2021 {

extern std::array<const char*, 14> parameterNames;
extern std::array<const char*, 11> fixedParameterNames;
extern std::array<const char*, 9> stateNames;
extern std::array<const char*, 2> observableNames;
extern std::array<const ObservableScaling, 2> observableScalings;
extern std::array<const char*, 22> expressionNames;
extern std::array<const char*, 14> parameterIds;
extern std::array<const char*, 11> fixedParameterIds;
extern std::array<const char*, 9> stateIds;
extern std::array<const char*, 2> observableIds;
extern std::array<const char*, 22> expressionIds;
extern std::array<int, 9> stateIdxsSolver;
extern std::array<bool, 4> rootInitialValues;

extern void Jy_Oliveira_NatCommun2021(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_Oliveira_NatCommun2021(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_Oliveira_NatCommun2021(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &rowvals, int index);






extern void root_Oliveira_NatCommun2021(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void dwdp_Oliveira_NatCommun2021(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static);
extern void dwdp_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &colptrs);
extern void dwdp_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &rowvals);
extern void dwdx_Oliveira_NatCommun2021(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static);
extern void dwdx_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &colptrs);
extern void dwdx_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &rowvals);
extern void dwdw_Oliveira_NatCommun2021(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static);
extern void dwdw_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &colptrs);
extern void dwdw_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &rowvals);
extern void dxdotdw_Oliveira_NatCommun2021(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &colptrs);
extern void dxdotdw_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &rowvals);



extern void dxdotdx_explicit_Oliveira_NatCommun2021(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdx_explicit_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &colptrs);
extern void dxdotdx_explicit_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &rowvals);
extern void dydx_Oliveira_NatCommun2021(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);





extern void sigmay_Oliveira_NatCommun2021(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);

extern void dsigmaydp_Oliveira_NatCommun2021(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip);


extern void w_Oliveira_NatCommun2021(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static);
extern void x0_Oliveira_NatCommun2021(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void x0_fixedParameters_Oliveira_NatCommun2021(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs);
extern void sx0_Oliveira_NatCommun2021(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);
extern void sx0_fixedParameters_Oliveira_NatCommun2021(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs);
extern void xdot_Oliveira_NatCommun2021(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_Oliveira_NatCommun2021(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);


extern void stau_Oliveira_NatCommun2021(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie);

extern void deltasx_Oliveira_NatCommun2021(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl);

extern void x_solver_Oliveira_NatCommun2021(realtype *x_solver, const realtype *x_rdata);












extern std::vector<HermiteSpline> create_splines_Oliveira_NatCommun2021(const realtype *p, const realtype *k);



/**
 * @brief AMICI-generated model subclass.
 */
class Model_Oliveira_NatCommun2021 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Oliveira_NatCommun2021()
        : amici::Model_ODE(
              amici::ModelDimensions(
                  9,                            // nx_rdata
                  9,                        // nxtrue_rdata
                  9,                           // nx_solver
                  9,                       // nxtrue_solver
                  0,                    // nx_solver_reinit
                  14,                                  // np
                  11,                                  // nk
                  2,                                  // ny
                  2,                              // nytrue
                  0,                                  // nz
                  0,                              // nztrue
                  4,                              // nevent
                  4,                       // nevent_solver
                  0,                                // nspl
                  1,                          // nobjective
                  22,                                  // nw
                  14,                               // ndwdx
                  12,                               // ndwdp
                  10,                               // ndwdw
                  25,                            // ndxdotdw
                  std::vector<int>{1,1},                              // ndjydy
                  0,                    // ndxrdatadxsolver
                  0,                        // ndxrdatadtcl
                  0,                        // ndtotal_cldx_rdata
                  0,                                       // nnz
                  9,                                 // ubw
                  9,                                 // lbw
                  true,                                    // pythonGenerated
                  0,                   // ndxdotdp_explicit
                  1,                   // ndxdotdx_explicit
                  1                    // w_recursion_depth
              ),
              amici::SimulationParameters(
                  std::vector<realtype>{12098867.0, 0.20000000000000001, 0.25, 0.28571400000000002, 0.25, 0.53000000000000003, 0.14999999999999999, 0.40000000000000002, 0.14000000000000001, 0.28999999999999998, 1.0}, // fixedParameters
                  std::vector<realtype>{0.69513610540000004, 0.0, 0.131351606436, 0.18702017558850001, 0.80719409095956296, 0.56064148877088604, 0.66000000000000003, 23.4530530945, 97.0, 2.2467317002508198e-6, 6.4828815498308996e-7, 1.9925913556037298e-6, 0.0, 0.0}        // dynamic parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},   // idlist
              std::vector<int>{},               // z2events
              {}               // state-independent events
          ) {
                 root_initial_values_ = std::vector<bool>(
                     rootInitialValues.begin(), rootInitialValues.end()
                 );
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_Oliveira_NatCommun2021(*this);
    }

    void fJrz(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_Oliveira_NatCommun2021(Jy, iy, p, k, y, sigmay, my);
    }


    void fJz(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_Oliveira_NatCommun2021(dJydsigma, iy, p, k, y, sigmay, my);
    }


    void fdJzdsigma(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz) override {}


    /**
     * @brief model specific implementation of fdeltasx
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     */
    void fdeltaqB(realtype *deltaqB, const realtype t,
                  const realtype *x, const realtype *p,
                  const realtype *k, const realtype *h, const int ip,
                  const int ie, const realtype *xdot,
                  const realtype *xdot_old,
                  const realtype *xB) override {}

    void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl) override {
        deltasx_Oliveira_NatCommun2021(deltasx, t, x, p, k, h, w, ip, ie, xdot, xdot_old, sx, stau, tcl);
    }


    void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {}


    /**
     * @brief model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     */
    void fdeltaxB(realtype *deltaxB, const realtype t,
                  const realtype *x, const realtype *p,
                  const realtype *k, const realtype *h, const int ie,
                  const realtype *xdot, const realtype *xdot_old,
                  const realtype *xB) override {}

    void fdrzdp(realtype *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdrzdx(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {
        dsigmaydp_Oliveira_NatCommun2021(dsigmaydp, t, p, k, y, ip);
    }


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_Oliveira_NatCommun2021(dJydy, iy, p, k, y, sigmay, my);
    }

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {
        dJydy_colptrs_Oliveira_NatCommun2021(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {
        dJydy_rowvals_Oliveira_NatCommun2021(rowvals, index);
    }


    std::vector<HermiteSpline> fcreate_splines(const realtype *p, const realtype *k) override {
        return create_splines_Oliveira_NatCommun2021(p, k);
    }

    void fdspline_valuesdp(realtype *dspline_valuesdp, const realtype *p, const realtype *k, const int ip) override {}

    void fdspline_slopesdp(realtype *dspline_slopesdp, const realtype *p, const realtype *k, const int ip) override {}


    void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static) override {
        dwdp_Oliveira_NatCommun2021(dwdp, t, x, p, k, h, w, tcl, dtcldp, spl, sspl, include_static);
    }

    void fdwdp_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdp_colptrs_Oliveira_NatCommun2021(colptrs);
    }

    void fdwdp_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdp_rowvals_Oliveira_NatCommun2021(rowvals);
    }


    void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static) override {
        dwdx_Oliveira_NatCommun2021(dwdx, t, x, p, k, h, w, tcl, spl, include_static);
    }

    void fdwdx_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdx_colptrs_Oliveira_NatCommun2021(colptrs);
    }

    void fdwdx_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdx_rowvals_Oliveira_NatCommun2021(rowvals);
    }


    void fdwdw(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static) override {
        dwdw_Oliveira_NatCommun2021(dwdw, t, x, p, k, h, w, tcl, include_static);
    }

    void fdwdw_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdw_colptrs_Oliveira_NatCommun2021(colptrs);
    }

    void fdwdw_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdw_rowvals_Oliveira_NatCommun2021(rowvals);
    }


    void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Oliveira_NatCommun2021(dxdotdw, t, x, p, k, h, w);
    }

    void fdxdotdw_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdw_colptrs_Oliveira_NatCommun2021(colptrs);
    }

    void fdxdotdw_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdw_rowvals_Oliveira_NatCommun2021(rowvals);
    }


    void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdx_explicit_Oliveira_NatCommun2021(dxdotdx_explicit, t, x, p, k, h, w);
    }

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdx_explicit_colptrs_Oliveira_NatCommun2021(colptrs);
    }

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdx_explicit_rowvals_Oliveira_NatCommun2021(rowvals);
    }


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        dydx_Oliveira_NatCommun2021(dydx, t, x, p, k, h, w, dwdx);
    }


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl) override {}


    void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        root_Oliveira_NatCommun2021(root, t, x, p, k, h, tcl);
    }


    void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_Oliveira_NatCommun2021(sigmay, t, p, k, y);
    }


    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) override {}


    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {
        stau_Oliveira_NatCommun2021(stau, t, x, p, k, h, tcl, sx, ip, ie);
    }

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {
        sx0_Oliveira_NatCommun2021(sx0, t, x, p, k, ip);
    }

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {
        sx0_fixedParameters_Oliveira_NatCommun2021(sx0_fixedParameters, t, x0, p, k, ip,  reinitialization_state_idxs);
    }


    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static) override {
        w_Oliveira_NatCommun2021(w, t, x, p, k, h, tcl, spl, include_static);
    }


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_Oliveira_NatCommun2021(x0, t, p, k);
    }


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {
        x0_fixedParameters_Oliveira_NatCommun2021(x0_fixedParameters, t, p, k,  reinitialization_state_idxs);
    }


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_Oliveira_NatCommun2021(xdot, t, x, p, k, h, w);
    }


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_Oliveira_NatCommun2021(y, t, x, p, k, h, w);
    }


    void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_Oliveira_NatCommun2021(x_solver, x_rdata);
    }


    void ftotal_cl(realtype *total_cl, const realtype *x_rdata, const realtype *p, const realtype *k) override {}


    void fdx_rdatadx_solver(realtype *dx_rdatadx_solver, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k) override {}

    void fdx_rdatadx_solver_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdx_rdatadx_solver_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdx_rdatadp(realtype *dx_rdatadp, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k, const int ip) override {}


    void fdx_rdatadtcl(realtype *dx_rdatadtcl, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k) override {}

    void fdx_rdatadtcl_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdx_rdatadtcl_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdtotal_cldp(realtype *dtotal_cldp, const realtype *x_rdata, const realtype *p, const realtype *k, const int ip) override {}


    void fdtotal_cldx_rdata(realtype *dtotal_cldx_rdata, const realtype *x_rdata, const realtype *p, const realtype *k, const realtype *tcl) override {}

    void fdtotal_cldx_rdata_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdtotal_cldx_rdata_rowvals(SUNMatrixWrapper &rowvals) override {}


    std::string getName() const override {
        return "Oliveira_NatCommun2021";
    }

    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>(parameterNames.begin(),
                                        parameterNames.end());
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>(stateNames.begin(), stateNames.end());
    }

    /**
     * @brief Get names of the solver states
     * @return the names
     */
    std::vector<std::string> getStateNamesSolver() const override {
        std::vector<std::string> result;
        result.reserve(stateIdxsSolver.size());
        for(auto idx: stateIdxsSolver) {
            result.push_back(stateNames[idx]);
        }
        return result;
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    std::vector<std::string> getFixedParameterNames() const override {
        return std::vector<std::string>(fixedParameterNames.begin(),
                                        fixedParameterNames.end());
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    std::vector<std::string> getObservableNames() const override {
        return std::vector<std::string>(observableNames.begin(),
                                        observableNames.end());
    }

    /**
     * @brief Get names of model expressions
     * @return Expression names
     */
    std::vector<std::string> getExpressionNames() const override {
        return std::vector<std::string>(expressionNames.begin(),
                                        expressionNames.end());
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>(parameterIds.begin(),
                                        parameterIds.end());
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>(stateIds.begin(), stateIds.end());
    }

    /**
     * @brief Get ids of the solver states
     * @return the ids
     */
    std::vector<std::string> getStateIdsSolver() const override {
        std::vector<std::string> result;
        result.reserve(stateIdxsSolver.size());
        for(auto idx: stateIdxsSolver) {
            result.push_back(stateIds[idx]);
        }
        return result;
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    std::vector<std::string> getFixedParameterIds() const override {
        return std::vector<std::string>(fixedParameterIds.begin(),
                                        fixedParameterIds.end());
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    std::vector<std::string> getObservableIds() const override {
        return std::vector<std::string>(observableIds.begin(),
                                        observableIds.end());
    }

    /**
     * @brief Get IDs of model expressions
     * @return Expression IDs
     */
    std::vector<std::string> getExpressionIds() const override {
        return std::vector<std::string>(expressionIds.begin(),
                                        expressionIds.end());
    }

    /**
     * @brief function indicating whether reinitialization of states depending
     * on fixed parameters is permissible
     * @return flag indicating whether reinitialization of states depending on
     * fixed parameters is permissible
     */
    bool isFixedParameterStateReinitializationAllowed() const override {
        return true;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    std::string getAmiciVersion() const override {
        return "0.31.0";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    std::string getAmiciCommit() const override {
        return "unknown";
    }

    bool hasQuadraticLLH() const override {
        return true;
    }

    ObservableScaling getObservableScaling(int iy) const override {
        return observableScalings.at(iy);
    }
};


} // namespace model_Oliveira_NatCommun2021

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
