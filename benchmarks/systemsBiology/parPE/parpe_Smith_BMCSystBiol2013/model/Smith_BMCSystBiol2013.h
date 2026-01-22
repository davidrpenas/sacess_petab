#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_ode.h"
#include "amici/splinefunctions.h"

namespace amici {

class Solver;

namespace model_Smith_BMCSystBiol2013 {

extern std::array<const char*, 28> parameterNames;
extern std::array<const char*, 370> fixedParameterNames;
extern std::array<const char*, 109> stateNames;
extern std::array<const char*, 9> observableNames;
extern std::array<const ObservableScaling, 9> observableScalings;
extern std::array<const char*, 396> expressionNames;
extern std::array<const char*, 28> parameterIds;
extern std::array<const char*, 370> fixedParameterIds;
extern std::array<const char*, 109> stateIds;
extern std::array<const char*, 9> observableIds;
extern std::array<const char*, 396> expressionIds;
extern std::array<int, 109> stateIdxsSolver;
extern std::array<bool, 5> rootInitialValues;

extern void Jy_Smith_BMCSystBiol2013(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_Smith_BMCSystBiol2013(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_Smith_BMCSystBiol2013(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Smith_BMCSystBiol2013(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_Smith_BMCSystBiol2013(SUNMatrixWrapper &rowvals, int index);






extern void root_Smith_BMCSystBiol2013(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void dwdp_Smith_BMCSystBiol2013(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static);
extern void dwdp_colptrs_Smith_BMCSystBiol2013(SUNMatrixWrapper &colptrs);
extern void dwdp_rowvals_Smith_BMCSystBiol2013(SUNMatrixWrapper &rowvals);
extern void dwdx_Smith_BMCSystBiol2013(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static);
extern void dwdx_colptrs_Smith_BMCSystBiol2013(SUNMatrixWrapper &colptrs);
extern void dwdx_rowvals_Smith_BMCSystBiol2013(SUNMatrixWrapper &rowvals);
extern void dwdw_Smith_BMCSystBiol2013(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static);
extern void dwdw_colptrs_Smith_BMCSystBiol2013(SUNMatrixWrapper &colptrs);
extern void dwdw_rowvals_Smith_BMCSystBiol2013(SUNMatrixWrapper &rowvals);
extern void dxdotdw_Smith_BMCSystBiol2013(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Smith_BMCSystBiol2013(SUNMatrixWrapper &colptrs);
extern void dxdotdw_rowvals_Smith_BMCSystBiol2013(SUNMatrixWrapper &rowvals);






extern void dydx_Smith_BMCSystBiol2013(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void dydp_Smith_BMCSystBiol2013(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl);




extern void sigmay_Smith_BMCSystBiol2013(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);




extern void w_Smith_BMCSystBiol2013(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static);
extern void x0_Smith_BMCSystBiol2013(realtype *x0, const realtype t, const realtype *p, const realtype *k);

extern void sx0_Smith_BMCSystBiol2013(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);

extern void xdot_Smith_BMCSystBiol2013(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_Smith_BMCSystBiol2013(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);


extern void stau_Smith_BMCSystBiol2013(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie);
extern void deltax_Smith_BMCSystBiol2013(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old);
extern void deltasx_Smith_BMCSystBiol2013(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl);

extern void x_solver_Smith_BMCSystBiol2013(realtype *x_solver, const realtype *x_rdata);












extern std::vector<HermiteSpline> create_splines_Smith_BMCSystBiol2013(const realtype *p, const realtype *k);



/**
 * @brief AMICI-generated model subclass.
 */
class Model_Smith_BMCSystBiol2013 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Smith_BMCSystBiol2013()
        : amici::Model_ODE(
              amici::ModelDimensions(
                  109,                            // nx_rdata
                  109,                        // nxtrue_rdata
                  109,                           // nx_solver
                  109,                       // nxtrue_solver
                  0,                    // nx_solver_reinit
                  28,                                  // np
                  370,                                  // nk
                  9,                                  // ny
                  9,                              // nytrue
                  0,                                  // nz
                  0,                              // nztrue
                  5,                              // nevent
                  3,                       // nevent_solver
                  0,                                // nspl
                  1,                          // nobjective
                  396,                                  // nw
                  876,                               // ndwdx
                  17,                               // ndwdp
                  5,                               // ndwdw
                  685,                            // ndxdotdw
                  std::vector<int>{1,1,1,1,1,1,1,1,1},                              // ndjydy
                  0,                    // ndxrdatadxsolver
                  0,                        // ndxrdatadtcl
                  0,                        // ndtotal_cldx_rdata
                  0,                                       // nnz
                  109,                                 // ubw
                  109,                                 // lbw
                  true,                                    // pythonGenerated
                  0,                   // ndxdotdp_explicit
                  0,                   // ndxdotdx_explicit
                  1                    // w_recursion_depth
              ),
              amici::SimulationParameters(
                  std::vector<realtype>{1.2e-5, 90000.0, 2500.0, 1.9999999999999999e-6, 0.033333333999999999, 0.29999999999999999, 0.0020999999999999999, 0.00021000000000000001, 4.6099999999999999e-6, 1.0, 1.0, 8.2999999999999998e-12, 6.4000000000000005e-14, 1.6500000000000001e-11, 6.02e+23, 0.13145000000000001, 0.0054999999999999997, 2.7000000000000002, 0.0014, 2.5000000000000001e-5, 1.1878000000000001e-6, 3.4999999999999997e-5, 1.2500000000000001e-6, 0.00033300000000000002, 9.9999999999999995e-7, 0.014999999999999999, 7.5000000000000002e-6, 0.16700000000000001, 110.88, 0.001155, 1000.0, 200.0, 0.00027, 0.002, 2.0000000000000002e-5, 0.001, 0.25, 0.001, 0.12, 180.0, 740000000.0, 6.5000000000000003e-9, 2.0, 0.050000000000000003, 1000.0, 90000.0, 602000000.0, 0.00022000000000000001, 3.0, 100.0, 260.0, 0.001, 0.01, 0.87, 6.9500000000000002, 100.0, 100.0, 5.0, 2.0, 3.0, 15.0, 0.0, 0.0, 0.0054999999999999997, 0.181818181818, 0.055, 0.25, 0.125, 0.181818181818, 0.055, 0.25, 0.125, 1.81818181818, 0.0054999999999999997, 0.25, 0.125, 1.81818181818, 0.0054999999999999997, 0.25, 0.125, 0.090909090909100002, 0.55000000000000004, 0.25, 0.125, 0.090909090909100002, 0.55000000000000004, 0.25, 0.125, 0.90909090909099999, 0.055, 0.25, 0.125, 0.90909090909099999, 0.055, 0.25, 0.125, 0.018181818181800001, 0.55000000000000004, 0.125, 0.125, 0.018181818181800001, 0.55000000000000004, 0.125, 0.125, 0.181818181818, 0.055, 0.125, 0.125, 0.181818181818, 0.055, 0.125, 0.125, 0.0090909090909099995, 5.5, 0.125, 0.125, 0.0090909090909099995, 5.5, 0.125, 0.125, 0.090909090909100002, 0.55000000000000004, 0.125, 0.125, 0.090909090909100002, 0.55000000000000004, 0.125, 0.125, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 0.00029999999999999997, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 5.0000000000000002e-5, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 9.9999999999999995e-7, 2.1999999999999999e-5, 2.1999999999999999e-5, 2.1999999999999999e-5, 2.1999999999999999e-5, 2.1999999999999999e-5, 2.1999999999999999e-5, 3.0000000000000001e-6, 3.0000000000000001e-6, 3.0000000000000001e-6, 3.0000000000000001e-6, 3.0000000000000001e-6, 3.0000000000000001e-6, 6.6000000000000005e-5, 6.6000000000000005e-5, 6.6000000000000005e-5, 6.6000000000000005e-5, 6.6000000000000005e-5, 6.6000000000000005e-5, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 0.23999999999999999, 0.94999999999999996, 5.0, 0.22, 5.6219999999999999, 2.46, 0.0044000000000000003, 15.0, 0.22, 5.6219999999999999, 1.23, 0.0019}, // fixedParameters
                  std::vector<realtype>{2.0000000000000002e-5, 12072.0, 0.080000000000000002, 0.0050000000000000001, 450.0, 5.7999999999999998, 8.7499999999999999e-5, 2.7999999999999999e-6, 2.6000000000000001e-6, 1.55, 5.0000000000000002e-5, 4.9999999999999998e-7, 5.0000000000000002e-5, 4.9999999999999998e-7, 0.00059999999999999995, 0.00040000000000000002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}        // dynamic parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},   // idlist
              std::vector<int>{},               // z2events
              {{2880.0, {3}}, {2895.0, {4}}}               // state-independent events
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
        return new Model_Smith_BMCSystBiol2013(*this);
    }

    void fJrz(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_Smith_BMCSystBiol2013(Jy, iy, p, k, y, sigmay, my);
    }


    void fJz(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_Smith_BMCSystBiol2013(dJydsigma, iy, p, k, y, sigmay, my);
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
        deltasx_Smith_BMCSystBiol2013(deltasx, t, x, p, k, h, w, ip, ie, xdot, xdot_old, sx, stau, tcl);
    }


    void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {
        deltax_Smith_BMCSystBiol2013(deltax, t, x, p, k, h, ie, xdot, xdot_old);
    }


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


    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {}


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_Smith_BMCSystBiol2013(dJydy, iy, p, k, y, sigmay, my);
    }

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {
        dJydy_colptrs_Smith_BMCSystBiol2013(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {
        dJydy_rowvals_Smith_BMCSystBiol2013(rowvals, index);
    }


    std::vector<HermiteSpline> fcreate_splines(const realtype *p, const realtype *k) override {
        return create_splines_Smith_BMCSystBiol2013(p, k);
    }

    void fdspline_valuesdp(realtype *dspline_valuesdp, const realtype *p, const realtype *k, const int ip) override {}

    void fdspline_slopesdp(realtype *dspline_slopesdp, const realtype *p, const realtype *k, const int ip) override {}


    void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static) override {
        dwdp_Smith_BMCSystBiol2013(dwdp, t, x, p, k, h, w, tcl, dtcldp, spl, sspl, include_static);
    }

    void fdwdp_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdp_colptrs_Smith_BMCSystBiol2013(colptrs);
    }

    void fdwdp_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdp_rowvals_Smith_BMCSystBiol2013(rowvals);
    }


    void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static) override {
        dwdx_Smith_BMCSystBiol2013(dwdx, t, x, p, k, h, w, tcl, spl, include_static);
    }

    void fdwdx_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdx_colptrs_Smith_BMCSystBiol2013(colptrs);
    }

    void fdwdx_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdx_rowvals_Smith_BMCSystBiol2013(rowvals);
    }


    void fdwdw(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static) override {
        dwdw_Smith_BMCSystBiol2013(dwdw, t, x, p, k, h, w, tcl, include_static);
    }

    void fdwdw_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdw_colptrs_Smith_BMCSystBiol2013(colptrs);
    }

    void fdwdw_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdw_rowvals_Smith_BMCSystBiol2013(rowvals);
    }


    void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Smith_BMCSystBiol2013(dxdotdw, t, x, p, k, h, w);
    }

    void fdxdotdw_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdw_colptrs_Smith_BMCSystBiol2013(colptrs);
    }

    void fdxdotdw_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdw_rowvals_Smith_BMCSystBiol2013(rowvals);
    }


    void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        dydx_Smith_BMCSystBiol2013(dydx, t, x, p, k, h, w, dwdx);
    }


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl) override {
        dydp_Smith_BMCSystBiol2013(dydp, t, x, p, k, h, ip, w, tcl, dtcldp, spl, sspl);
    }


    void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        root_Smith_BMCSystBiol2013(root, t, x, p, k, h, tcl);
    }


    void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_Smith_BMCSystBiol2013(sigmay, t, p, k, y);
    }


    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) override {}


    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {
        stau_Smith_BMCSystBiol2013(stau, t, x, p, k, h, tcl, sx, ip, ie);
    }

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {
        sx0_Smith_BMCSystBiol2013(sx0, t, x, p, k, ip);
    }

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {}


    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static) override {
        w_Smith_BMCSystBiol2013(w, t, x, p, k, h, tcl, spl, include_static);
    }


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_Smith_BMCSystBiol2013(x0, t, p, k);
    }


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {}


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_Smith_BMCSystBiol2013(xdot, t, x, p, k, h, w);
    }


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_Smith_BMCSystBiol2013(y, t, x, p, k, h, w);
    }


    void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_Smith_BMCSystBiol2013(x_solver, x_rdata);
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
        return "Smith_BMCSystBiol2013";
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


} // namespace model_Smith_BMCSystBiol2013

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
