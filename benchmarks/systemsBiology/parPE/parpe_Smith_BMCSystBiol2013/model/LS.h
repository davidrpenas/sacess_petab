#ifndef LS_H
#define LS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Call ParPE optimization application.
 * 
 * @param name_input Input HDF5 file path.
 * @param name_output Output directory or file path.
 * @param id Local ID for parallel execution or logging.
 * @param sol_vect Pointer to the initial solution vector.
 * @param fx Pointer to the variable where the objective function value will be stored.
 * @param N Size of the solution vector.
 */
int call_parPE(int id, double *fx, int N, int *nevals, int iters, double init_time, double max_time, double *sol_vect_ptr, const char *name_input, const char *name_output);

#ifdef __cplusplus
}
#endif

#endif // LS_H
