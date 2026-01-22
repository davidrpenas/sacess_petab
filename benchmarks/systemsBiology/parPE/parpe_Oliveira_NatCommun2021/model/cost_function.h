#ifdef __cplusplus
extern "C" {
#endif

extern void load_proplem(const char* input_string);
extern void return_bounds(double* ub, double* lb);
extern int return_size_problem();
extern void eval_cost(double* solution, double* fx, double **res, int *size_res, double* grad);


#ifdef __cplusplus
}
#endif
