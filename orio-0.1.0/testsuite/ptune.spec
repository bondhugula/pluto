
let REGS = 32;

spec unrolljam_mm_mul {
 def build 
 {
   arg command = 'gcc';
   arg options = '-O3';
 }

 def performance_counter 
 {
   arg method = 'basic timer';
   arg repetitions = 5;
 }

 def performance_params
 {
   param Ui[] = range(1,9);
   param Uj[] = range(1,9);
   param Uk[] = range(1,9);
   constraint reg_capacity = Ui*Uj + Ui*Uk + Uk*Uj <= REGS;
 }

 def input_params
 {
   param MVAL[] = [100,200];
   param NVAL[] = [100,200];
   param OVAL[] = [100,200];
   constraint rectangular_shape = (MVAL == NVAL == OVAL);
   decl in int M = MVAL;
   decl in int N = NVAL;
   decl in int O = OVAL;
   decl out static double X[MVAL][NVAL] = 0;
   decl in static double A[MVAL][OVAL] = random;
   decl in static double B[OVAL][NVAL] = random;
   #decl out dynamic double X[MVAL*NVAL] = 0;
   #decl in dynamic double A[MVAL*OVAL] = random;
   #decl in dynamic double B[OVAL*NVAL] = random;
 }

 def search
 {
   arg algorithm = 'Simplex';  
   arg time_limit = 0.1;
   #arg total_runs = 1;

   #arg simplex_local_distance = 1;
   #arg simplex_reflection_coef = 1.0;
   #arg simplex_expansion_coef = [1.5, 2.0];
   #arg simplex_contraction_coef = [0.5, 0.75];
   #arg simplex_shrinkage_coef = 0.7;
 }
}
