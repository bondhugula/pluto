let REGS=32; 

spec bratu {
 def build   
 {  
   arg command = '/Users/norris/software/openmpi-1.1.4/bin/mpicc ';  
   arg options = '-O3 -lm';  
 }  
  
 def performance_counter   
 {  
   arg method = 'basic timer';  
   arg repetitions = 10; 
 }  
  
 def performance_params  
 {  
   param Ui[] = range(1,16);
   param Uj[] = range(1,2);
   constraint reg_capacity = Ui * Uj <= REGS;
 }  
  
 def input_params  
 {  
   param SIZE = 512;
   param lambda = 6;
   decl int jl = 0;
   decl int il = 0;
   decl int jh = SIZE; 
   decl int ih = SIZE; 
   decl double x[ih][ih] = random;    
   decl double f[jh][jh] = 0;
 }  

 def search
 {
   arg algorithm = 'Exhaustive';
   arg time_limit = 10; 
   arg run_command = '/Users/norris/software/openmpi-1.1.4/bin/mpirun -np 4';
   arg num_processes = 4;
 }
}

