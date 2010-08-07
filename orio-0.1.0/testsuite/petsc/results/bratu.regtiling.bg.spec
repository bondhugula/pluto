let REGS=32; 

spec bratu {
 def build   
 {  
   arg command = 'mpixlc ';
   arg options = '-O3 -qstrict -lm';  
 }  
  
 def performance_counter   
 {  
   arg method = 'basic timer';  
   arg repetitions = 1000; 
 }  
  
 def performance_params  
 {  
   param Ui[] = range(1,32);
   param Uj[] = range(1,4);
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
   arg run_command = 'cqsub -n 64 -t 10 -q short ';
   arg num_processes = 64;
 }
}

