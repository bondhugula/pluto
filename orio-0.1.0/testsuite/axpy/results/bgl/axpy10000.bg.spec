
spec align_unroll {
 def build   
 {  
   arg command = 'mpixlc ';
   arg options = '-O3 -qstrict -lm';  
 }  
  
 def performance_counter   
 {  
   arg method = 'basic timer';  
   arg repetitions = 10000; 
 }  
  
 def performance_params  
 {  
   param UF[] = range(1,20);
 }  
  
 def input_params  
 {  
   param SIZE = 10000;
   decl int n = SIZE; 
   decl double a1 = 4.232;
   decl double a2 = 134531.2145;
   decl double a3 = 43.24141;
   decl double a4 = 241.24314;
   decl double x1[n] = random;    
   decl double x2[n] = random;    
   decl double x3[n] = random;    
   decl double x4[n] = random;    
   decl double y[n] = 0;
 }  

 def search
 {
   arg algorithm = 'Exhaustive';
   arg time_limit = 20; 
   arg run_command = 'cqsub -n 64 -t 10 -q short ';
   arg num_processes = 64;
 }
}

