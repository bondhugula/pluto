
/*@ begin PerfTuning (
  import spec bratu;
) @*/

int i, j, k;

double hx     = 1.0/(ih-1);
double hy     = 1.0/(jh-1);
double sc     = hx*hy*lambda;
double hxdhy  = hx/hy;
double hydhx  = hy/hx;

/*@ begin Loop (
transform UnrollJam(ufactor=Uj)
for (j = jl; j <= jh-1; j++)
  transform UnrollJam(ufactor=Ui)
  for (i = il; i <= ih-1; i++) {
      f[j][i] = (2.0*x[j][i] - x[j][i-1] - x[j][i+1]) * hydhx + (2.0*x[j][i] - x[j-1][i] - x[j+1][i])*hxdhy - sc*exp(x[j][i]);
  }
) @*/

for (j=jl; j<jh; j++) {
    for (i=il; i<ih; i++) {
        f[j][i] = (2.0*x[j][i] - x[j][i-1] - x[j][i+1])*hydhx + (2.0*x[j][i] - x[j-1][i] - x[j+1][i])*hxdhy - sc*exp(x[j][i]);
    }
}
/*@ end @*/


/*@ end @*/
