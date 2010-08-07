
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
  ) @*/  {
    for (j=jl; j<=jh-3; j=j+3) {
      for (i=il; i<=ih-8; i=i+8) {
        f[j][i]=(2*x[j][i]-x[j][i-1]-x[j][i+1])*hydhx+(2*x[j][i]-x[j-1][i]-x[j+1][i])*hxdhy-sc*exp(x[j][i]);
        f[(j+1)][i]=(2*x[(j+1)][i]-x[(j+1)][i-1]-x[(j+1)][i+1])*hydhx+(2*x[(j+1)][i]-x[j][i]-x[j+2][i])*hxdhy-sc*exp(x[(j+1)][i]);
        f[(j+2)][i]=(2*x[(j+2)][i]-x[(j+2)][i-1]-x[(j+2)][i+1])*hydhx+(2*x[(j+2)][i]-x[j+1][i]-x[j+3][i])*hxdhy-sc*exp(x[(j+2)][i]);
        f[j][(i+1)]=(2*x[j][(i+1)]-x[j][i]-x[j][i+2])*hydhx+(2*x[j][(i+1)]-x[j-1][(i+1)]-x[j+1][(i+1)])*hxdhy-sc*exp(x[j][(i+1)]);
        f[(j+1)][(i+1)]=(2*x[(j+1)][(i+1)]-x[(j+1)][i]-x[(j+1)][i+2])*hydhx+(2*x[(j+1)][(i+1)]-x[j][(i+1)]-x[j+2][(i+1)])*hxdhy-sc*exp(x[(j+1)][(i+1)]);
        f[(j+2)][(i+1)]=(2*x[(j+2)][(i+1)]-x[(j+2)][i]-x[(j+2)][i+2])*hydhx+(2*x[(j+2)][(i+1)]-x[j+1][(i+1)]-x[j+3][(i+1)])*hxdhy-sc*exp(x[(j+2)][(i+1)]);
        f[j][(i+2)]=(2*x[j][(i+2)]-x[j][i+1]-x[j][i+3])*hydhx+(2*x[j][(i+2)]-x[j-1][(i+2)]-x[j+1][(i+2)])*hxdhy-sc*exp(x[j][(i+2)]);
        f[(j+1)][(i+2)]=(2*x[(j+1)][(i+2)]-x[(j+1)][i+1]-x[(j+1)][i+3])*hydhx+(2*x[(j+1)][(i+2)]-x[j][(i+2)]-x[j+2][(i+2)])*hxdhy-sc*exp(x[(j+1)][(i+2)]);
        f[(j+2)][(i+2)]=(2*x[(j+2)][(i+2)]-x[(j+2)][i+1]-x[(j+2)][i+3])*hydhx+(2*x[(j+2)][(i+2)]-x[j+1][(i+2)]-x[j+3][(i+2)])*hxdhy-sc*exp(x[(j+2)][(i+2)]);
        f[j][(i+3)]=(2*x[j][(i+3)]-x[j][i+2]-x[j][i+4])*hydhx+(2*x[j][(i+3)]-x[j-1][(i+3)]-x[j+1][(i+3)])*hxdhy-sc*exp(x[j][(i+3)]);
        f[(j+1)][(i+3)]=(2*x[(j+1)][(i+3)]-x[(j+1)][i+2]-x[(j+1)][i+4])*hydhx+(2*x[(j+1)][(i+3)]-x[j][(i+3)]-x[j+2][(i+3)])*hxdhy-sc*exp(x[(j+1)][(i+3)]);
        f[(j+2)][(i+3)]=(2*x[(j+2)][(i+3)]-x[(j+2)][i+2]-x[(j+2)][i+4])*hydhx+(2*x[(j+2)][(i+3)]-x[j+1][(i+3)]-x[j+3][(i+3)])*hxdhy-sc*exp(x[(j+2)][(i+3)]);
        f[j][(i+4)]=(2*x[j][(i+4)]-x[j][i+3]-x[j][i+5])*hydhx+(2*x[j][(i+4)]-x[j-1][(i+4)]-x[j+1][(i+4)])*hxdhy-sc*exp(x[j][(i+4)]);
        f[(j+1)][(i+4)]=(2*x[(j+1)][(i+4)]-x[(j+1)][i+3]-x[(j+1)][i+5])*hydhx+(2*x[(j+1)][(i+4)]-x[j][(i+4)]-x[j+2][(i+4)])*hxdhy-sc*exp(x[(j+1)][(i+4)]);
        f[(j+2)][(i+4)]=(2*x[(j+2)][(i+4)]-x[(j+2)][i+3]-x[(j+2)][i+5])*hydhx+(2*x[(j+2)][(i+4)]-x[j+1][(i+4)]-x[j+3][(i+4)])*hxdhy-sc*exp(x[(j+2)][(i+4)]);
        f[j][(i+5)]=(2*x[j][(i+5)]-x[j][i+4]-x[j][i+6])*hydhx+(2*x[j][(i+5)]-x[j-1][(i+5)]-x[j+1][(i+5)])*hxdhy-sc*exp(x[j][(i+5)]);
        f[(j+1)][(i+5)]=(2*x[(j+1)][(i+5)]-x[(j+1)][i+4]-x[(j+1)][i+6])*hydhx+(2*x[(j+1)][(i+5)]-x[j][(i+5)]-x[j+2][(i+5)])*hxdhy-sc*exp(x[(j+1)][(i+5)]);
        f[(j+2)][(i+5)]=(2*x[(j+2)][(i+5)]-x[(j+2)][i+4]-x[(j+2)][i+6])*hydhx+(2*x[(j+2)][(i+5)]-x[j+1][(i+5)]-x[j+3][(i+5)])*hxdhy-sc*exp(x[(j+2)][(i+5)]);
        f[j][(i+6)]=(2*x[j][(i+6)]-x[j][i+5]-x[j][i+7])*hydhx+(2*x[j][(i+6)]-x[j-1][(i+6)]-x[j+1][(i+6)])*hxdhy-sc*exp(x[j][(i+6)]);
        f[(j+1)][(i+6)]=(2*x[(j+1)][(i+6)]-x[(j+1)][i+5]-x[(j+1)][i+7])*hydhx+(2*x[(j+1)][(i+6)]-x[j][(i+6)]-x[j+2][(i+6)])*hxdhy-sc*exp(x[(j+1)][(i+6)]);
        f[(j+2)][(i+6)]=(2*x[(j+2)][(i+6)]-x[(j+2)][i+5]-x[(j+2)][i+7])*hydhx+(2*x[(j+2)][(i+6)]-x[j+1][(i+6)]-x[j+3][(i+6)])*hxdhy-sc*exp(x[(j+2)][(i+6)]);
        f[j][(i+7)]=(2*x[j][(i+7)]-x[j][i+6]-x[j][i+8])*hydhx+(2*x[j][(i+7)]-x[j-1][(i+7)]-x[j+1][(i+7)])*hxdhy-sc*exp(x[j][(i+7)]);
        f[(j+1)][(i+7)]=(2*x[(j+1)][(i+7)]-x[(j+1)][i+6]-x[(j+1)][i+8])*hydhx+(2*x[(j+1)][(i+7)]-x[j][(i+7)]-x[j+2][(i+7)])*hxdhy-sc*exp(x[(j+1)][(i+7)]);
        f[(j+2)][(i+7)]=(2*x[(j+2)][(i+7)]-x[(j+2)][i+6]-x[(j+2)][i+8])*hydhx+(2*x[(j+2)][(i+7)]-x[j+1][(i+7)]-x[j+3][(i+7)])*hxdhy-sc*exp(x[(j+2)][(i+7)]);
      }
      for (; i<=ih-1; i=i+1) {
        f[j][i]=(2*x[j][i]-x[j][i-1]-x[j][i+1])*hydhx+(2*x[j][i]-x[j-1][i]-x[j+1][i])*hxdhy-sc*exp(x[j][i]);
        f[(j+1)][i]=(2*x[(j+1)][i]-x[(j+1)][i-1]-x[(j+1)][i+1])*hydhx+(2*x[(j+1)][i]-x[j][i]-x[j+2][i])*hxdhy-sc*exp(x[(j+1)][i]);
        f[(j+2)][i]=(2*x[(j+2)][i]-x[(j+2)][i-1]-x[(j+2)][i+1])*hydhx+(2*x[(j+2)][i]-x[j+1][i]-x[j+3][i])*hxdhy-sc*exp(x[(j+2)][i]);
      }
    }
    for (; j<=jh-1; j=j+1) {
      for (i=il; i<=ih-8; i=i+8) {
        f[j][i]=(2*x[j][i]-x[j][i-1]-x[j][i+1])*hydhx+(2*x[j][i]-x[j-1][i]-x[j+1][i])*hxdhy-sc*exp(x[j][i]);
        f[j][(i+1)]=(2*x[j][(i+1)]-x[j][i]-x[j][i+2])*hydhx+(2*x[j][(i+1)]-x[j-1][(i+1)]-x[j+1][(i+1)])*hxdhy-sc*exp(x[j][(i+1)]);
        f[j][(i+2)]=(2*x[j][(i+2)]-x[j][i+1]-x[j][i+3])*hydhx+(2*x[j][(i+2)]-x[j-1][(i+2)]-x[j+1][(i+2)])*hxdhy-sc*exp(x[j][(i+2)]);
        f[j][(i+3)]=(2*x[j][(i+3)]-x[j][i+2]-x[j][i+4])*hydhx+(2*x[j][(i+3)]-x[j-1][(i+3)]-x[j+1][(i+3)])*hxdhy-sc*exp(x[j][(i+3)]);
        f[j][(i+4)]=(2*x[j][(i+4)]-x[j][i+3]-x[j][i+5])*hydhx+(2*x[j][(i+4)]-x[j-1][(i+4)]-x[j+1][(i+4)])*hxdhy-sc*exp(x[j][(i+4)]);
        f[j][(i+5)]=(2*x[j][(i+5)]-x[j][i+4]-x[j][i+6])*hydhx+(2*x[j][(i+5)]-x[j-1][(i+5)]-x[j+1][(i+5)])*hxdhy-sc*exp(x[j][(i+5)]);
        f[j][(i+6)]=(2*x[j][(i+6)]-x[j][i+5]-x[j][i+7])*hydhx+(2*x[j][(i+6)]-x[j-1][(i+6)]-x[j+1][(i+6)])*hxdhy-sc*exp(x[j][(i+6)]);
        f[j][(i+7)]=(2*x[j][(i+7)]-x[j][i+6]-x[j][i+8])*hydhx+(2*x[j][(i+7)]-x[j-1][(i+7)]-x[j+1][(i+7)])*hxdhy-sc*exp(x[j][(i+7)]);
      }
      for (; i<=ih-1; i=i+1) {
        f[j][i]=(2.0*x[j][i]-x[j][i-1]-x[j][i+1])*hydhx+(2.0*x[j][i]-x[j-1][i]-x[j+1][i])*hxdhy-sc*exp(x[j][i]);
      }
    }
  }
/*@ end @*/


/*@ end @*/
