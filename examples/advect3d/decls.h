#define nx 300
#define ny 300
#define nz 300
// #define T 20
// #define 9 9

#define f60 0.2 
#define f61 0.5
#define f62 0.3

#define halfdtbydx 0.5
#define thirddtbydz 0.3
#define thirddtbydx 0.3
#define thirddtbydy 0.3

double a[nx+10][ny+10][nz+10];
double af[nx+10][ny+10][nz+10];
double ab[nx+10][ny+10][nz+10];
double al[nx+10][ny+10][nz+10];
double athird[nx+10][ny+10][nz+10];
double uxl[nx+10][ny+10][nz+10];
double uzf[nx+10][ny+10][nz+10];
double uyb[nx+10][ny+10][nz+10];
