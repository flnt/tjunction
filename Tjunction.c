// Velocity based inflow
// Non-wetting boundary condition
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "tag.h"
#include "view.h"
#include "profiling.h"

// Case B by default
char *str = "A";
double wi = 30e-6 [1], wo = 14e-6 [1], hc = 37e-6 [1];

// Drop length to width ratio
double l0 = 8. [0];    

// Capillary number based on inflow velocity
double Ca = 0.01 [0];

// Repulsive force coefficient (zero by default)
double C_f = 0. [0];

// Refinement thresholds
double f_err = 1e-3, u_err = 1e-3;

// Refinement levels
int minlevel = 7, maxlevel = 7;

// Avoid division by zero
double EPS = 1e-3 [0];

// Physical parameters
double rho_1 = 1e3 [1, -3, 0]; // Water density
double mu_1 = 1e-3 [1, -1, -1]; // Water viscosity
double rho_2 = 1e3 [1, -3, 0]; // Fluorinated oil density
double mu_2 = 2.3e-3 [1, -1, -1]; // Fluorinated oil viscosity
double sig = 1.6e-3 [1, 0, -2]; // Water to oil surface tension

#define u0 (Ca*sig/mu_2) // Scaling velocity in m/s
#define Re (rho_2*u0*wo/mu_2) // Reynolds number
#define We (Re*Ca) // Webber number
#define T_unit (wo/u0)

// T-junction geometry and initial drop shape
#define WI (wi/wo) // Dimensionless inlet width
#define WO (wo/wo) // Dimensionless outlet width
#define HC (hc/wo) // Dimensionless height

#define LD (l0*WI) // Dimensionless initial drop length
#define L1 (LD + WI) // Dimensionless inlet length
#define L0 (L1*WI/WO) // Dimensionless outlet length

#define INLET intersection((WI/2. - y - EPS)*(WI/2. + y - EPS), (HC/2. - z - EPS)*(HC/2. + z - EPS))
#define OUTLET intersection(((x - L1 - WO) - EPS)*(-(x - L1) - EPS), (HC/2. - z - EPS)*(HC/2. + z - EPS))
#define TJUNC max(intersection(INLET, (-(x - L1) - EPS)), OUTLET)

#define CYLINDER_LEFT intersection((- sq(y/(0.95*WI/2.)) - sq(z/(0.95*HC/2.)) + sq(1.)), x - WI)
#define CYLINDER intersection(-(x - LD), CYLINDER_LEFT)
#define ELLIPSOID_LEFT max(CYLINDER, (-sq((x - WI)/(WI/2.)) - sq(y/(0.95*WI/2.)) - sq(z/(0.95*HC/2.)) + sq(1.)))
#define ELLIPSOID max(ELLIPSOID_LEFT, (-sq((x - LD)/(WI/2.)) - sq(y/(0.95*WI/2.)) - sq(z/(0.95*HC/2.)) + sq(1.)))

// #define ELLIPSOID (-sq((x - LI)/(LD/2.)) - sq(y/(0.95*WI/2.)) - sq(z/(0.95*HC/2.)) + sq(1.))

// Final and output times
#define TEND (L1 + L0/2.)
#define TOUT (TEND/100.)
#define TMOVIE (TEND/4000.)

// Function for geometrical factor of flow rate in rectangular cross section
double beta_function(double a, double b) {
    double beta = 0., sum = 0., term = 0.;
    for (int n = 1; n < 60; n += 2) {
        term = tanh(n*pi*b/(2.*a))/pow(n, 5);
        sum += term;
        if (fabs(term) < 1e-15)
            break;
    }
    beta = (4.*pow(a, 3)/3.)*(b - (192.*a/(pow(pi, 5))*sum));
    return beta;
}

// Boundary conditions
// Inlet
double P_in = 0.; // Inlet pressure
p[left] = neumann(0);
u.n[left] = dirichlet(1.); 
u.t[left] = neumann(0);
u.r[left] = neumann(0);
f[left] = 0.;

// Outlet = outflow
p[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[top] = neumann(0);
u.r[top] = neumann(0);
p[bottom] = dirichlet(0);
u.n[bottom] = neumann(0);
u.t[bottom] = neumann(0);
u.r[bottom] = neumann(0);

// Solid no-slip + repulsive force
vector repf[];
u.n[embed] = dirichlet(repf.x[]);
u.t[embed] = dirichlet(repf.y[]);
u.r[embed] = dirichlet(repf.z[]);
f[embed] = 0.;

// Global scalar for curvature
scalar kappa[];

// Output files
char dataname1[200], dataname2[200], dataname3[200], moviename[200], snapname[100];

int main(int argc, char *argv[])
{
    if (argc > 1)
        str = argv[1];
    if (argc > 2)
        l0 = atof(argv[2]);
    if (argc > 3)
        Ca = atof(argv[3]);
    if (argc > 4)
        C_f = atof(argv[4]);
    if (argc > 5)
        maxlevel = atoi(argv[5]);
    if (argc > 6)
        minlevel = atoi(argv[6]);

    if (strcmp(str, "A") == 0)
    {
        wi = 30e-6 [1], wo = 14e-6 [1], hc = 62e-6 [1];
    }
    else if (strcmp(str, "B") == 0)
    {
        wi = 30e-6 [1], wo = 14e-6 [1], hc = 37e-6 [1];
    }
    else if (strcmp(str, "C") == 0)
    {
        wi = 30e-6 [1], wo = 14e-6 [1], hc = 5e-6 [1];
    }
    else if (strcmp(str, "D") == 0)
    {
        wi = 60e-6 [1], wo = 14e-6 [1], hc = 37e-6 [1];
    }
    else if (strcmp(str, "E") == 0)
    {
        wi = 45e-6 [1], wo = 14e-6 [1], hc = 37e-6 [1];
    }
    else if (strcmp(str, "F") == 0)
    {
        wi = 30e-6 [1], wo = 30e-6 [1], hc = 80e-6 [1];
    }
    else if (strcmp(str, "G") == 0)
    {
        wi = 100e-6 [1], wo = 12e-6 [1], hc = 37e-6 [1];
    }
    else if (strcmp(str, "H") == 0)
    {
        wi = 30e-6 [1], wo = 12e-6 [1], hc = 85e-6 [1];
    }

    size(L0);
    origin(0., -L0/2., -L0/2.);
    stokes = true;

    init_grid(1 << minlevel);

    DT = 1. [0];

    rho1 = (rho_1/rho_2), mu1 = (mu_1/mu_2)/Re;
    rho2 = 1. [0], mu2 = 1./Re;

    f.sigma = 1./We;

    P_in = WI*HC*((L1)/(Re*beta_function(WI/2.,HC/2.)) + (L0/2.)/(Re*beta_function(WO,HC/2.)));

    for (scalar s in {u})
        s.third = true;

    run();
}

event init(t = 0)
{
    sprintf(dataname1, "data_case%s_l%g_Ca%.3f_f%g_LEVEL%d", str, l0, Ca, C_f, maxlevel);
    sprintf(dataname2, "inlet_case%s_l%g_Ca%.3f_f%g_LEVEL%d", str, l0, Ca, C_f, maxlevel);
    sprintf(dataname3, "outlet_case%s_l%g_Ca%.3f_f%g_LEVEL%d", str, l0, Ca, C_f, maxlevel);
    sprintf(snapname, "snap_case%s_l%g_Ca%.3f_f%g_LEVEL%d", str, l0, Ca, C_f, maxlevel);

    if (!restore(file = "restart"))
    {
        if (maxlevel > minlevel)
        {
            scalar f0[], f1[];
            for (int lvl = minlevel; lvl <= maxlevel; lvl += 1)
            {
                fraction(f0, TJUNC);
                fraction(f1, ELLIPSOID);
                adapt_wavelet({f0, f1}, (double[]){0., 0.}, maxlevel = lvl);
            }
        }
        solid(cs, fs, TJUNC);
        fraction(f, ELLIPSOID);

        char name[200];
        sprintf(name, "properties_case%s_l%g_Ca%.3f_f%g_LEVEL%d", str, l0, Ca, C_f, maxlevel);
        FILE *fp = fopen(name, "w");

        fprintf(fp, "wi (m) = %g wo (m) = %g hc (m) = %g\n", wi, wo, hc);
        fprintf(fp, "Re = %g Ca = %g We = %g P_in = %g\n", Re, Ca, We, P_in);
        fprintf(fp, "Grid = %g Grid (m) = %g PPD = %g\n", L0/(pow(2, maxlevel)), wo*L0/(pow(2, maxlevel)), min(WO, HC)/(L0/(pow(2, maxlevel))));
        fprintf(fp, "T_unit (s) = %g u0 (m/s) = %g\n", T_unit, u0);
        fprintf(fp, "tend (s) = %g tout (s) = %g\n", TEND*T_unit, TOUT*T_unit);
        fprintf(fp, "TEND = %g TOUT = %g\n", TEND, TOUT);

        fclose(fp);

        sprintf(moviename, "movie_case%s_l%g_Ca%.3f_f%g_LEVEL%d.mp4", str, l0, Ca, C_f, maxlevel);
    }
    else
    {
        if (maxlevel > minlevel)
        {
            scalar f0[];
            for (int lvl = minlevel; lvl <= maxlevel; lvl += 1)
            {
                fraction(f0, TJUNC);
                adapt_wavelet({f0, f, u}, (double[]){f_err, f_err, u_err, u_err, u_err}, maxlevel = lvl);
            }
        }
        solid(cs, fs, TJUNC);

        sprintf(moviename, "restart_movie_case%s_l%g_Ca%.3f_f%g_LEVEL%d.mp4", str, l0, Ca, C_f, maxlevel);
    }

}

event adapt(i++)
{
    adapt_wavelet({cs, f, u}, (double[]){f_err, f_err, u_err, u_err, u_err}, maxlevel);
}

event repulsive_force(i++){
    foreach(){
        foreach_dimension(){
            repf.x[] = 0.;
        }
        if ((cs[] < 1.) && (f[] > 1e-3)){
            coord n = interface_normal (point, cs);
            foreach_dimension(){
                repf.x[] = -C_f*n.x*f[];
            }
            f[] = f[]*cs[]; // non-wetting condition
        }
    }
}

event output(i++)
{
    static FILE *fp1 = fopen(dataname1, "a");
    static FILE *fp2 = fopen(dataname2, "a");
    static FILE *fp3 = fopen(dataname3, "a");

    if (t == 0)
    {
        fprintf(fp1, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n", "#1 time", "#2 time (s)", "#3 x_position", "#4 Ca_drop", "#5 Ca", "#6 k_min", "#7 k_max", "#8 n_drops");
        fprintf(fp2, "%-12s %-12s %-12s %-12s %-12s %-12s\n\n", "#1 time", "#2 time (s)", "#3 x_min", "#4 x_max", "#5 x_min (m)", "#6 x_max (m)");
        fprintf(fp3, "%-12s %-12s %-12s %-12s %-12s %-12s\n\n", "#1 time", "#2 time (s)", "#3 y_min", "#4 y_max", "#5 y_min (m)", "#6 y_max (m)");
    }

    scalar m[];
    foreach ()
        m[] = f[] > 1e-3;
    int n = tag(m);

    double uc = 0., xc = 0., vd = 0.;
    double x_min = 1e30, x_max = 0., y_min = 0., y_max = 0.;

    curvature(f, kappa);
    stats s = statsf(kappa);
    vector h[];
    heights (f, h);

    foreach (reduction(+:uc) reduction(+:xc) reduction(+:vd) reduction(min:x_min) reduction(max:x_max) reduction(min:y_min) reduction(max:y_max))
    {
        double dv1 = f[]*dv();
        xc += x*dv1;
        uc += u.x[]*dv1;
        vd += dv1;
        if ((h.x[] != nodata) & (f[] < (1. - 1e-6)) & (f[] > 1e-6) & (fabs(y) < 2.*Delta) & (fabs(z) < 2.*Delta))
        {
            double pos_x = x +  Delta*height(h.x[]);
            if (pos_x > x_max){
                x_max = pos_x;
            }
            if (pos_x < x_min){
                x_min = pos_x;
            }
        }
        if ((h.y[] != nodata) & (f[] < (1. - 1e-6)) & (f[] > 1e-6) & (x > (L1 + WO/2. - 2.*Delta)) & (x < (L1 + WO/2. + 2.*Delta)) & (fabs(z) < 2.*Delta))
        {
            double pos_y = y +  Delta*height(h.y[]);
            if (pos_y > y_max){
                y_max = pos_y;
            }
            if (pos_y < y_min){
                y_min = pos_y;
            }
        }
    }


    fprintf(fp1, "%-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12d\n", t, t*T_unit, xc/vd, (mu2/f.sigma)*uc/vd, (mu2/f.sigma), s.min, s.max, n);
    fprintf(fp2, "%-12g %-12g %-12g %-12g %-12g %-12g\n", t, t*T_unit, x_min, min(x_max, L1), x_min*wo, min(x_max, L1)*wo);
    fprintf(fp3, "%-12g %-12g %-12g %-12g %-12g %-12g\n", t, t*T_unit, y_min, y_max, y_min*wo, y_max*wo);

    fflush(fp1);
    fflush(fp2);
    fflush(fp3);
}

event movie(t = TMOVIE, t += TMOVIE)
{
    scalar pid[];
    foreach ()
        pid[] = fmod(pid()*(npe() + 37), npe());
    view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.436, ty = 0.087, tz = -3.882);
    box ();
    draw_vof (c = "f");
    squares (color = "cs", alpha = -0.001);
    save(moviename);
}

event snapshot(t = TOUT, t += TOUT, t <= TEND)
{
    char name[400];
    sprintf(name, "%s_t%g", snapname, t/TOUT);
    dump(name);
}

// #if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
// #endif