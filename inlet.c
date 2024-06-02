#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "navier-stokes/conserving.h" // not compatible with embed.h ??
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "tag.h"
#include "view.h"
#include "profiling.h"

// Viscosity ratio
double ratio_mu = 100. [0];

// Density ratio
double ratio_rho = 1000. [0];

// Reynolds number
double Re = 5. [0];

// Capillary number
double Ca = 0.070 [0];

// Drop initial length
double LI = 15. [0];    

// Channel length
double L1 = 40. [0];

// Refinement thresholds
double f_err = 1e-3, u_err = 1e-3;

// Refinement levels
int minlevel = 8, maxlevel = 11;

// Avoid division by zero
double EPS = 1e-3 [0];

// Channel geometry definition
#define INLET intersection((1. - y - EPS)*(1. + y - EPS), (1. - z - EPS)*(1. + z - EPS))

// Initial bubble geometry
#define CYLINDER_LEFT intersection((- sq(y/(0.9)) - sq(z/(0.9)) + sq(1.)), x - 2.)
#define CYLINDER intersection(-(x - LI - 1.), CYLINDER_LEFT)
#define ELLIPSOID_LEFT max(CYLINDER, (-sq((x - 2.)/(0.9)) - sq(y/(0.9)) - sq(z/(0.9)) + sq(1.)))
#define ELLIPSOID max(ELLIPSOID_LEFT, (-sq((x - LI - 1.)/(0.9)) - sq(y/(0.9)) - sq(z/(0.9)) + sq(1.)))

// Final and output times
#define TEND (L1)
#define TOUT (TEND/100.)
#define TDUMP (TEND/100.)
#define TMOVIE (TEND/400.)

// Boundary conditions
// Inflow
p[left] = neumann(0);
u.n[left] = dirichlet((1. - sq(y))*(1. - sq(z)));
u.t[left] = neumann(0);
u.r[left] = neumann(0);
f[left] = 0.;

// Outflow
p[right] = dirichlet(0);
u.n[right] = neumann(0);
u.t[right] = neumann(0);
u.r[right] = neumann(0);

// Solid no-slip
u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);
u.r[embed] = dirichlet(0);
f[embed] = 0.;

// Global scalar for curvature
scalar kappa[];

// Global double for the maximum position in x
double x_max_global = 0.;

// Output files
char dataname1[200], dataname2[200], moviename[200], snapname[100], facetname[100];

int main(int argc, char *argv[])
{
    if (argc > 1)
        Re = atof(argv[1]);
    if (argc > 2)
        Ca = atof(argv[2]);
    if (argc > 3)
        LI = atof(argv[3]);
    if (argc > 4)
        L1 = atof(argv[4]);
    if (argc > 5)
        maxlevel = atoi(argv[5]);
    if (argc > 6)
        minlevel = atoi(argv[6]);

    size(L1);
    origin(0., -L1/2., -L1/2.);
    stokes = true;

    init_grid(1 << minlevel);

    DT = 1. [0];

    rho1 = 1./ratio_rho, mu1 = (1./ratio_mu)/Re;
    rho2 = 1. [0], mu2 = 1./Re;

    f.sigma = 1./(Re*Ca);

    for (scalar s in {u})
        s.third = true;

    run();
}

event init(t = 0)
{
    sprintf(dataname1, "data_Re%.2f_Ca%.3f_LI%.1f_L%.1f_LEVEL%d", Re, Ca, LI, L1, maxlevel);
    sprintf(dataname2, "inlet_Re%.2f_Ca%.3f_LI%.1f_L%.1f_LEVEL%d", Re, Ca, LI, L1, maxlevel);
    sprintf(snapname, "snap_Re%.2f_Ca%.3f_LI%.1f_L%.1f_LEVEL%d", Re, Ca, LI, L1, maxlevel);
    sprintf(facetname, "facet_Re%.2f_Ca%.3f_LI%.1f_L%.1f_LEVEL%d", Re, Ca, LI, L1, maxlevel);

    if (!restore(file = "restart"))
    {
        if (maxlevel > minlevel)
        {
            scalar f0[], f1[];
            for (int lvl = minlevel; lvl <= maxlevel; lvl += 1)
            {
                fraction(f0, INLET);
                fraction(f1, ELLIPSOID);
                adapt_wavelet({f0, f1}, (double[]){0., 0.}, maxlevel = lvl);
            }
        }
        solid(cs, fs, INLET);
        fraction(f, ELLIPSOID);

        char name[200];
        sprintf(name, "properties_Re%.2f_Ca%.3f_LI%.1f_L%.1f_LEVEL%d", Re, Ca, LI, L1, maxlevel);
        FILE *fp = fopen(name, "w");

        fprintf(fp, "Grid = %g PPD = %g\n", L1/(pow(2, maxlevel)), 1./(L1/(pow(2, maxlevel))));
        fclose(fp);

        sprintf(moviename, "movie_Re%.2f_Ca%.3f_LI%.1f_L%.1f_LEVEL%d.mp4", Re, Ca, LI, L1, maxlevel);
    }
    else
    {
        if (maxlevel > minlevel)
        {
            scalar f0[];
            for (int lvl = minlevel; lvl <= maxlevel; lvl += 1)
            {
                fraction(f0, INLET);
                adapt_wavelet({f0, f, u}, (double[]){f_err, f_err, u_err, u_err, u_err}, maxlevel = lvl);
            }
        }
        solid(cs, fs, INLET);

        sprintf(moviename, "restart_movie_Re%.2f_Ca%.3f_LI%.1f_L%.1f_LEVEL%d.mp4", Re, Ca, LI, L1, maxlevel);
    }

}

event adapt(i++)
{
    adapt_wavelet({cs, f, u}, (double[]){f_err, f_err, u_err, u_err, u_err}, maxlevel);
}

event repulsive_force(i++){
    foreach(){
        if ((cs[] < 1.) && (f[] > 1e-3)){
            f[] = f[]*cs[]; // non-wetting condition
        }
    }
}

event output(i++)
{
    static FILE *fp1 = fopen(dataname1, "a");
    static FILE *fp2 = fopen(dataname2, "a");

    if (t == 0)
    {
        fprintf(fp1, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n", "#1 time", "#2 x_position", "#3 u_drop", "#4 Ca_drop", "#5 k_min", "#6 k_max", "#7 n_drops");
        fprintf(fp2, "%-12s %-12s %-12s\n\n", "#1 time", "#2 x_min", "#3 x_max");
    }

    scalar m[];
    foreach ()
        m[] = f[] > 1e-3;
    int n = tag(m);

    double uc = 0., xc = 0., vd = 0.;
    double x_min = 1e30, x_max = 0.;

    curvature(f, kappa);
    stats s = statsf(kappa);
    vector h[];
    heights (f, h);

    foreach (reduction(+:uc) reduction(+:xc) reduction(+:vd) reduction(min:x_min) reduction(max:x_max))
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
    }

    x_max_global = x_max;

    fprintf(fp1, "%-12g %-12g %-12g %-12g %-12g %-12g %-12d\n", t, xc/vd, uc/vd, (mu2/f.sigma)*uc/vd, s.min, s.max, n);
    fprintf(fp2, "%-12g %-12g %-12g\n", t, x_min, min(x_max, L1));

    fflush(fp1);
    fflush(fp2);
}

event output_interfaces(t = 0, t += TOUT){
    char name[400];
    sprintf(name, "interfaces_t%g_pid%d", t/TOUT, pid());
    FILE * fp3 = fopen (name, "w");

    foreach(){
        if (f[] > 1e-6 && f[] < 1. - 1e-6 && x < x_max_global - 11. && x > x_max_global - 11. - Delta) {
            face vector s = {{-1}};
            coord n = facet_normal (point, f, s);
            double alpha = plane_alpha (f[], n);
            coord v[12];
            int m = facets (n, alpha, v, 1.);
            for (int i = 0; i < m; i++)
	            fprintf (fp3, "%g %g %g\n", x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
            if (m > 0)
	            fputc ('\n', fp3);
        }
    }

    fclose(fp3);

    MPI_Barrier(MPI_COMM_WORLD);

    if (pid() == 0) {
        char command[1000];
        sprintf(command, "cat interfaces_t%g_pid* > %s_t%g", t/TOUT, facetname, t/TOUT);
        system(command);

        char command2[200];
        sprintf(command2, "rm interfaces_t%g_pid*", t/TOUT);
        system(command2);
    }
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

event snapshot(t = TDUMP, t += TDUMP, t <= TEND)
{
    char name[400];
    sprintf(name, "%s_t%g", snapname, t/TDUMP);
    dump(name);
}

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif