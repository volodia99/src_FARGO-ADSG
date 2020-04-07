/** \file Force.c

Contains the function used to evaluate the %force due to the disk, and
the function that writes the 'tqwk' log files. We calculate the force
on the planet by smoothly excluding a certain amount of the planet's
Hill Radius R_H. The exclusion function is equal to 0 below the
exclusion distance, equal to 1 beyond R_H, and rises as a sinus
squared in between. In absence of self-gravity, 11 exclusion distances
are considered (dimfxy is set by default to 11 in main.c), that is we
exclude 0/10, 1/10, 2/10,..., 10/10 times R_H. Each exclusion factor
is related to a unique file, e.g. tqwk0_4.dat for an exclusion of 4/10
R_H. Although the planet mass is given as an argument to
ComputeForce(), this mass is used only to specify the distance cutoff
in the case of the Hill sphere avoidance.  The force returned is a
specific force. It has therefore the dimension of an acceleration
(LT^-2).

October 2014: we no longer assign dimfxy to 11. It is now set to 2 to
save computing time (15 to 20%). Only the case 0 and 1 times R_H are
accounted for.
*/

#include "mp.h"

extern boolean OpenInner, NonReflecting, BarMass08;
extern Pair DiskOnPrimaryAcceleration;

Force *AllocateForce ()
{
  int i;
  Force *force;
  real *globalforce;
  force  = (Force *) prs_malloc (sizeof(Force));
  /* dimfxy is a global integer defined in global.h. It is set to 2
     (oct. 2014) */
  globalforce = (real *) prs_malloc (sizeof(real) * 4 * dimfxy);
  for (i = 0; i < 4*dimfxy; i++)
    globalforce[i] = 0.;
  force->GlobalForce = globalforce;
  return force;
}

void FreeForce (force)
     Force *force;
{
  free (force->GlobalForce);
}

void ComputeForce (force, Rho, x, y, rsmoothing, mass, sys, index)
     Force *force;
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
     PlanetarySystem *sys;
     int index;
{
  int i, j, l, ns, k;
  real xc, yc, cellmass, dx, dy, distance, d2, dist2, rh, a;
  int dimexclude;
  real x0, x1, y0, y1, m0, m1, xb, yb;
  real planet_distance, cutoff;
  real InvDist3, hill_cut, hillcutfactor;
  real *fxi, *fxo, *fyi, *fyo;
  real *localforce, *globalforce;
  real *dens, *abs, *ord, *azidens;
  real cutoffdist, raddistfrompla;
  /* dimfxy is a global integer defined in global.h. It is set to 2
     (oct. 2014) */
  fxi = (real *) prs_malloc (sizeof(real) * dimfxy);
  fxo = (real *) prs_malloc (sizeof(real) * dimfxy);
  fyi = (real *) prs_malloc (sizeof(real) * dimfxy);
  fyo = (real *) prs_malloc (sizeof(real) * dimfxy);
  localforce = (real *) prs_malloc (sizeof(real) * 4 * dimfxy);
  globalforce = force->GlobalForce;
  azidens = (real *) prs_malloc (sizeof(real) * GLOBALNRAD);
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  a = sqrt(x*x+y*y);  // star-planet distance
  rh = pow(mass/3., 1./3.)*a+1e-15;
  for ( k = 0; k < dimfxy; k++ ) {
    fxi[k] = 0.;
    fxo[k] = 0.;
    fyi[k] = 0.;
    fyo[k] = 0.;
  }
  for ( k = 0; k < 4*dimfxy; k++ ) {
    localforce[k] = 0.;
    globalforce[k] = 0.;
  }
  if (FakeSequential && (CPU_Rank > 0)) {
    MPI_Recv (&globalforce[0], 4*dimfxy, MPI_DOUBLE, CPU_Rank-1, 27, MPI_COMM_WORLD, &fargostat);
    for ( k = 0; k < dimfxy; k++ ) {
      fxi[k] = globalforce [k];
      fxo[k] = globalforce [k + dimfxy];
      fyi[k] = globalforce [k + 2*dimfxy];
      fyo[k] = globalforce [k + 3*dimfxy];
    }
  }
  if (sys->Binary[0] == YES) {
    /* Case of a binary-star system. Exclusion is done with
       respect to the barycenter of the two stars, the position of
       which is determined below */
    x0 = sys->x[0];
    x1 = sys->x[1];
    y0 = sys->y[0];
    y1 = sys->y[1];
    m0 = sys->mass[0];
    m1 = sys->mass[1];
    xb = (m0*x0 + m1*x1) / (m0+m1);
    yb = (m0*y0 + m1*y1) / (m0+m1);
  }

if(BarMass08) {
  for (i = Zero_or_active; i < Max_or_active; i++) {
    azidens[i]=0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      azidens[i] += dens[l];
    }
    azidens[i] /= (real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      dens[l] -= azidens[i];
    }
  }
}

#pragma omp parallel for private(j,hill_cut,cellmass,l,xc,yc,dist2,distance,InvDist3,dx,dy) shared(fxi,fyi,fxhi,fyhi,fxo,fyo,fxho,fyho)
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      xc = abs[l];
      yc = ord[l];
      cellmass = Surf[i]*dens[l];
      dx = xc-x;
      dy = yc-y;
      d2 = dx*dx+dy*dy;
      planet_distance = sqrt(d2);
      dist2 = d2 + rsmoothing*rsmoothing;
      distance = sqrt(dist2);
      InvDist3 = 1.0/dist2/distance;
      if (sys->Binary[0] == YES) {
	/* Remember if a binary-star is assumed, exclusion is
	   effective from the center-of-mass of the binary-star
	   system, instead of the location of each star:
	   planet_distance holds here the distance from a cell to the
	   barycenter of the two satellites */
	d2 = (xc-xb)*(xc-xb) + (yc-yb)*(yc-yb);
	planet_distance = sqrt(d2);
      }
      /* --------------- */
      for ( k = 0; k < dimfxy; k++ ) {
	if (dimfxy != 2)
	  hillcutfactor = (real)k / (real)(dimfxy-1);
	else
	  hillcutfactor = EXCLUDEHILLFACTOR;
	if ( k != 0 ) {
	  cutoff = hillcutfactor * rh;
	  /* New default exclusion function */
	  if (planet_distance/cutoff < 0.5)
	    hill_cut = 0.0;
	  else {
	    if (planet_distance > cutoff)
	      hill_cut = 1.0;
	    else
	      hill_cut = pow(sin((planet_distance/cutoff-.5)*M_PI),2.);
	  }
	  /* Old default exclusion function */
	  //hill_cut = 1.-exp(-d2/(cutoff*cutoff));
	}
	else
	  hill_cut = 1.; // if k=0
	if (Rmed[i] < a) {
#pragma omp atomic
	  fxi[k] += G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	  fyi[k] += G*cellmass*dy*InvDist3*hill_cut;
	} else {
#pragma omp atomic
	  fxo[k] += G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	  fyo[k] += G*cellmass*dy*InvDist3*hill_cut;
	}
      }
    }
  }

  if (BarMass08) {
    for (i = Zero_or_active; i < Max_or_active; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        dens[l] += azidens[i];
      }
    }
  }

  if (FakeSequential) {
    for ( k = 0; k < dimfxy; k++ ) {
      globalforce [k]            = fxi[k];
      globalforce [k + dimfxy]   = fxo[k];
      globalforce [k + 2*dimfxy] = fyi[k];
      globalforce [k + 3*dimfxy] = fyo[k];
    }
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&globalforce[0], 4*dimfxy, MPI_DOUBLE, CPU_Rank+1, 27, MPI_COMM_WORLD);
  } else {
    for ( k = 0; k < dimfxy; k++ ) {
      localforce [k]            = fxi[k];
      localforce [k + dimfxy]   = fxo[k];
      localforce [k + 2*dimfxy] = fyi[k];
      localforce [k + 3*dimfxy] = fyo[k];
    }
    MPI_Allreduce (localforce, globalforce, 4*dimfxy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential)
    MPI_Bcast (globalforce, 4*dimfxy, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
  if (dimfxy != 2)
    dimexclude = (int)(EXCLUDEHILLFACTOR*(dimfxy-1));
  else
    dimexclude = 1;
  /* New (july 2012): option with two planets to switch off outer
     force for inner planet, and inner force for outer planet */
  force->fx_inner    = globalforce[0];
  force->fx_ex_inner = globalforce[dimexclude];
  force->fx_outer    = globalforce[dimfxy];
  force->fx_ex_outer = globalforce[dimfxy+dimexclude];
  force->fy_inner    = globalforce[2*dimfxy];
  force->fy_ex_inner = globalforce[2*dimfxy+dimexclude];
  force->fy_outer    = globalforce[3*dimfxy];
  force->fy_ex_outer = globalforce[3*dimfxy+dimexclude];
  force->GlobalForce = globalforce;
  free (azidens);
  free (localforce);
  free (fxi);
  free (fxo);
  free (fyi);
  free (fyo);
}

real compute_smoothing (r)
     real r;
{
  real smooth;
  smooth = THICKNESSSMOOTHING * AspectRatio(r) * pow(r, 1.0+FLARINGINDEX);
  return smooth;
}


void UpdateLog (fc, psys, Rho, Energy, outputnb, time)
     Force *fc;
     PolarGrid *Rho, *Energy;
     PlanetarySystem *psys;
     int outputnb;
     real time;
{
  int i, nb;
  real x, y, r, m, vx, vy, smoothing;
  real *globalforce;
  FILE *out;
  char filename[MAX1D];
  nb = psys->nb;
  for (i = 0; i < nb; i++) {
    x = psys->x[i];
    y = psys->y[i];
    vx = psys->vx[i];
    vy = psys->vy[i];
    r = sqrt(x*x+y*y);
    m = psys->mass[i];
    if (RocheSmoothing)
      smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
    else
      smoothing = compute_smoothing(r);
    /* dimfxy is a global integer defined in global.h. It is set to 2
       (oct. 2014) */
    ComputeForce (fc, Rho, x, y, smoothing, m, psys, i);
    globalforce = fc->GlobalForce;
    if (CPU_Rank == CPU_Number-1) {
      sprintf (filename, "%stqwk%d.dat", OUTPUTDIR, i);
      out = fopen (filename, "a");
      if (out == NULL) {
	fprintf (stderr, "Can't open %s\n", filename);
	fprintf (stderr, "Aborted.\n");
      }
      fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb, \
	       x*fc->fy_inner-y*fc->fx_inner,				\
	       x*fc->fy_outer-y*fc->fx_outer,				\
	       x*fc->fy_ex_inner-y*fc->fx_ex_inner,			\
	       x*fc->fy_ex_outer-y*fc->fx_ex_outer,			\
	       vx*fc->fx_inner+vy*fc->fy_inner,				\
	       vx*fc->fx_outer+vy*fc->fy_outer,				\
	       vx*fc->fx_ex_inner+vy*fc->fy_ex_inner,			\
	       vx*fc->fx_ex_outer+vy*fc->fy_ex_outer, time);
      fclose (out);
    }
  }
}
