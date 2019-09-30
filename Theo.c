/** \file Theo.c

A few functions that manipulate the surface density, internal energy
and cooling time profiles.

*/

#include "mp.h"

extern real ScalingFactor; 

/* Surface density */
real Sigma(r)
     real r;
{
  real cavity = 1.0;
  real sigmabg;
  extern boolean TailOffGauss, TailOffIn, TailOffABA, TailOffOwen, ExponentialCutoff;
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  /* This is *not* a steady state */
  /* profile, if a cavity is defined. It first needs */
  /* to relax towards steady state, on a viscous time scale */
  sigmabg = cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
  if (TailOffGauss) {
    /* We take a Gaussian initial profile */
    sigmabg = SIGMA0*exp(-0.5*pow(r-1.0,2.0)*pow(2.0*ASPECTRATIO,-2.0)) + DENSITYJUMP*SIGMA0;
  }
  if (TailOffABA) {
    /* We take a Gaussian initial profile */
    if (r <= 1)
      sigmabg = SIGMA0*exp(-0.5*pow(r-1.0,2.0)*pow(0.1,-2.0)) + 1e-2*SIGMA0;
    else
      sigmabg = SIGMA0*exp(-0.5*pow(r-1.0,2.0)*pow(0.4,-2.0)) + 1e-2*SIGMA0;
  }
  if (TailOffIn) {
    /* We take here the prescription in Heemskrek, Papaloizou & Savonije 1992 */
    if (r < 1.4)
      sigmabg *= pow(r-GlobalRmed[0]+0.1,2.0);
  }
  if (TailOffOwen) {
    /* We take here the prescription in Owen et al. 2011b */
    sigmabg *= exp(-r/18.0);
  }
  if (ExponentialCutoff) {
    /* We take here the prescription in Owen et al. 2011b */
    sigmabg *= exp(-r*FACTORUNITLENGTH/CUTDIST); // CUTDIST is in AU in .par file
  }
  return sigmabg;
}

real DSigma(r)
real r;
{
  real cavity=1.0;
  real sigmabg;
  extern boolean TailOffGauss, RestartWithNewDust;
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  /* This is *not* a steady state */
  /* profile, if a cavity is defined. It first needs */
  /* to relax towards steady state, on a viscous time scale */
  sigmabg = cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE)*DUSTTOGASDENSITYRATIO;
  if (TailOffGauss) {
    /* We take a Gaussian initial profile */
    sigmabg = SIGMA0*DUSTTOGASDENSITYRATIO*exp(-0.5*pow(r-1.0,2.0)*pow(2.0*ASPECTRATIO,-2.0)) + DENSITYJUMP*DUSTTOGASDENSITYRATIO*SIGMA0;
  }
  if (RestartWithNewDust) {
    if ( (r >= RMINDUST) && (r <= RMAXDUST) )
      sigmabg = DSIGMA0*pow(r,-DSIGMASLOPE);
    else {
      if (r < RMINDUST)
	sigmabg =  DSIGMA0*pow(RMINDUST,-DSIGMASLOPE) * exp(-0.5*pow(r-RMINDUST,2.0)*pow(5.0*DAspectRatio(RMINDUST)*pow(RMINDUST,1.0+DFLARINGINDEX),-2.0));
      if (r > RMINDUST)
	sigmabg =  DSIGMA0*pow(RMAXDUST,-DSIGMASLOPE) * exp(-0.5*pow(r-RMAXDUST,2.0)*pow(5.0*DAspectRatio(RMAXDUST)*pow(RMAXDUST,1.0+DFLARINGINDEX),-2.0));
    }
    sigmabg += floordens;
  }
  return sigmabg;
}

void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void FillDSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    DSigmaMed[i] = DSigma(Rmed[i]);
    DSigmaInf[i] = DSigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

void RefillDSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    DSigmaMed[i] = moy;
  }
  DSigmaInf[0] = DSigmaMed[0];
  for (i = 1; i < nr; i++) {
    DSigmaInf[i] = (DSigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   DSigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

/* Thermal energy */
real Energy(r)
     real r;
{
  real energy0;
  real cavity = 1.0;
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  }
  else {
    energy0 = R/MU/(ADIABATICINDEX-1.0)*Sigma(r)*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX);
    //energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
  }
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  return cavity*ScalingFactor*energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    EnergyMed[i] = Energy(Rmed[i]);
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

/* Temperature prescription time */
real PrescTime(r)
     real r;
{
  real pt0;
  pt0 = PRESCTIME0*pow(r,2.0+2.0*FLARINGINDEX);
  return pt0;
}

void FillPrescTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    PrescTimeMed[i] = PrescTime(Rmed[i]);
}


/* Beta cooling prescription time */
real ComputeBetaCooling(r)
     real r;
{
  real pt0;
  pt0 = BETACOOLINGTIME*pow(r,-BETACOOLINGSLOPE);
  return pt0;
}
