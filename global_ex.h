/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varparser.pl" for details   */
/*                              */
/********************************/
extern int CPU_Rank;
extern int CPU_Number;
extern boolean CPU_Master;
extern int CPU_Next, CPU_Prev, CPU_Highest;
extern real unit_mass, unit_length, unit_temperature, unit_time, mmw, sigma_SB;
extern int DimToNext, DimToPrev, DimFromPrev, DimFromNext, intfoo;
extern int dimfxy;
extern int CPU_Friend, CPU_NoFriend;
extern real *dens_friend;
extern real *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
extern real *ffttohydro_transfer, *ffttohydro_transfer_friend;
extern int local_Nx, local_i_start, local_i_start_friend, total_local_size_friend, local_Nx_friend;
extern int local_Ny_after_transpose, local_j_start_after_transpose;
extern int total_local_size, ifront, Zero_or_active_friend;
extern int transfer_size, transfer_size_friend;
extern int hydro_totalsize, active_hydro_totalsize, active_hydro_totalsize_friend;  
extern int IMIN;
extern int IMAX;
extern int Zero_or_active;
extern int Max_or_active;
extern int One_or_active;
extern int MaxMO_or_active;		/* MO: Minus One */
extern int GLOBALNRAD;
extern int ievaporation;
extern real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
extern real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
extern real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D], Azimuth[MAX1D], CosAzimuth[MAX1D], SinAzimuth[MAX1D], AziInf[MAX1D], AziSup[MAX1D];
extern real SigmaMed[MAX1D], SigmaInf[MAX1D], MassTaper, DustMassTaper;
extern real DSigmaMed[MAX1D], DSigmaInf[MAX1D];
extern real EnergyMed[MAX1D], PrescTimeMed[MAX1D];
extern real InitialPlanetMass[MAX1D], FinalPlanetMass[MAX1D];
extern real VMed[MAX1D];
extern real GLOBAL_SoundSpeed[MAX1D], GLOBAL_DustSoundSpeed[MAX1D];
extern real OmegaFrame, PhysicalTime, PhysicalTimeInitial;
extern int TimeStep;
extern real HillRadius, mdcp, mdcp0, exces_mdcp;
extern real GLOBAL_bufarray[MAX1D];
extern real Particles_Mass, Particles_Mass_Initial;
extern boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
extern boolean	GotoNextOutput, StoreSigma, StoreEnergy, ViscosityAlpha, DViscosityAlpha, RocheSmoothing;
extern boolean DustDiffusion;
extern boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
extern MPI_Status fargostat;
extern PolarGrid *CellAbscissa, *CellOrdinate;
extern PolarGrid *RhoStar, *RhoInt, *Potential, *IndPotential, *TurbPotential, *Pressure, *SoundSpeed, *Temperature, *RadIndAcc, *AziIndAcc, *RadSGAcc, *AziSGAcc, *RadFBAcc, *AziFBAcc, *FBedot, *RadGradP, *AziGradP, *RadiativeKCoeff, *RadiativeChiCoeff;
extern PolarGrid *DRhoStar, *DRhoInt, *TStop, *Diag1, *Diag2, *Diag3, *Fdiffrp, *Fdifftp, *DPressure, *DSoundSpeed;
extern PolarGrid *DivergenceVelocity, *TAURR, *TAUPP, *TAURP, *ViscHeat, *EntropyDiff, *RadiativeDiff, *ThermCool, *Opacity;
extern PolarGrid *DDivergenceVelocity, *DTAURR, *DTAUPP, *DTAURP;
extern PolarGrid *gr, *gtheta;
extern PolarGrid *Test;
extern PolarGrid *a_SORarray, *b_SORarray, *c_SORarray, *d_SORarray, *e_SORarray, *f_SORarray;
extern PolarGrid *Global_a_SORarray, *Global_b_SORarray, *Global_c_SORarray, *Global_d_SORarray, *Global_e_SORarray, *Global_f_SORarray, *Global_tempint;
extern real *Radii1D, *Rmed1D, *Rsup1D, *Rinf1D, *Sigma1D, *Vrad1D, *cs1D, *viscosity1D, *term1, *term2, *term3;
extern boolean LogGrid;
extern boolean OverridesOutputdir;
extern char NewOutputdir[1024];
extern real *potturb;
extern real *GLOBAL_Axidens_Evap;
extern real VthetaMed[MAX1D], VradMed[MAX1D];
extern real DVthetaMed[MAX1D], DVradMed[MAX1D];
extern real Minimum_Stopping_Time;
