/** \file var.c

Contains the function that connects the string of the parameter file
to global variables.  The var() function is found in Interpret.c

*/

#define __LOCAL
#include "mp.h"
#undef __LOCAL

void
InitVariables()
{
  var("DT", &DT, REAL, YES, "1.");
  var("SIGMA0", &SIGMA0, REAL, YES, "173.");
  var("NINTERM", &NINTERM, INT, YES, "10.");
  var("NTOT", &NTOT, INT, YES, "1501.");
  var("OUTPUTDIR", OUTPUTDIR, STRING, YES, "~masset");
  var("INNERBOUNDARY", OPENINNERBOUNDARY, STRING, NO, "WALL");
  var("INNERBOUNDARYDUST", OPENINNERBOUNDARYDUST, STRING, NO, "WALL");
  var("LABELADVECTION", ADVLABEL, STRING, NO, "NO");
  var("TRANSPORT", TRANSPORT, STRING, NO, "FAST");
  var("PLANETCONFIG", PLANETCONFIG, STRING, NO, "Systems/SolarSystem.cfg");
  var("MASSTAPER", &MASSTAPER, REAL, NO, "0.0000001");
  var("RADIALSPACING", GRIDSPACING, STRING, NO, "ARITHMETIC");
  var("NRAD", &NRAD, INT, YES, "64.0");
  var("NSEC", &NSEC, INT, YES, "64.0");
  var("RMIN", &RMIN, REAL, YES, "1.0");
  var("RMAX", &RMAX, REAL, YES, "1.0");
  var("THICKNESSSMOOTHING", &THICKNESSSMOOTHING, REAL, NO, "0.3");
  var("SGTHICKNESSSMOOTHING", &SGTHICKNESSSMOOTHING, REAL, NO, "0.3");
  var("ROCHESMOOTHING", &ROCHESMOOTHING, REAL, NO, "0.0");
  var("ASPECTRATIO", &ASPECTRATIO, REAL, YES, "0.05");
  var("VISCOSITY", &VISCOSITY, REAL, NO, "0.0");
  var("RELEASEVISCOSITY", &RELEASEVISCOSITY, REAL, NO, "0.0");
  var("RELEASEDATEVISCOSITY", &RELEASEDATEVISCOSITY, REAL, NO, "0.0");
  var("ALPHAVISCOSITY", &ALPHAVISCOSITY, REAL, NO, "0.0");  
  var("SIGMASLOPE", &SIGMASLOPE, REAL, YES, "0.0");  
  var("RELEASERADIUS", &RELEASERADIUS, REAL, NO, "0.0");  
  var("RELEASEDATE", &RELEASEDATE, REAL, NO, "0.0");  
  var("OMEGAFRAME", &OMEGAFRAME, REAL, NO, "0.0");
  var("DISK", DISK, STRING, NO, "YES");
  var("FRAME", FRAME, STRING, NO, "FIXED");
  var("OUTERSOURCEMASS", OUTERSOURCEMASS, STRING, NO, "NO");
  var("WRITEDENSITY", WRITEDENSITY, STRING, NO, "YES");
  var("WRITEVELOCITY", WRITEVELOCITY, STRING, NO, "YES");
  var("WRITEENERGY", WRITEENERGY, STRING, NO, "NO");
  var("WRITETEMPERATURE", WRITETEMPERATURE, STRING, NO, "NO");
  var("WRITEDIVV", WRITEDIVV, STRING, NO, "NO");
  var("WRITEPOTENTIAL", WRITEPOTENTIAL, STRING, NO, "NO");
  var("WRITEVISCHEAT", WRITEVISCHEAT, STRING, NO, "NO");
  var("WRITETHERDIFF", WRITETHERDIFF, STRING, NO, "NO");
  var("WRITERADDIFF", WRITERADDIFF, STRING, NO, "NO");
  var("WRITETHERCOOL", WRITETHERCOOL, STRING, NO, "NO");
  var("WRITEDUSTDENSITY", WRITEDUSTDENSITY, STRING, NO, "NO");
  var("WRITEDUSTSTOP", WRITEDUSTSTOP, STRING, NO, "NO");
  var("WRITETEST", WRITETEST, STRING, NO, "NO");
  var("INDIRECTTERM", INDIRECTTERM, STRING, NO, "YES");
  var("EXCLUDEHILL", EXCLUDEHILL, STRING, NO, "NO");
  var("EXCLUDEHILLFACTOR", &EXCLUDEHILLFACTOR, REAL, NO, "1.0");
  var("IMPOSEDDISKDRIFT", &IMPOSEDDISKDRIFT, REAL, NO, "0.0");
  var("FLARINGINDEX", &FLARINGINDEX, REAL, NO, "0.0");
  var("ECCENTRICITY", &ECCENTRICITY, REAL, NO, "0.0");
  var("CAVITYRADIUS", &CAVITYRADIUS, REAL, NO, "0.0");
  var("CAVITYRATIO", &CAVITYRATIO, REAL, NO, "1.0");
  var("CAVITYWIDTH", &CAVITYWIDTH, REAL, NO, "1.0");
  var("TRANSITIONRADIUS", &TRANSITIONRADIUS, REAL, NO, "0.0");
  var("TRANSITIONRATIO", &TRANSITIONRATIO, REAL, NO, "1.0");
  var("TRANSITIONWIDTH", &TRANSITIONWIDTH, REAL, NO, "1.0");
  var("LAMBDADOUBLING", &LAMBDADOUBLING, REAL, NO, "0.0");
  var("SELFGRAVITY", SELFGRAVITY, STRING, NO, "NO");
  var("CICPLANET", CICPLANET, STRING, NO, "NO");
  var("FORCEDCIRCULAR", FORCEDCIRCULAR, STRING, NO, "NO");
  var("FORCEDINNERCIRCULAR", FORCEDINNERCIRCULAR, STRING, NO, "NO");
  var("ZMPLUS", ZMPLUS, STRING, NO, "NO");
  var("ENERGYEQUATION", ENERGYEQUATION, STRING, NO, "NO");
  var("ADIABATICINDEX", &ADIABATICINDEX, REAL, YES, "1.4");
  var("PLANETASPECTRATIO", &PLANETASPECTRATIO, REAL, NO, "0.5");
  var("ENTROPYDIFFUSION", ENTROPYDIFFUSION, STRING, NO, "NO");
  var("RADIATIVEDIFFUSION", RADIATIVEDIFFUSION, STRING, NO, "NO");
  var("THERMALCOOLING", THERMALCOOLING, STRING, NO, "NO");
  var("STELLARIRRADIATION", STELLARIRRADIATION, STRING, NO, "NO");
  var("STELLARLUMINOSITY", &STELLARLUMINOSITY, REAL, NO, "1.0");
  var("BACKGROUNDTEMPERATURE", &BACKGROUNDTEMPERATURE, REAL, NO, "50.0");
  var("SLOPEBACKGROUNDTEMPERATURE", &SLOPEBACKGROUNDTEMPERATURE, REAL, NO, "-0.5");
  var("VISCOUSHEATING", VISCOUSHEATING, STRING, NO, "YES");
  var("DIFFUSIVITY", &DIFFUSIVITY, REAL, NO, "0.000001");
  var("TEMPPRESC", TEMPPRESC, STRING, NO, "NO");
  var("PRESCTIME0", &PRESCTIME0, REAL, NO, "6.28");
  var("MHD", MHD, STRING, NO, "NO");
  var("SOFTWRITING", SOFTWRITING, STRING, NO, "NO");
  var("RETROGRADEPLANET", RETROGRADEPLANET, STRING, NO, "NO");
  var("GAMMATURB", &GAMMATURB, REAL, NO, "0.00001");
  var("LSAMODESPEEDUP", &LSAMODESPEEDUP, REAL, NO, "0.1");
  var("NBTURBMODES", &NBTURBMODES, INT, NO, "50");
  var("HIGHMCUTOFF", HIGHMCUTOFF, STRING, NO, "NO");
  var("PMIN", &PMIN, REAL, NO, "0.0");
  var("PMAX", &PMAX, REAL, NO, "6.2831853071795864");
  var("ADDMASS", ADDMASS, STRING, NO, "NO");
  var("BINARYSEPARATION", &BINARYSEPARATION, REAL, NO, "0.3");
  var("BINARYECCENTRICITY", &BINARYECCENTRICITY, REAL, NO, "0.0");
  var("RETROGRADEBINARY", RETROGRADEBINARY, STRING, NO, "NO");
  var("IMPOSEDDENSITY", IMPOSEDDENSITY, STRING, NO, "NO");
  var("FACTORUNITMASS", &FACTORUNITMASS, REAL, NO, "1.0");
  var("FACTORUNITLENGTH", &FACTORUNITLENGTH, REAL, NO, "1.0");
  var("FACTORMMW", &FACTORMMW, REAL, NO, "1.0");
  var("BETACOOLING", BETACOOLING, STRING, NO, "NO");
  var("BETACOOLINGTIME", &BETACOOLINGTIME, REAL, NO, "10.0");
  var("BETACOOLINGSLOPE", &BETACOOLINGSLOPE, REAL, NO, "0.0");
  var("WRITEGR", WRITEGR, STRING, NO, "NO");
  var("WRITEGTHETA", WRITEGTHETA, STRING, NO, "NO");
  var("ADDNOISE", ADDNOISE, STRING, NO, "NO");
  var("WKZRMIN", &WKZRMIN, REAL, NO, "0.0");
  var("WKZRMAX", &WKZRMAX, REAL, NO, "0.0");
  var("READPLANETFILEATRESTART", READPLANETFILEATRESTART, STRING, NO, "YES");
  var("DONTAPPLYSUBKEPLERIAN", DONTAPPLYSUBKEPLERIAN, STRING, NO, "NO");
  var("DENSDAMPRAD", &DENSDAMPRAD, REAL, NO, "0.0");
  var("DAMPTOINI", DAMPTOINI, STRING, NO, "NO");
  var("DAMPTOAXI", DAMPTOAXI, STRING, NO, "YES");
  var("DAMPTOVISCOUS", DAMPTOVISCOUS, STRING, NO, "NO");
  var("COROTATEWITHOUTERPLANET", COROTATEWITHOUTERPLANET, STRING, NO, "NO");
  var("DISCEVAPORATION", DISCEVAPORATION, STRING, NO, "NO");
  var("TEVAP", &TEVAP, REAL, NO, "0.0");
  var("CUSTIT", CUSTIT, STRING, NO, "NO");
  var("TURBRMIN", &TURBRMIN, REAL, NO, "0.0");
  var("TURBRMAX", &TURBRMAX, REAL, NO, "0.0");
  var("ADDFLOORS", ADDFLOORS, STRING, NO, "YES");      
  var("SIZEMINPART", &SIZEMINPART, REAL, NO, "0.001");
  var("SIZEMAXPART", &SIZEMAXPART, REAL, NO, "0.001");
  var("SIZEPARTSLOPE", &SIZEPARTSLOPE, REAL, NO, "3.0");
  var("RHOPART", &RHOPART, REAL, NO, "1.0");
  var("NBPART", &NBPART, INT, NO, "0");
  var("DUSTSLOPE", &DUSTSLOPE, REAL, NO, "0.5");
  var("DUSTFEELDISK", DUSTFEELDISK, STRING, NO, "YES");
  var("DUSTFEELSG", DUSTFEELSG, STRING, NO, "YES");
  var("DUSTFEELPLANETS", DUSTFEELPLANETS, STRING, NO, "YES");
  var("DUSTFEELTURB", DUSTFEELTURB, STRING, NO, "NO");
  var("DUSTFLUID", DUSTFLUID, STRING, NO, "NO");
  var("DUSTTIMESTEP", &DUSTTIMESTEP, REAL, NO, "0.0");
  var("MINDT", &MINDT, REAL, NO, "100.0");
  var("MAXDT", &MAXDT, REAL, NO, "0.00001");
  var("WRITEJACOBI", WRITEJACOBI, STRING, NO, "NO");
  var("WRITEDUSTSYSTEM", WRITEDUSTSYSTEM, STRING, NO, "YES");
  var("RMINDUST", &RMINDUST, REAL, NO, "0.0");
  var("RMAXDUST", &RMAXDUST, REAL, NO, "0.0");
  var("RESTARTWITHNEWDUST", RESTARTWITHNEWDUST, STRING, NO, "NO");
  var("ADDM1", ADDM1, STRING, NO, "NO"); 
  var("ADDM1TOM10", ADDM1TOM10, STRING, NO, "NO"); 
  var("TAILOFF", TAILOFF, STRING, NO, "NO");
  var("DENSITYJUMP", &DENSITYJUMP, REAL, NO, "1e-2");
  var("ZZINTEGRATOR", ZZINTEGRATOR, STRING, NO, "NO");
  var("NODTCONSTRAINTBYPCS", NODTCONSTRAINTBYPCS, STRING, NO, "YES");
  var("MDOTTIME", &MDOTTIME, REAL, NO, "1e5");
  var("PHOTOEVAPORATION", PHOTOEVAPORATION, STRING, NO, "NO");
  var("LX", &LX, REAL, NO, "1e30");
  var("DECINNER", DECINNER, STRING, NO, "NO");
  var("INTERPOLATION", INTERPOLATION, STRING, NO, "TSC");
  var("DUSTFEEDBACK", DUSTFEEDBACK, STRING, NO, "NO");
  var("SFTAPPROX", SFTAPPROX, STRING, NO, "YES");
  var("DUSTTOGASMASSRATIO", &DUSTTOGASMASSRATIO, REAL, NO, "1e-2");
  var("DUSTGROWTH", DUSTGROWTH, STRING, NO, "NO");
  var("REMOVEDUSTFROMPLANETSHILLRADIUS", REMOVEDUSTFROMPLANETSHILLRADIUS, STRING, NO, "YES");
  var("DUSTGROWTHPARAMETER", &DUSTGROWTHPARAMETER, REAL, NO, "0.0000001");
  var("DUSTMASSTAPER", &DUSTMASSTAPER, REAL, NO, "0.0000001");
  var("DASPECTRATIO", &DASPECTRATIO, REAL, YES, "0.05");
  var("DVISCOSITY", &DVISCOSITY, REAL, NO, "0.0");
  var("DALPHAVISCOSITY", &DALPHAVISCOSITY, REAL, NO, "0.0");
  var("DFLARINGINDEX", &DFLARINGINDEX, REAL, NO, "0.0");
  var("DUSTTOGASDENSITYRATIO", &DUSTTOGASDENSITYRATIO, REAL, NO, "0.01");
  var("SIZEPART", &Sizepart, REAL, NO, "1e-3");
  var("DUSTDIFFUSION", DUSTDIFFUSION, STRING, NO, "NO");
  var("DUSTDIFFINN", &DUSTDIFFINN, REAL, NO, "0.0");
  var("NELSONBOUND", &NELSONBOUND, INT, NO, "1");
  var("NELSONBOUNDD", &NELSONBOUNDD, INT, NO, "1");
  var("DENSITYFLOOR", &floordens, REAL, NO, "1e-9");
  var("DSIGMA0", &DSIGMA0, REAL, NO, "173.");
  var("DSIGMASLOPE", &DSIGMASLOPE, REAL, NO, "0.0");  
  var("NRAD1D", &NRAD1D, INT, NO, "100");
  var("RMIN1D", &RMIN1D, REAL, NO, "0.01");
  var("RMAX1D", &RMAX1D, REAL, NO, "100.0");
  var("BOUNDARY1DGRID", BOUNDARY1DGRID, STRING, NO, "T");
  var("COMPUTECPDMASS", COMPUTECPDMASS, STRING, NO, "NO");
  var("CUTDIST", &CUTDIST, REAL, NO, "100.0");
  var("FACTOROPACITIES", &FACTOROPACITIES, REAL, NO, "1.0");
  var("SETCONSTANTOPACITY", SETCONSTANTOPACITY, STRING, NO, "NO");
  var("IMPOSEDCONSTANTOPACITY", &IMPOSEDCONSTANTOPACITY, REAL, NO, "0.1");
}