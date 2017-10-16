//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2017
//  Amanda Hart (altered from Gavin Fay/Steve Cadrin)
//
//  Lab 3: Age Structured Production Model
//  
//////////////////////////////////////////////
DATA_SECTION
  init_int fyear;   // first year of catch data, general so it can be changed by altering data file rather than line-by-line code
  init_int lyear;   // last year of catch data
  init_int fage; // Youngest age
  init_int lage; // Oldest age
  init_vector mat(fage,lage); // Proportion mature, define dimensions
  init_vector sel(fage,lage); // Selectivity
  init_vector wt(fage,lage); // Mean weight
  init_int ncpue;   // number of years of CPUE data
  init_matrix cat(fyear,lyear,1,2);   //matrix of catch
  init_ivector cpueyear(1,ncpue);   // years for cpue data points
  init_vector cpue(1,ncpue);   // cpue data   

  //number eps;

  int year;     // looping variables instead of for(i in 1:3) you write for(year in 1:3)
  int age;      // looping variable
  int icpue;    // looping variables

  // print some things to screen to check data read in correctly
  !!cout << nyear << endl; // number of years
  !!cout << ncpue << endl; // Length of CPUE series
  !!cout << cat(fyear,1) << " " << cat(lyear,2) << endl; // First and last year of catch
  !!cout << cpue(1) << endl; // first
  !!cout << cpue(ncpue) << endl;

PARAMETER_SECTION // Define all variables dimensions/type not identified in DATA SECTION
  // Transformed variables
  init_number dummy(-1);    // dummy parameter for debugging. -1 phase to turn off estimation

  init_number M(-1); // Rather than setting M as constant in data, M is added to Parameter section but not estimated(-1)
  init_number logit_h(2); // Logit of steepness(h), logit so between 0 and 1,
  // (2) above indicates it should be estimated in second phase (after other parameters estimated 1st then estimate steepness together with all the parameters)
  init_number logRzero(1); // log of asymptotic recruitment, 
  init_number logit_fracinit(1); // logit of the fraction of initial recruitment (don't assume Biomass in 1st yr =K and R 1st yr = Rzero)
  // init_number logq; // log of catchability
  init_number logsigma(1); // log of sigma

  // Variables back in meaningful units
  number h; // Steepness
  number Rzero; // Asymptotic recruitment
  number fracinit; // Initial year recruitment
  number q; // catchability
  number sigma; // observation error std. deviation
  number SSBzero;

  sdreport_vector biomass(fyear,lyear+1);    //model predictions of biomass

  // temporary variables
  number survival;  //survival from exploitation
  number fpen;  // Gradual penalty to prevent catch > biomass, if catch > biomass then penalty added so ADMB doesn't search that parameter space

  // Population dynamics variables
  vector Fful(fyear,lyear); // Fishing mortality for fully vulnerable ages (vulnerability determined by selectivity)
  matrix F(fyear,lyear,fage,lage); // Fishing mortality
  matrix N(fyear,lyear+1,fage,lage); // Abundance
  vector SSB(fyear,lyear+1); // Spawning biomass

  // Performance metrics
  sdreport_number BMSY;    //biomass at MSY
  sdreport_number MSY;     // maximum sustainable yield
  sdreport_number UMSY;    // exploitation rate at MSY

  objective_function_value obj_fun;    //value to be minimized


PROCEDURE_SECTION
  //dummy for testing
  obj_fun = square(dummy);

  //parameter transformations
  Rzero = mfexp(logRzero);
  //q = mfexp(logq);
  sigma = mfexp(logsigma);
  h = 0.2 + 0.8*(mfexp(logit_h)/(1.+mfexp(logit_h))); // logit_h can be 0-1, if =0 then this code adds 0.2 so that steepness reflect realistic negative density dependence (as pop increase recruitment levels off)
  fracinit = mfexp(logit_fracinit)/(1.+mfexp(logit_fracinit)); // retransform of logit


  //cout << logit_p << " " << p << endl;

  // initialize biomass in 1935
  // biomass(fyear) = p*K;

  // Initialize recruitment in 1935 as R0 (top left of matrix)
  N(fyear,fage) = Rzero*fracinit; // Rzero*fracinit since we are not assuming that Recruitment in first year = R0 (biomass doesn't = unfished)

  // loop for first year (first row of matrix) abundance to derive equilibrium abundance-at-age at F=0
  for (age=fage+1;age<=lage;age++) {
    N(fyear,age) = N(fyear,age-1)*mfexp(-1.*M);  // Assume no fishing mortality in the first year, survivors from previous age go to next year
  } 
   N(fyear,lage) = N(fyear,lage)/(1.-mfexp(-1.*M)); // Fill in last column of final age

  biomass(fyear) = sum(elem_prod(elem_prod(N(fyear)*mfexp(-0.5*M),sel),wt)); // Need to calculate fishing mortality (F = Yield/Biomass)
  // Biomass = sum((Survivor abundance(in fyear) * selectivity (-0.5*M))* Weight)
  // element products elem_prod() multiply age 10 survivor* age 10 selectivity* age 10 weight... for all ages
  SSB(fyear) = sum(elem_prod(elem_prod(N(fyear),mat),wt));
  // SSB (in fyear) = sum((Survivor abundance(in fyear) * maturity)*Weight)
  SSBzero = SSB(fyear)/fracinit; // Need starting SSB to get recruitment in year 2(stock recruit parameter), divide by fracinit since assuming initial recruitment isn't Rzero, biomass not at unfished
  // Similar to Bzero in Schaefer model

  // loop for series recruitment, older abundance, biomass, SSB, and F
  for (year=fyear+1;year<=lyear+1;year++){
    // Get F rate
    fpen = 0.;
    Ffull(year-1) = -1.*log(posfun((1.-cat(year-1,2)/biomass(year-1)),0.01,fpen));
    // Fishing mortality(year -1) 
    // = Survival (1-catch)/biomass
    obj_fun += 1000*fpen; // Big penalty if the value is negative

    // F-at-age 
    F(year-1) = Ffull(year-1)*sel; // F in year-1 * Selectivity gives F-at-age

    // Within each year loop over age
    for (age=fage+1;age<=lage;age++) {
      N(year,age) = N(year-1,age-1)*mfexp(-1.*(M+F(year-1,age-1)));
    }
    N(year,lage) += N(year-1,lage)*mfexp(-1.*(M+F(year-1,lage))); // Fill in bottom right corner, add remaining survival to this group

    biomass(ear) = sum(elem_prod(elem_prod(N(year)*mfexp(-0.5*M),sel),wt));

    // Recruitment
    N(year,fage) = 0.;
    SSB(year) = sum(elem_prod(elem_prod(N(year),mat,wt)));
    N(year,fage) = 4.*h*Rzero*SSB(year-1)/(SSBzero*(1.-h)+SSB(year-1)*(5.*h-1.)); // Recruitment depend on SSB in previous year
  }

  // MLE for catchability
    
///////////////////////////////////////////////// Everything past this may be wrong, finish typing/////////////////////

  // observation model
  // constant added to cpue predictions
  //eps = 1.e-07;

  //MLE for catchability
  q = 0.;
  for (icpue=1;icpue<=ncpue;icpue++) {
   q += log(cpue(icpue))-log(biomass(cpueyear(icpue)));
  }
  q = mfexp(q/ncpue);
  
  // compare model predictions to data
  for (icpue=1;icpue<=ncpue;icpue++) {
   obj_fun += logsigma + square(log(cpue(icpue))-log(q)-log(biomass(cpueyear(icpue))))/(2*square(sigma));
  }

  cout << obj_fun << endl;

  // calculate biological reference points
  BMSY = K/2;
  MSY = r*K/4;
  UMSY = r/2;
  
GLOBALS_SECTION
  #include <admodel.h>

REPORT_SECTION
//  report model estimates
  report << "biomass" << endl;
  for (year=fyear;year<=(lyear+1);year++)
   report << year << " " << biomass(year) << endl;
   
  report << "parameters" << endl;
  report << "r" << " " << r << endl;
  report << "K" << " " << K << endl;
  report << "p" << " " << p << endl;
  report << "sigma" << " " << sigma << endl;
  report << "q" << " " << q << endl;



