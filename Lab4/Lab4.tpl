//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2017
//  Steve Cadrin & Gavin Fay
//
//  Age Structured Production Model lab
//  time series observation error model assuming Beverton-Holt stock-recruitment
//////////////////////////////////////////////
DATA_SECTION
  init_int diag; // print diagnostics?
  init_int fyear;   // first year of catch data
  init_int lyear;   // last year of catch data
  init_int fage; // youngest age
  init_int lage; // oldest age
  init_vector mat(fage,lage); // proportion mature
  init_vector sel(fage,lage); // selectivity
  init_vector wt(fage,lage); // mean weight
  init_int frecdev; // first year recruitment deviation
  init_int  lrecdev; // last year recruitment deviation
  int nrecdevs; // number of years of recruitment deviation
  !! nrecdevs = lrecdev-frecdev+1; 
  init_int ncpue;   // number of years of CPUE data
  init_matrix cat(fyear,lyear,1,2);   //matrix of catch
  init_ivector cpueyear(1,ncpue);   // years for cpue data points
  init_vector cpue(1,ncpue);   // cpue data
  init_int nsurvey;
  init_matrix survey(1,nsurvey,1,2); // matrix of survey data
     // Index survey by number of observations in survey data set (don't need a data point for each year)
  init_int nage_survey;
  int nage;
   !! nage = lage-fage+1;
   // Or use syntax in next three lines
   // LOCAL_CALCS
   // nage = lage-fage+1;
   // END_CALCS
   init_matrix agecomp_s(1,nage_survey,1,nage+2); // Survey age composition, matrix(1, nrow, ncol)
   init_int nage_fishery;
   init_matrix agecomp_f(1,nage_fishery,1,nage+2); // Fishery age composition

   ivector agecols(1,nage);
   ivector ages(1,nage);
   !! for (age=1;age<=nage;age++) agecols(age) = 2+ age; // makes index so it is easy to pull out columns with corresponding age (since first two columns are age and number samples)
   !! for (age=1;age<=nage;age++) ages(age) = fage + age -1;

  int year;     // looping variable
  int age;    // looping variable
  int icpue;    // looping variable

  number eps; // Small constant to add to stuff (arbitrary name)
  !! eps = 1.e-07; // Exclamation points indicate this is C++ code, don't search data file for it

  // Check that we read things in correctly
  !! cout << agecomp_f(nage_fishery) << endl;
  //number M;     // natural mortality
  // print some things to screen to check data read in correctly
//  cout << ncpue << endl;
//  cout << cat(fyear,1) << " " << cat(lyear,2) << endl;
//  cout << cpue(1) << endl;
//  cout << cpue(ncpue) << endl;

PARAMETER_SECTION
  init_number dummy(-1);    // dummy parameter for debugging. -1 phase to turn off estimation
  // -1 phase turned off
  init_number M(-1);
  init_number logit_h(-3);  //logit of steepness
  init_number logRzero(-1);  //log of asymptotic recruitment
  init_number logit_fracinit(-1);

  // init_number ldgq; // log of catchability
  // init_number logsigma(1); // log of sigma
  init_number log_fA50(-1);  //log of age at 50% selectivity for fishery
  init_number log_fbeta(-1); // log fo slope in selectivity for fishery
  init_number log_sA50(-1); // log of age at 50% selectivity for survey
  init_number log_sbeta(-1); //log of slope in selectivity for survey
  init_number log_sigmaR(-1); //log of standard deviation of recruitment deviations
  init_bounded_vector recdevs(frecdev,lrecdev,-10,10,-2); // recruitment deviations, bounded
     // recdevs (first year, last year, )

  vector f_sel(fage,lage); // fishery selectivity
  vector s_sel(fage,lage); // survey selectivity

  number h;    //steepness
  number Rzero;    //asymptotic recruitment
  number fracinit; //initial yr recruitment
  number s_q; // Survey catchability
  number s_sigma; // survey observation error std. deviation
  // number q; //test selectivity
  // number sigma; // test fishery catchability
  // number logsima; //test
  number f_q; // fishery ctchability
  number f_sigma; // fishery observation error std. deviation
  number SSBzero; // unfished spawning biomass
  number fA50; // age at 50% selectivity for the fishery
  number sA50; // age at 50% selectivity for the survey
  number fbeta; // sloe in selectivity for the fishery
  number sbeta; // slope in selectiv ity for the survey
  number sigmaR; // standard deviation in sigma R

  sdreport_vector biomass(fyear,lyear+1); // model estimate of exploitable biomass

  // Temporary Variables
  number fpen;  // penalty to prevent catch > biomass
  number survival; // survival from removing catch for penalty
  // population dynamics variables
  vector Ffull(fyear,lyear); // fishing mortality for fully vunerable ages
  matrix F(fyear,lyear,fage,lage); // fishing mortality

  matrix N(fyear,lyear+1,fage,lage); // abundance
  vector SSB(fyear,lyear+1); //spawning biomass
  //Storage matrices
  matrix pred_survey(1,nsurvey,fage,lage); // Predicted survey numbers at age
  matrix pred_age_fishery(1,nage_fishery,fage,lage); // Predicted fishery catch at age

  objective_function_value obj_fun;    //value to be minimized

PROCEDURE_SECTION
  //dummy for testing
  obj_fun = square(dummy);

  // Parameter Transformations
  transform_params(); // no arguments passed, computer will automatically use global variables from PARAMETER Section

  // Selectivity
  get_selectivity();

  // Initial Conditions
  get_init();

  // Population Projection
  get_pop();

  // Objective Function
  nll_cpue();
  nll_survey();
  nll_ages();
  penalties_recruits();

 cout<< "total" << obj_fun << endl;

FUNCTION transform_params
  // Transform Parameters
  Rzero = mfexp(logRzero);
  //q = mfexp(logq);
  //sigma = mfexp(logsigma);
  h = 0.2 + 0.8*mfexp(logit_h)/(1.+mfexp(logit_h));  //your h was allowing range of 0.-1.0
  fracinit = mfexp(logit_fracinit)/(1.+mfexp(logit_fracinit));
  fA50 = mfexp(log_fA50);
  sA50 = mfexp(log_sA50);
  fbeta = mfexp(log_fbeta);
  sbeta = mfexp(log_sbeta);
  sigmaR = mfexp(log_sigmaR);

  // Diagonstics may be turned on
  if (diag==1) {
    cout << "Rzero" << Rzero << endl;
    cout << "fA50" << fA50 << endl;
    cout << "fbeta" << fbeta << endl;
    cout << "sA50" << sA50 << endl;
    cout << "sbeta" << sbeta << endl;
  }

FUNCTION  get_selectivity
  // Determine Fishery  Selectivity
  f_sel = 1./(1.+mfexp(-1.*fbeta*(ages-fA50)));
  f_sel /= max(f_sel); // Divide vector by max value so standardized from zero to 1, improves interpretation
  // Determine Survey Selectivity
  s_sel = 1./(1.+mfexp(-1.*sbeta*(ages-sA50)));
  s_sel /= max(s_sel);  // Divide vector by max value so standardized from zero to 1, improves interpretation

  if (diag==1) {
    cout << "fishery selectivity" << f_sel << endl;
    cout << "survery selectivity" << s_sel << endl;
  }

FUNCTION get_init
  // Initial recruitment as R0
  N(fyear,fage) = Rzero*fracinit;

  //loop for first year abundance to derive equilibrium abundnace at age at F=0
  for (age=fage+1;age<=lage;age++) {
   N(fyear,age) = N(fyear,age-1)*mfexp(-1.*M);
  }
  // Plus group
  N(fyear,lage) = N(fyear,lage)/(1.-mfexp(-1.*M));

  // Exploitable Biomass
  biomass(fyear) = sum(elem_prod(elem_prod(N(fyear)*mfexp(-0.5*M),f_sel),wt));
  // Spawning Biomass
  SSB(fyear) = sum(elem_prod(elem_prod(N(fyear),mat),wt));
  // Unfished Spawning Biomass
  SSBzero = SSB(fyear)/fracinit;

  if (diat==2) {
    cout<< fyear << " " << SSB(fear) << " " << biomass (fyear) << " " << N(fyear) << endl;
  }

FUNCTION get_pop
  // Loop for series recruitment, older abundance, biomass, SSB and F
  for (year=fyear+1;year<=lyear+1;year++) {
    get_N_at_age();
    if (diat==1) {
      cout << "bio" << " " << obj_fun << endl;
    }
  }

  if (diag == 1) {
    cout << "bio " << obj_fun << endl;
  }

// ??????????????????????????????????????????????????? Start Here to fill in holes
FUNCTION get_N_at_age
  // Annual Population Update
  
   //get the F rate
   fpen = 0.;
   Ffull(year-1) = -1.*log(posfun((1.-cat(year-1,2)/biomass(year-1)),0.01,fpen));
   obj_fun += 1000*fpen;
   F(year-1) = Ffull(year-1)*sel;  //not sure why this is needed
   for (age=fage+1;age<=lage;age++) {
    N(year,age) = N(year-1,age-1)*mfexp(-1.*(M+F(year-1,age-1)));
   }  //
   N(year,lage) += N(year-1,lage)*mfexp(-1.*(M+F(year-1,lage)));
   //biomass(year) = N * sel * wt;    /// matrix * vector * vector !!!!!
   //SSB(year) = N * mat * wt;
   biomass(year) = sum(elem_prod(elem_prod(N(year)*mfexp(-0.5*M),sel),wt));
   //recruitment
   N(year,fage) = 0.;
   SSB(year) = sum(elem_prod(elem_prod(N(year),mat),wt));
   N(year,fage) = 4.*h*Rzero*SSB(year-1)/(SSBzero*(1.-h)+SSB(year-1)*(5.*h-1.));
   //cout << year << " " << SSB(year) << " " << biomass(year) << " " << N(year) << endl;
  }   


FUNCTION nll_cpue
  //MLE for catchability
  f_q = 0.;
  for (icpue=1;icpue<=ncpue;icpue++) {
   f_q += log(cpue(icpue))-log(biomass(cpueyear(icpue))*mfexp(-0.5*Ffull(cpueyear(icpue))));
  }
  f_q = mfexp(f_q/ncpue);

  // MLE for sigma
  f_sigma = 0.;
  for(icpue=1;icpue<=ncpue;icpue++) {
    f_sigma += square(log(cpue(icpue))-log(f-q)-log(biomass(cpueyear(icpue))*mfexp(-0.5*Ffull(cpueyear(icpue)))));
  }
  f_sigma = sqrt(f_sigma/ncpue);

  obj_fun += ncpue*log(f_sigma) + 0.5*ncpue;

  if (diag==1) {
    cout << "cpue " << obj_fun << endl;
  }

FUNCTION nll_survey
  // Vulnerable biomass for survey (predicted value for survey)
  for (icpue=1;icpue<=nsurv;icpue++) {
    pred_survey(icpue) = elem_prod(elem_prod(N(agecomp_s(icpue,1)),mfexp(-0.5*(M+F(survey(icpue,1))))),s_sel);
     // 0.5 bit mean survey takes place in middle of year
     // Multiplied lots of things as vector calculations (first item in vector A * first item in vector B...), alternatively could loop over age
    pred_survey_bio(icpue) = sum(elem_prod(pred_survey(icpue),wt));
  }

// ??????????????????????????????????????????????????? Start Here to fill in holes

  // MLE for catchability
  s_q = 0.;
  for (icpue=1;icpue<=nsurvey;icpue++) {
    s_q +=

  obj_fun += nsurvey*log(s_sigma
  

//////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compare model predictions to data
  for (icpue=1;icpue<=ncpue;icpue++) {
   obj_fun += logsigma + square(log(cpue(icpue))-log(q)-log(biomass(cpueyear(icpue))))/(2*square(sigma));
  }

  cout << obj_fun << endl;
  
GLOBALS_SECTION
  #include <admodel.h>

REPORT_SECTION
//  report model estimates
  report << "biomass" << endl;
  for (year=fyear;year<=(lyear+1);year++)
   report << year << " " << biomass(year) << endl;
   
  report << "parameters" << endl;
  report << "h" << " " << h << endl;
  report << "R0" << " " << Rzero << endl;
  report << "fracinit" << " " << fracinit << endl;
  report << "sigma" << " " << sigma << endl;
  report << "q" << " " << q << endl;



