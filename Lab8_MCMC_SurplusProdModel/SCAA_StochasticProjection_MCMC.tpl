//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2017
//  Steve Cadrin & Gavin Fay
//
//  Statistical catch at age model with bayesian
//  Projections in this model give a dynamics equilibrium
//  time series observation error model assuming Beverton-Holt stock-recruitment
//////////////////////////////////////////////
DATA_SECTION
  init_int diag;  //print diagnostics?
  init_int fyear;   // first year of catch data
  init_int lyear;   // last year of catch data
  init_int fage; // youngest age
  init_int lage; // oldest age
  init_vector mat(fage,lage); // proportion mature
  init_vector initsel(fage,lage); // selectivity
  init_vector wt(fage,lage); // mean weight
  init_int frecdev;
  init_int lrecdev;
  int nrecdevs;
  !! nrecdevs = lrecdev-frecdev+1;
  init_int ncpue;   // number of years of CPUE data
  init_matrix cat(fyear,lyear,1,2);   //matrix of catch
  init_ivector cpueyear(1,ncpue);   // years for cpue data points
  init_vector cpue(1,ncpue);   // cpue data
  init_int nsurvey;       // number of years of survey index data
  init_matrix survey(1,nsurvey,1,2);  //survey index
  init_int nage_survey;  //number of years of age composition info
  int nage;   //number of ages
  !! nage = lage-fage+1;

  init_matrix agecomp_s(1,nage_survey,1,nage+2);   // matrix of agecomp info
  init_int nage_fishery;                           // number of yrs of fishery age comp data
  init_matrix agecomp_f(1,nage_fishery,1,nage+2);  // matrix of fishery age comps

  init_int nprojyears; // number of years in projection
  init_number projFfull; // projection F at full recruitment

  ivector agecols(1,nage);    //ivector for storing the columns that contain the age data
  ivector ages(1,nage);       // index vector for computing things like selectivity
  !! for (age=1;age<=nage;age++) agecols(age) = 2 + age;
  !! for (age=1;age<=nage;age++) ages(age) = fage + age - 1;
  

  int year;     // looping variable
  int age;    // looping variable
  int icpue;    // looping variable
  
  number eps;   // small constant
  !! eps = 1.e-07;
  // print some things to screen to check data read in correctly
//  cout << ncpue << endl;
//  cout << cat(fyear,1) << " " << cat(lyear,2) << endl;
//  cout << cpue(1) << endl;
//  cout << cpue(ncpue) << endl;
  vector v(lyear+1,lyear+nprojyears+1); // This is the vector for random numbers

   !! cout << agecomp_f(nage_fishery) << endl ;

  !!CLASS ofstream MCMCreport("SCAA_MCMCreport.out", ios::trunc); // C++ statement, object of class: ofstream, named SCAA_MCMCreport
  

 // LOCAL_CALCS
 // This does C++ calculations locally, LOCAL_CALCS & calculations in this section must be indented 1 space
 // random_number_generator rng(911); // 911 is the seed number
 // v.fill_randn(rng); // function fills vector v. with standard normal deviates (mean = 0, std dev =1) with seed rng
 //END_CALCS 

PARAMETER_SECTION
  init_number dummy(-1);    // dummy parameter for debugging. -1 phase to turn off estimation

  init_number M(2);      
  init_number logit_h(-3);  //logit of steepness
  init_number logRzero(1);  //log of asymptotic recruitment
  init_number logit_fracinit(-1);
  //init_number logq;  //log of catchability
  //init_number logsigma(1); //log of sigma
  init_number log_fA50(1);   //log of age at 50% selex for fishery
  init_number log_fbeta(1);  //log of slope in selectivity for fishery
  init_number log_sA50(1);   // log of age at 50% selex for survey
  init_number log_sbeta(1);  // log of slope in selectivity for survey
  init_number log_sigmaR(-1); // log of standard deviation of recruitment devs
  init_bounded_vector recdevs(frecdev,lrecdev,-10,10,2);  //recruitment deviations

  vector f_sel(fage,lage); // selectivity of fishery
  vector s_sel(fage,lage); // selectivity of survey

  number h;    //steepness
  number Rzero;    //asymptotic recruitment
  number fracinit; //initial yr recruitment
  number s_q;    //catchability
  number s_sigma;   //observation error std. deviation
  number f_q;    //catchability
  number f_sigma;   //observation error std. deviation
  number SSBzero;  //unfished spawning biomass
  number fA50;   //age at 50% selectivity for the fishery
  number sA50;   // age at 50% selectivity for the survey
  number fbeta;  // slope in selectivity for the fishery
  number sbeta;  // slope in selectivity for the survey
  number sigmaR;  // standard deviation in sigma R

  sdreport_vector biomass(fyear,lyear+1); //model predictions of exploitable biomass

  // temporary variables
  number fpen;  // penalty to prevent catch > biomass
  number survival; // survival from removing catch for penalty
  
// population dynamics variables
  vector Ffull(fyear,lyear+1); // fishing mortality for fully vunerable ages through the first projected year (2017)
  matrix F(fyear,lyear,fage,lage); // fishing mortality
  matrix N(fyear,lyear+nprojyears+1,fage,lage); // abundance over time series and projected years, could make separate matrices
  sdreport_vector SSB(fyear,lyear+nprojyears+1); //spawning biomass over time series and projected years
  vector projF(fage,lage); // Projected F at age
  vector projY(lyear+1,lyear+1+nprojyears); // projected yield

  matrix pred_survey(1,nsurvey,fage,lage);  //predicted survey numbers at age
  vector pred_survey_bio(1,nsurvey);   //predicted survey biomass

  matrix pred_age_survey(1,nage_survey,fage,lage);  // model predicted survey catch at age
  matrix pred_age_fishery(1,nage_fishery,fage,lage); //model predicted fishery catch at age
  matrix proj_age_fishery(lyear+1,lyear+nprojyears,fage,lage); // Could have combined with pred_age_fishery matrix, chose to separate 

  sdreport_number SSBrel; //SSB in 2017 relative to SSBzero, report std deviation

  objective_function_value obj_fun;    //value to be minimized

PROCEDURE_SECTION
  //dummy for testing
  obj_fun = square(dummy);


  //parameter transformations
  transform_params();

  //selectivity
  get_selex();

  // initial conditions
  get_init();

  // population projection
  get_pop();

  //objective function;
  nll_cpue();
  nll_survey();
  nll_ages();
  pen_recruits();
 
  cout << "total " << obj_fun << endl;

  SSBrel =  SSB(2017)/SSBzero; // calculate SSBrel so its std dev is reported

  // Set up what we want reported in mceval_phase() as part of the analysis of MCMC
  if (mceval_phase() ) {
    // Bayesian analysis
       // Make sure the code is compiled before trying to do the MCMC lines below
       // Run the following two lines as part of ADMB -> Run with Args   to produce MCMC draws
       // -mcmc #draws -mcsave #saved (every ith draw saved)

       // THIS CODE
       // -mcmc 4500000 -mcsave 3000
       // in order to discard the first 1,500,000 cycles, only the last 1000 samples are saved


    // Projection calculations included in normal report (similar equations to earlier estimates, just use projected F (projF) instead
       // This doesn't alter the objective function so we can do projections in last phase
       // This code doesn't add to the objective function
       // Could loop over lyear+nprojyears+1 for N and SSB and loop over lyear+nprojyears for Yield and Catch
    //projF = projFfull*f_sel; // projected full F * selectivity at age, f_sel is an estimated vector
    projF = M*f_sel; // projected F = natural mortality * F selectivity at age
    for(year=lyear+2;year<=lyear+nprojyears;year++){
       for(age=fage+1;age<=lage;age++){ // loop over age
          N(year,age) = N(year-1,age-1)*mfexp(-1.*(M+projF(age-1))); // abundance at age = survival from previous age * survival in current age
          proj_age_fishery(year,age) = N(year,age)*projF(age)*(1.-mfexp(-1.*(M+projF(age))))/(M+projF(age));
       }
       // Create random number, starting seed number is just the year
       random_number_generator rng(year); // 911 is the seed number
       v.fill_randn(rng); // function fills vector v. with standard normal deviates (mean = 0, std dev =1) with seed rng

       N(year,lage) += N(year-1,lage)*mfexp(-1.*(M+projF(lage))); // Cumulatively adds to plus group abundance at age
       N(year,fage) = 4.*h*Rzero*SSB(year-1)/(SSBzero*(1.-h)+SSB(year-1)*(5.*h-1))*mfexp(v(year)*sigmaR-0.5*sigmaR*sigmaR); // Calculate first age group (recruits) using beverton-holt equation and the parameter estimates from above
          // Line above is vector of R * std deviates - bias correction, where std deviates = v(year) * sigma R = vector of std normal deviates times the sigmaR (std deviates for Recruitment)
          // bias correction is half of sigma squared (sigma^2)/2, need this if we want the mean rather than the median recruitment in the projection
          // the model biomass uses the mean rather than median so we make this correction for consistency
       projY(year) = sum(elem_prod(proj_age_fishery(year),wt));
       SSB(year) = sum(elem_prod(elem_prod(N(year),mat),wt));
    }

    // Report objective function, and all the estimated parameters in the MCMCreport
    MCMCreport << "objective function " << obj_fun << " M " << M  << " Rzero " << Rzero  << " SSB2017/SSBzero " << SSB(2017)/SSBzero  << " Ffull2017 " <<  Ffull(2016) << " AvgYield " << sum(projY(lyear+1,lyear+1+nprojyears))/(nprojyears)  << endl;    // ?? The above is wrong because r K and MSY are not calculated. What needs to be reported here ??
    // ?? I think this needs to be unobserved/estimated things (treated as random variables, but I am not sure ??
 ///////////////////////// START HERE AND COMMENT OUT LINE 191-192 ///////////////////////
   // Report 2017-2031 values for SSB2017 rel SSB0, F, Y, M
   for(year=lyear+1;year<=2031;year++){
     MCMCreport << year << " SSBrel " << SSB(year)/SSBzero << " Yield " << projY(year) << " F=M " << M << endl;
   }
   ////////////////////// THIS PART USE FOR QUESTION 4 /////////////////////////////

 // Report projected Yield
 /* MCMCreport << "year projected_F projected_SSB projected_Yield" << endl;
  for (year=lyear+2;year<=lyear+nprojyears+1;year++){
   MCMCreport << year << " " << projFfull << " " << SSB(year) << " " << projY(year) << endl; 
   // In last year SSB reported only if projection looped over nproj+1, no projFfull, no projY reported
  } */


  }
 

FUNCTION transform_params
  //transform parameters
  Rzero = mfexp(logRzero);
  //sigma = mfexp(logsigma);
  h = 0.2 + 0.8*mfexp(logit_h)/(1.+mfexp(logit_h));  //your h was allowing range of 0.-1.0
  fracinit = mfexp(logit_fracinit)/(1.+mfexp(logit_fracinit));
  fA50 = mfexp(log_fA50);
  sA50 = mfexp(log_sA50);
  fbeta = mfexp(log_fbeta);
  sbeta = mfexp(log_sbeta);
  sigmaR = mfexp(log_sigmaR);
  //cout << logit_p << " " << p << endl;
  if (diag==1) cout << "Rzero " << Rzero << endl;
  if (diag==1) cout << "fA50 " << fA50 << endl;
  if (diag==1) cout << "fbeta " << fbeta << endl;
  if (diag==1) cout << "sA50 " << sA50 << endl;
  if (diag==1) cout << "sbeta " << sbeta << endl;

FUNCTION get_selex
///selex = 1./(1+e^(-1.*beta*(ages-A50)))
  f_sel = 1./(1.+mfexp(-1.*fbeta*(ages-fA50)));
  s_sel = 1./(1.+mfexp(-1.*sbeta*(ages-sA50)));
  f_sel /= max(f_sel);
  s_sel /= max(s_sel);  
  if (diag==1) {
    cout << "fishery selex " << f_sel << endl;
    cout << "survey selex " << s_sel << endl;
  }

FUNCTION get_init
  //initial recruitment
  N(fyear,fage) = Rzero*fracinit; 

  //loop for first year abundance to derive equilibrium abundance at age at F=0
  for (age=fage+1;age<=lage;age++)
    N(fyear,age) = N(fyear,age-1)*mfexp(-1.*M);
  //plus group
  N(fyear,lage) = N(fyear,lage)/(1.-mfexp(-1.*M));

  //exploitbale biomass
  biomass(fyear) = sum(elem_prod(elem_prod(N(fyear)*mfexp(-0.5*M),f_sel),wt));
  //spawningbiomass
  SSB(fyear) = sum(elem_prod(elem_prod(N(fyear),mat),wt));
  //unfished spawning biomass
  SSBzero = SSB(fyear)/fracinit;

   if (diag==2) cout << fyear << " " << SSB(fyear) << " " << biomass(fyear) << " " << N(fyear) << endl;


FUNCTION get_pop
  // loop for series recruitment, older abundance, biomass, SSB and F
  for (year=fyear+1;year<=lyear+1;year++) {
     get_Natage();
     if (diag==1) cout << year << " " << N(year) << endl;
  }
  
  //cout << biomass << endl;
  if (diag==1) cout << "bio " << obj_fun << endl;

FUNCTION get_Natage
   // annual population update
   //get the F rate
   //if (diag==1) cout << year << " " << -1.*log(posfun((1.-cat(year-1,2)/biomass(year-1)),0.01,fpen)) << " " << 1.-cat(year-1,2)/biomass(year-1)<< " " << cat(year-1,2) << " " << biomass(year-1) << endl;
   fpen = 0.;
   Ffull(year-1) = cat(year-1,2)/(eps+biomass(year-1));
   Ffull(year-1) = -1.*log(posfun(1.-Ffull(year-1),0.01,fpen));
   obj_fun += 1000*fpen;

   F(year-1) = Ffull(year-1)*f_sel;
   
   if (diag==1) cout << cat(year-1,2) << " " << biomass(year-1) << " " << Ffull(year-1) << endl;
   if (diag==1) cout << F(year-1) << endl;

   for (age=fage+1;age<=lage;age++)
     N(year,age) = N(year-1,age-1)*mfexp(-1.*(M+F(year-1,age-1)));
   N(year,lage) += N(year-1,lage)*mfexp(-1.*(M+F(year-1,lage)));

   //recruitment
   N(year,fage) = 4.*h*Rzero*SSB(year-1)/(SSBzero*(1.-h)+SSB(year-1)*(5.*h-1.));
   if (year>=frecdev & year<=lrecdev) N(year,fage) *= mfexp(recdevs(year)-0.5*sigmaR*sigmaR);

   //spawning biomass 
   SSB(year) = sum(elem_prod(elem_prod(N(year),mat),wt));

   if (diag==2) cout << year << " " << SSB(year) << " " << biomass(year) << " " << N(year) << endl;
     
   biomass(year) = sum(elem_prod(elem_prod(N(year)*mfexp(-0.5*M),f_sel),wt));
   //cout << year << " " << biomass(year) << " " << N(year) << " " << f_sel << endl;

  
FUNCTION nll_cpue
  //MLE for catchability
  f_q = 0.;
  for (icpue=1;icpue<=ncpue;icpue++) {
   f_q += log(cpue(icpue))-log(biomass(cpueyear(icpue))*mfexp(-0.5*Ffull(cpueyear(icpue))));
  }
  f_q = mfexp(f_q/ncpue);
  
  //MLE for sigma
  f_sigma = 0.;
  for (icpue=1;icpue<=ncpue;icpue++) {
   f_sigma += square(log(cpue(icpue))-log(f_q)-log(biomass(cpueyear(icpue))*mfexp(-0.5*Ffull(cpueyear(icpue)))));
  }
  f_sigma = sqrt(f_sigma/ncpue);
  
  obj_fun += ncpue*log(f_sigma) + 0.5*ncpue;
  if (diag==1) cout << "cpue " << obj_fun << endl;


FUNCTION nll_survey
  //vulnerable biomass for survey
  for (icpue=1;icpue<=nsurvey;icpue++) {
    pred_survey(icpue) = elem_prod(elem_prod(N(agecomp_s(icpue,1)),mfexp(-0.5*(M+F(survey(icpue,1))))),s_sel);
    pred_survey_bio(icpue) = sum(elem_prod(pred_survey(icpue),wt));
  }

  //MLE for catchability
  s_q = 0.;
  for (icpue=1;icpue<=nsurvey;icpue++)
    s_q += log(survey(icpue,2))-log(pred_survey_bio(icpue));
  s_q = mfexp(s_q/nsurvey);
  
  //MLE for sigma
  s_sigma = 0.;
  for (icpue=1;icpue<=nsurvey;icpue++)
    s_sigma += square(log(survey(icpue,2))-log(s_q)-log(pred_survey_bio(icpue)));
  s_sigma = sqrt(s_sigma/nsurvey);

  //cout << "s " << s_sigma << " " << nsurvey*log(s_sigma) + 0.5*nsurvey << endl; 

  obj_fun += nsurvey*log(s_sigma) + 0.5*nsurvey;
  if (diag==1) cout << "survey index " << obj_fun << endl;


FUNCTION nll_ages
  dvar_vector obs_age(1,nage);
  // fishery ages 
  for (int iobs=1;iobs<=nage_fishery;iobs ++) {
    year = agecomp_f(iobs,1);
    for (age=1;age<=nage;age++) {
      pred_age_fishery(iobs,age) = N(year,age)*F(year,age)*(1.-mfexp(-1.*(M+F(year,age))))/(M+F(year,age));   //Catch numbers at age, Baranov catch equation
      //pred_age_fishery(iobs,age) = N(year,age)*mfexp(-0.5*(F(year,age)+M))*f_sel(age);
    }
    pred_age_fishery(iobs) /= sum(pred_age_fishery(iobs));
    obs_age = agecomp_f(iobs)(agecols);
    obs_age /= sum(obs_age);

    for (age=1;age<=nage;age++) {    
      if (agecomp_f(iobs,2+age)>0)
//        obj_fun += -1.*agecomp_f(iobs,2)*sum(elem_prod(obs_age(age),log(pred_age_fishery(iobs,age))));
//        obj_fun += -1.*agecomp_f(iobs,2)*obs_age(age)*log(pred_age_fishery(iobs,age));
        obj_fun += -1.*agecomp_f(iobs,2)*obs_age(age)*log(pred_age_fishery(iobs,age)/obs_age(age));        
    }
  }
  if (diag==1) cout << "fishery ages " << obj_fun << endl;


  // survey ages 
  for (int iobs=1;iobs<=nage_survey;iobs ++) {
    year = agecomp_s(iobs,1);
    pred_survey(iobs) /= sum(pred_survey(iobs));
    //dvar_vector 
    obs_age = agecomp_s(iobs)(agecols);
    obs_age /= sum(obs_age);

    for (age=1;age<=nage;age++) {    
        if (agecomp_s(iobs,2+age)>0)
//        obj_fun += -1.*agecomp_s(iobs,2)*sum(elem_prod(obs_age(age),log(pred_survey(iobs,age))));
//        obj_fun += -1.*agecomp_s(iobs,2)*obs_age(age)*log(pred_survey(iobs,age));
        obj_fun += -1.*agecomp_s(iobs,2)*obs_age(age)*log(pred_survey(iobs,age)/obs_age(age));        
    }
  }
  if (diag==1) cout << "survey ages " << obj_fun << endl;



FUNCTION pen_recruits
  obj_fun += nrecdevs*log(sigmaR) + 0.5*norm2(recdevs)/square(sigmaR);
  if (diag==1) cout << "recruitment penalty  " << obj_fun << endl;


GLOBALS_SECTION
  #include <admodel.h>

REPORT_SECTION
  // Report estimated SSB
  report << "Year estimated_SSB" << endl;
    for (year=fyear;year<=lyear+1;year++){
     report << year << " " << SSB(year) << endl;
    }

  report << "parameters" << endl;
  report << "h" << " " << h << endl;
  report << "R0" << " " << Rzero << endl;
  report << "fracinit" << " " << fracinit << endl;
  report << "f_sigma" << " " << f_sigma << endl;
  report << "f_q" << " " << f_q << endl;
  report << "s_sigma" << " " << s_sigma << endl;
  report << "s_q" << " " << s_q << endl;

  report << "f_sel" << endl;
  report << f_sel << endl;
  report << "s_sel" << endl;
  report << s_sel << endl;
