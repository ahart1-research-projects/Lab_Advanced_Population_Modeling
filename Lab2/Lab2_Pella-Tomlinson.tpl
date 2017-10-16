//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling Homework 1
//  Fall 2017
//  Amanda Hart
//
//  Biomass Dynamics: Pella-Tomlinson production 
//////////////////////////////////////////////
DATA_SECTION // Read in Data 
  init_int fyear; // First year of catch
  init_int lyear; // Last year of catch
  !! int nyear = lyear - fyear + 1; // Number of years
  init_int ncpue;
  init_matrix cat(fyear,lyear,1,2); // Matrix of Catch
  init_ivector cpueyear(1,ncpue); //index of biomasses for which we have CPUE
  init_vector cpue(1,ncpue);

  vector logcpue(1,ncpue);   

  number eps;

  int year // looping variables
  int icpue // looping variables
  
  !!cout << nyear << endl; // Prints something to screen when code is run
  !!cout << ncpue << endl;
  !!cout << cat(fyear,1) << " " << cat(lyear,2) << endl;
  !!cout << cpue(1) << endl;
  !!cout << cpue(ncpue) << endl;

PARAMETER_SECTION
 
  init_number dummy(-1); // -1 means use initial value and don't change Estimated parameter that is not used in real objective function, but is estimated so code can be compiled and tested

  init_number logr;     // Log of growth rate
  init_number logK;     // Log of carrying capacity
  // init_number logq;     // Log of catchability, not necessary to estimate since max likelihood of q calculated
  // init_number logsigma; //Log of sigma
  init_number logit_p;  // Logit of initial biomass, logit requires that number be between 0 and 1
  //init_number logm+1; // This is the shape parameter for the production model m > 1 required
  init_number transform_m(-1); 
  // log normal transformation can not be less than 0

  number r;
  number K;
  number q;
  number sigma;
  number p;
  // number m;
  //likeprof_number m; // Attempt to make likelihood distribution of m
 
  sdreport_vector biomass(fyear,lyear+1); //sdreport_ means this(the estimated biomass) will be stored as output vector

  number bpen;
  number survival;
  number fpen;

  sdreport_number BMSY; //sdreport_ means this will be stored as output
  sdreport_number MSY;
  sdreport_number UMSY;
  sdreport_number BKratio; // report standard deviation of B in last year:K ratio
  sdreport_number m; // report std dev of m

  objective_function_value obj_fun;




PROCEDURE_SECTION
  //dummy for testing
  obj_fun = square(dummy);


  //parameter transformations
  r = mfexp(logr);  // mfexp is a function similar to exp but is more stable
  K = mfexp(logK);
  // q = mfexp(logq);  // Not needed since max likelihood of q is calculated
  // sigma = mfexp(logsigma); // Not needed since max likelihood of sigma is calculated
  p = mfexp(logit_p) / (1.+mfexp(logit_p)); // Constrains value of p, 1. ensures computer knows it is a real number (and not an integer)
  //m = mfexp(logm-1); // Constrains the value of m to be greater than 1
  m = (mfexp (transform_m)+1); // logm allow 0 to 1, adding 1 allows m to be 1 to infinity

  // cout << logit_p << " " << p << endl;
  
  // initialize biomass in 1935
  biomass(fyear) = p*K; // Biomass in first year is p * Carrying capacity

  // loop over years to do the biomass dynamics 
  for (year = fyear+1; year <= lyear+1; year++){ // Loop from first year +1 to last year +1 loop by year each time
    bpen = 0.;
    // Below is the pella-tomlinson production model
    biomass(year) = biomass(year-1) + posfun(biomass(year-1) * r * (1. - pow((biomass(year-1)/K),(m-1.))),0.1,bpen);
    // posfun is a function which keeps everything within the function positive ?? where do you find help/ function descriptions
    obj_fun += 1000*bpen;
    survival = 1.-cat(year-1,2)/biomass(year); // survival = 1 - exploitation rate(proportion of total biomass which was caught), pick second column of cat
    fpen = 0.;
    biomass(year) *= posfun(survival, 0.01, fpen);
    obj_fun += 1000*fpen; // Some sort of equivalency thing ???? look again

  }
  //cout << biomass << endl; // check biomass correct

  // observation model

  // constant added to cpue predictions
  eps = 1.e-07;

  //MLE for catchability (max likelihood estimate)
  q = 0.;
  for (icpue=1;icpue<=ncpue;icpue++){
    q += log(cpue(icpue))-log(biomass(cpueyear(icpue))); // cpueyears gives index of biomass
  } // += adds the result to the previously recorded value of q
  q = mfexp(q/ncpue); // This divides the above by n and makes it exponent

  //cout << q << endl; // check q 

  // MLE for standard deviation of observation error (sigma)
 sigma=0.;
  for (icpue=1;icpue<=ncpue;icpue++) {
  sigma += square(log(cpue(icpue))-log(q)-log(biomass(cpueyear(icpue))));
  }
  sigma= sqrt(sigma/ncpue); 


  // compare model predictions to data (this uses the negative log likelihood function)
  for (icpue=1;icpue<=ncpue;icpue++){
    obj_fun += log(sigma) + square(log(cpue(icpue))-log(q)-log(biomass(cpueyear(icpue))))/(2*square(sigma));
  }
     // Now that this is part of obj_fun, ADMB by default is minimizing obj_fun (the objective function)

  cout << obj_fun << endl; // Print objective function each time (shows improvement in model fit)


  // calculate biological reference points (meaningful output)  
  BMSY = K/2;   // BMSY
  MSY = r*K/4;  // MSY
  UMSY = r/2;   // Exploitation rate at MSY
  BKratio = biomass(lyear+1)/K; // biomass in last year/K ratio
  // don't forget to declare variables that are added here



GLOBALS_SECTION
  #include <admodel.h>

REPORT_SECTION
//  report model estimates
  report << "biomass" << endl;
  for( year=fyear;year<=(lyear+1);year++){
    report << year << " " << biomass(year) << endl; // print year, space, biomass
    }

  report << "parameters" << endl; // print parameters
  report << "r" << " " << r << endl; // print year, space, biomass  //Estimated
  report << "K" << " " << K << endl; // print year, space, biomass  //Estimated
  report << "p" << " " << p << endl; // print year, space, biomass  //Estimated
  report << "sigma" << " " << sigma << endl; // print year, space, biomass   // calculated using MLE
  report << "q" << " " << q << endl; // print year, space, biomass     // calculated using MLE
  report << "m" << " " << m << endl;   // Estimated
  
  report << "Perfomance Metrics" << endl;
  report << "BMSY" << " " << BMSY << endl;
  report << "MSY" << " " << MSY << endl;
  report << "Exploitation Rate UMSY" << " " << UMSY << endl;
  report << "Biomass:K ratio" << " " << BKratio << endl;
