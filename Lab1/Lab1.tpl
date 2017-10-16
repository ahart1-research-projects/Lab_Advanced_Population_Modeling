//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2017
//  Gavin Fay
//
//  Lab 1: Linear regression example
//////////////////////////////////////////////
DATA_SECTION
  init_int ndata;    //number of data points
  init_vector x(1,ndata);   // x data values
  init_vector y(1,ndata);   // y data values

PARAMETER_SECTION
  init_number a;     //intercept parameter
  init_number b;     //slope parameter 
  
  vector ypred(1,ndata);   //predicted values for y
  objective_function_value obj_fun; //objective function value: sum of squares

PROCEDURE_SECTION
  ypred = a+b*x;   //calculate the predictions for y given current parameter values
  obj_fun = norm2(y-ypred);   //work out the sum of squares
  
