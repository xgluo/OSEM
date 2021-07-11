#include <Rcpp.h>
#include <random>
#include "tree.h"
#include "bits.h"
using namespace Rcpp;
using namespace std;
using namespace pcart;


// This function calls the pcart optimization algorithm
// Predictor type can only be CAT
// Response type can be CAT, ORD, or BDEU_CAT
// [[Rcpp::export]]
List cli_pcart(NumericMatrix rdata, size_t predictorCount, size_t dataCount,
               std::vector<size_t> pred_catCounts, std::string res_type, size_t res_catCount,
               double alpha = 0.5, double kappa = 0.25, double ess = 1.0, bool useStructureScore = true) {
  
  vector<VarPtr> predictors;
  for(size_t i = 0; i < predictorCount; ++i) {
    stringstream varName;
    varName << i;
    predictors.push_back(createCatVar(varName.str(), i, vector<string>(pred_catCounts[i])));
  }

  VarPtr response;
  string type = res_type;
  if(type == "CAT") {
    response = createCatVar("response", predictorCount, vector<string>(res_catCount), alpha);
  } else if(type == "ORD") {
    response = createOrdVar("response", predictorCount, vector<string>(res_catCount), alpha, kappa);
  } else if(type == "BDEU_CAT") {
    response = createBDeuCatVar("response", predictorCount, vector<string>(res_catCount), ess);
  } else {
    fail("Unknown response variable type ", type);
  }

  vector< vector<double> > data(dataCount);
  for(size_t i = 0; i < dataCount; ++i) {
    data[i].resize(predictorCount + 1);
    NumericVector temp = rdata.row(i);
    vector<double> temp2 = as< vector<double> >(temp);
    data[i] = temp2;
  }
  
  TreeResult res = optimizeTree(predictors, response, data, useStructureScore);

  List L = List::create(Named("dataScore") = res.dataScore,
                        Named("structureScore") = res.structureScore);
  
  return L;
}
