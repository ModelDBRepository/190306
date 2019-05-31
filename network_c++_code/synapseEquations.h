#include <string>
//#include <unordered_map>
#include "cellVariables.h"
#include "cellParameters.h"
using namespace std;

//typedef unordered_map<string,double> Fmap;
double synapse_gating(cellVariables&, cellParameters&);
double synapse_gating_decay_part(cellVariables&, cellParameters&);
double vmd_gating(cellVariables&, cellParameters&);
double vmd_gating2(cellVariables&, cellParameters&);
