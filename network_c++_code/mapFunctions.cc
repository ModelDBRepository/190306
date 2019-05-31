#include "mapFunctions.h"
#include <string>
#include <map>

double assign_value(Fmap some_map, string sname, double default_value){
  Fmap::iterator found=some_map.find(sname);
  if (found != some_map.end()){
    return (*found).second;
  }
  else{
    return default_value;
  }
}

