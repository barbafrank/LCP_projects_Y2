#include <Python/Python.h>
#include <tuple>
#include <vector>
#include <queue>

// getDegree expects a list of tuples (int, double)
// we first define the c-types version and then the
// overloaded python function
double getDegree(std::vector<std::tuple<int, double>> edges){
  double out = 0;
  for(auto edge:edges){
    out += std::get<1>(edge);
  }
  return out;
}

// c-type implementation of the push operation
std::vector<bool> push(std::vector<double> *p, std::vector<double> *r,
    double alpha, int u, double du, std::vector<std::tuple<int, double>> neighbours,
        std::vector<double> neighbours_deg, double){

  //update p and r
  (*p)[u] += alpha*(*r)[u];
  (*r)[u] *= (1.0-alpha)/2.;

  // initialize r_above_th
  std::vector<bool> r_above_th = std::vector<bool>(neighbours.size(), false);

  for(int i=0; i<neighbours.size(); i++){
    int n = std::get<0>(neighbours[i]);
    (*r)[n] += (1.0-alpha)*(*r)[u]/(2.0*du);

    if(neighbours_deg[i] == 0 && (*r)[i] != 0)
  }

  return r_above_th;
}
