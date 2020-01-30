#include <Python/Python.h>
#include <tuple>
#include <vector>
#include <queue>

typedef std::tuple<int, double> edge_t;
typedef std::vector<edge_t> nodelist_t;
typedef std::vector<nodelist_t> edgelist_t;
typedef std::vector<double> double_vec_t;
typedef std::vector<bool> bool_vec_t;
typedef std::tuple<double_vec_t, double_vec_t> double_vec_tuple_t;

class queue_elem_t{
  public:
     int node;
     double weight;

     queue_elem_t(int node, double weight){
       this->node = node;
       this->weight = weight;
     }

     bool operator>(const queue_elem_t& rhs){
       return this->weight > rhs.weight;
     }
     bool operator>=(const queue_elem_t& rhs){
       return this->weight >= rhs.weight;
     }
     bool operator<(const queue_elem_t& rhs){
       return this->weight < rhs.weight;
     }
     bool operator<=(const queue_elem_t& rhs){
       return this->weight <= rhs.weight;
     }
     bool operator==(const queue_elem_t& rhs){
       return this->weight == rhs.weight;
     }
     bool operator!=(const queue_elem_t& rhs){
       return this->weight != rhs.weight;
     }
};
// getDegree expects a list of tuples (int, double)
// we first define the c-types version and then the
// overloaded python function
double getDegree(nodelist_t edges){
  double out = 0;
  for(edge_t edge:edges){
    out += std::get<1>(edge);
  }
  return out;
}

// c-type implementation of the push operation
bool_vec_t push(double_vec_t *p, double_vec_t *r,
                double alpha, int u, double du, nodelist_t neighbours,
                                  double_vec_t neighbours_deg, double epsilon){

  //update p and r
  (*p)[u] += alpha*(*r)[u];
  (*r)[u] *= (1.0-alpha)/2.;

  // initialize r_above_th
  bool_vec_t r_above_th = bool_vec_t(neighbours.size(), false);

  for(int i=0; i<neighbours.size(); i++){
    int n = std::get<0>(neighbours[i]);
    (*r)[n] += (1.0-alpha)*(*r)[u]/(2.0*du);

    if((neighbours_deg[i] == 0 && (*r)[n] != 0) || (*r)[n]/neighbours_deg[i] >= epsilon){
      r_above_th[i] = true;
    }
  }

  return r_above_th;
}

// c-type implementation of the approximateSimrank
double_vec_tuple_t approximateSimrank(edgelist_t A, int v, double alpha,
            double epsilon, int max_iters=200, bool return_only_neighbours=false){

  int N = A.size();

  double_vec_t p = double_vec_t(N, 0.);
  double_vec_t r = double_vec_t(N, 0.);
  r[v] = 1;

  // define the min-pq as per python
  std::priority_queue<queue_elem_t, std::vector<queue_elem_t>, std::less<queue_elem_t>> pq;
  pq.push(queue_elem_t(v, getDegree(A[v])));

  double inv_epsilon = 1.0/epsilon;

  std::vector<int> v_neighs = std::vector<int>(A[v].size());
  for(int i=0; i<A[v].size(); i++){
    v_neighs[i] = std::get<0>(A[v][i]);
  }

  int iters = 0;
  while(pq.size()>0 && pq.top().weight <= inv_epsilon && iters<max_iters){
    int u = pq.top().node;
    pq.pop();

    nodelist_t neigh = A[u];

    double_vec_t neigh_deg = double_vec_t(neigh.size(), 0);
    double du = 0;
    for(int i=0; i<neigh.size(); i++){
      int n = std::get<0>(neigh[i]);
      du += std::get<1>(neigh[i]);
      neigh_deg[i] = getDegree(A[n]);
    }

    bool_vec_t r_above_th = push(&p, &r, alpha, u, du, neigh, neigh_deg, epsilon);

    for(int i=0; i<r_above_th.size(); i++){
      if(r_above_th[i]){
        int n = std::get<0>(neigh[i]);
        pq.push(queue_elem_t(n, neigh_deg[i]/r[n]));
      }
    }

    iters++;
  }
  return double_vec_tuple_t(p, r);
}


edgelist_t localPageRank(edgelist_t A, double c, double epsilon=1e-5,
      int max_iters=200, bool return_only_neighbours=false){

  int N = A.size();
  edgelist_t L = edgelist_t(N);

  double alpha = 2*c/(1+c);
  // andersen's paper inverts alpha
  alpha = 1-alpha;

  for(int i=0; i<N; i++){
    // out = (p, r)
    double_vec_tuple_t out = approximateSimrank(A, i, alpha,  epsilon, max_iters, return_only_neighbours);

    // create the new nodelist
    L[i] = nodelist_t();

    for(int k=0; k<A[i].size(); k++){
      if(return_only_neighbours){
        L[i].push_back(edge_t(std::get<0>(A[i][k]), std::get<0>(out)[k]));
      }else{
        L[i].push_back(edge_t(std::get<0>(A[i][k]), std::get<0>(out)[std::get<0>(A[i][k])]));
      }
    }

  }
  return L;
}
