#define FDEEP_FLOAT_TYPE double

#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <fstream>
#include <climits>
#include <float.h>
#include <limits>
#include <queue>
#include <chrono>
#include <cmath>
#include <random>
#include <string>
#include <list>

#include <boost/timer/timer.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/binomial_heap.hpp>

#include "mcpu_timer.hpp"
#include "bucket_list.hpp"

#include "nn.hpp"
#include <fdeep/fdeep.hpp>


struct Edge {
  public:
  int u; // start
  int v; // end
  double l; // length

  friend std::ostream& operator<<(std::ostream& os, const Edge& e) {
    return os << e.u << " --> " << e.v << " (" << e.l << ")";
  }

  bool operator==(const Edge & rhs) const {
    return u == rhs.u && v == rhs.v && l == rhs.l; 
  }

  bool operator!=(const Edge & rhs) const { return !(*this == rhs); }
};

struct Path{
  std::vector<Edge> edges;

  bool operator==(const Path & rhs) const {
    if(edges.size() != rhs.edges.size()) return false;
    for(int i=0; i < edges.size(); i++){
      if(edges[i] != rhs.edges[i]) return false;
    }
    return true;
  }

  bool operator!=(const Path & rhs) const { return !(*this == rhs);}

  size_t size(){
    return edges.size();
  }
};

struct Features{
  double dist;
  std::vector<double> B;
  std::vector<double> D;
  double nr_its;
};



class Graph {
  public:
    int n; // number of nodes
    int m; // number of edges

    double max_edge_weight;
    double sp_upper_bound;

    std::vector<std::vector<Edge>> outgoing;
    std::vector<std::vector<Edge>> incoming;

    std::vector<bool> target_set;

    // operation timers and counters (q for priority queue, r for reserve list (if any))
    mcpu_timer qin_t, qrm_t, qis_t, qdp_t;
    mcpu_timer rin_t, rrm_t, ris_t, rdp_t;
    int        qin_c, qrm_c, qis_c, qdp_c;
    int        rin_c, rrm_c, ris_c, rdp_c;    
    // prediction timer
    mcpu_timer pre_t;

    // keep track of number of buckets used
    int buck_c;

    // keep track of cumulative queue size and number of trials
    int qcum_c, num_t;

    
    std::shared_ptr<fdeep::model> nn_model;
    std::shared_ptr<NN> nn_model_info;

    void read_from_stream(std::istream & stream);
    void reset_timers();

    Graph();
    // Graph(std::string fileName, std::string fileformat,std::shared_ptr<fdeep::model> nn_model,  std::shared_ptr<NN> nn_model_info);
    Graph(std::string fileName, std::string fileformat);
    Graph(std::istream & stream);
    Graph(const int n, const double c, const double f);
    Graph(const int n, const double r );

    double predict(int s , std::string pred_method, std::vector<std::pair<double, double>> raw_features);
    double predictBFS( int s, bool weighted );
    double bfsPathLength( int s, bool weighted );

    Path backtrack_path(const int s, const int t, const std::vector<double> &D);

    double dijkstra(const int s);
    double dijkstra_pruning(const int s, double B = DBL_MAX);
    double dijkstra_oracle(const int s, double B);
    double dijkstra_prediction(const int s, const int prediction_iteration, double alpha, double beta, double P, std::string pred_method);
    double dijkstra_prediction_oracle(const int s, double alpha, double beta, double sp_value);
    double dijkstra_prediction_ML(const int s, const int prediction_iteration, double alpha, double beta, std::shared_ptr<fdeep::model> nn_model,  std::shared_ptr<NN> nn_model_info);
    double dijkstra_prediction_BFS(const int s, const int prediction_iteration, double alpha, double beta);
    double dijkstra_prediction_BFS_w(const int s, const int prediction_iteration, double alpha, double beta);
};

// CPU times bookkeeping
struct time_stats { 

  // total times spent on priority queue 
  std::vector<double> q_t;

  // total times spent on reserve list 
  std::vector<double> r_t;

  // total time spent on prediction 
  double pre_t;

  // total number of priority queue operations counter
  std::vector<int> q_c;

  // total number of reserve list operations counter
  std::vector<int> r_c;

  // number of buckets used
  int buck_c;
  
  // total cumulative queue sizie
  int cum_q;
  
  // number of trials 
  int num_t;
  
  time_stats() {
    q_t = std::vector<double>(4, 0);
    r_t = std::vector<double>(4, 0);
    pre_t = 0;
    q_c = std::vector<int>(4,0);
    r_c = std::vector<int>(4,0);   
    buck_c = 0;
    cum_q = 0;
    num_t = 0;
  }

  // add measured times/counters in graph G
  void add(Graph &G) { 
    q_t[0] += G.qin_t.out_sec();
    q_t[1] += G.qrm_t.out_sec();
    q_t[2] += G.qis_t.out_sec();
    q_t[3] += G.qdp_t.out_sec();

    r_t[0] += G.rin_t.out_sec();
    r_t[1] += G.rrm_t.out_sec();
    r_t[2] += G.ris_t.out_sec();
    r_t[3] += G.rdp_t.out_sec();    

    pre_t += G.pre_t.out_sec();

    q_c[0] += G.qin_c;
    q_c[1] += G.qrm_c; 
    q_c[2] += G.qis_c;
    q_c[3] += G.qdp_c;

    r_c[0] += G.rin_c;
    r_c[1] += G.rrm_c; 
    r_c[2] += G.ris_c;
    r_c[3] += G.rdp_c;

    buck_c += G.buck_c;

    cum_q += G.qcum_c;
    num_t += G.num_t;
  }

};


