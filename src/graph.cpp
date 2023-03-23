
#include "graph.hpp"

unsigned int graph_seed = 1;
std::mt19937 graph_rng(graph_seed);


void Graph::reset_timers() {

  // initialize and start timers
  qin_t = mcpu_timer();
  qrm_t = mcpu_timer();
  qis_t = mcpu_timer();
  qdp_t = mcpu_timer();

  rin_t = mcpu_timer();
  rrm_t = mcpu_timer();
  ris_t = mcpu_timer();
  rdp_t = mcpu_timer();

  pre_t = mcpu_timer();

  // initialize counters
  qin_c = 0;
  qrm_c = 0;
  qis_c = 0;
  qdp_c = 0;

  rin_c = 0;
  rrm_c = 0;
  ris_c = 0;
  rdp_c = 0;  

  buck_c = 0;
  qcum_c = 0;
  num_t = 1;
}


Graph::Graph(){

  reset_timers();

  n = 0;
  m = 0;
  incoming = std::vector<std::vector<Edge>>();
  outgoing = std::vector<std::vector<Edge>>();
  target_set = std::vector<bool>();

  max_edge_weight = 0;
  sp_upper_bound = 0;

}

void Graph::read_from_stream(std::istream & stream){
  // std::cout << "Reading graph file... ";
  std::string s;
  stream >> s; // 'nodes:'
  stream >> n; 
  stream >> s; // 'directed'
  stream >> s; // '---'

  target_set = std::vector<bool>(n, false);
  outgoing = std::vector<std::vector<Edge>>(n, std::vector<Edge>());
  incoming = std::vector<std::vector<Edge>>(n, std::vector<Edge>());
  int u;
  for(int i=0; i<n; i++){
    stream >> u >> s;
    if(s == "True"){
      target_set[u] = true;
    }
  }
  stream >> s; // '---'
  int v; double l;
  while(stream >> u){
    stream >> v >> l;
    outgoing[u].push_back({u, v, l});
    incoming[v].push_back({u, v, l});
    max_edge_weight = (max_edge_weight < l ? l : max_edge_weight);
    m++;
  }
  sp_upper_bound = (n-1)*max_edge_weight;
  // std::cout << "Done reading graph with " << n << " nodes" << std::endl;
}

Graph::Graph(std::istream & stream){
  reset_timers();
  read_from_stream(stream);
}

// Graph::Graph(std::string fileName, std::string fileformat, std::shared_ptr<fdeep::model> nn_model_, std::shared_ptr<NN> nn_model_info_){
Graph::Graph(std::string fileName, std::string fileformat){

  reset_timers();

  std::ifstream filestream;
  filestream.open(fileName);
  if(fileformat == "TestSet"){
    read_from_stream(filestream);
  }
  filestream.close();

}

Graph::Graph(const int n_, const double r ){
  
  n = n_;
  int x = (int) (n * r);
  // std::cout << "x: " << x << std::endl;
  m = ( x-1 ) * ( 1 + n-x ); //number of nodes which have outgoing edges ( x-1 ) * number of outgoing edges ( 1+ n-x)

  incoming = std::vector<std::vector<Edge>>(n, std::vector<Edge>());
  outgoing = std::vector<std::vector<Edge>>(n, std::vector<Edge>());

  target_set = std::vector<bool>(n, false);
  target_set[x-1] = true; //Only x-1 is a target
  
  //std::uniform_real_distribution<double> Noise(0.0, 1/(double)(x-1));
  
  for( int u = 0; u < x-1; u++ ){
    Edge e_sp = {u, u+1, 1/(double)(x-1) };
    outgoing[u].push_back(e_sp);
    incoming[u+1].push_back(e_sp);
    for( int v = x; v < n; v++ ){
      double noise = 0.0; //Noise(graph_rng);
      Edge e_out = {u, v, 2 - 2*u/(double)(x-1) + noise };
      outgoing[u].push_back(e_out);
      incoming[v].push_back(e_out);
    }
  }
  //std::cout << "Generated graph with " << n << " nodes and " << m << " edges" << std::endl;
  
  //max_edge_weight = 2;
  sp_upper_bound = 3;
}


Graph::Graph(const int n_, const double avg_degree, const double q){
  // fast generation of G(n,c/n)
  // based on fast_gnp_random_graph from networkx (python library) 
  // https://networkx.org/documentation/stable/_modules/networkx/generators/random_graphs.html#fast_gnp_random_graph
  // runs in O(n + m) = O(n + c*n) time

  // GS
  // MODIFIED: slightly tweaked below to ensure that we get "suitable" random graphs: 
  //   - s should not be a target node (to avoid shortest path distance 0) 
  //   - target set should not be empty (choose a random node then)
  //   - at least one target node in T is reachable from s (to avoid shortest path distance infinity)
  //


  do { // repeat until random graph has been constructed in which s is connected to T

    reset_timers();

    n = n_;
    m = 0;

    incoming = std::vector<std::vector<Edge>>(n, std::vector<Edge>());
    outgoing = std::vector<std::vector<Edge>>(n, std::vector<Edge>());

    target_set = std::vector<bool>(n, false);

    double p = avg_degree/n;
    int w = -1;
    double lp = log(1.0-p);

    std::uniform_real_distribution<double> L(0.0, 1.0);

    int v = 0;
    while(v < n){
      double lr = log(1.0 - L(graph_rng));
      w = w + 1 + int(lr/lp);
    if(v==w) w++; // avoid selfloops
    while(v < n && n <= w){
      w = w-n;
      v = v+1;
      if(v==w) w++; // avoid selfloops
    }
    if(v < n){
      m++;
      double l = L(graph_rng);
      max_edge_weight = (max_edge_weight < l ? l : max_edge_weight);
      Edge e = {v, w, l};
      outgoing[v].push_back(e);
      incoming[w].push_back(e);
    }
  }
  
  sp_upper_bound = (n-1)*max_edge_weight;

  // generate target set
  int no_targets = 0;
  while (no_targets == 0) {
    for (int v = 0; v < n; v++) target_set[v] = false;
      for(int v=0; v < n; v++) {
        if(L(graph_rng) < q) {
          target_set[v] = true;
          no_targets++;
        }
      }
      if (target_set[0]) {
        target_set[0] = false;
        no_targets--;
      }
    }
  } while (bfsPathLength(0, 0) == -1); 
}


double Graph::predictBFS( int s, bool weighted ){
  
  double bfs_path_length = bfsPathLength( s, weighted );
  
  if( weighted ) {
    return bfs_path_length;
  } else {
    double average_edge_length = 0.5;
    return average_edge_length * bfs_path_length;
  }
}


double Graph::predict(int s, std::string pred_method, std::vector<std::pair<double, double>> raw_features){

  if( pred_method == "bfs" ) {
    return predictBFS( s, false );
  } 
  else if( pred_method == "bfs_weighted" ) {
    return predictBFS( s, true );
  }
  else {
    fdeep::tensor features(fdeep::tensor_shape(18), 0);
    for(int i=0; i < 9; i++){
      if(raw_features[i+1].second == raw_features[0].second){
        double second_value = (0 - nn_model_info->nn_info.sample_means[i+9]) / nn_model_info->nn_info.sample_stds[i+9];
        features.set(fdeep::tensor_pos(i+9), second_value); // filter out infinity
      }
      else{
        double second_value = (raw_features[i+1].second - nn_model_info->nn_info.sample_means[i+9]) / nn_model_info->nn_info.sample_stds[i+9];
        features.set(fdeep::tensor_pos(i+9), second_value);
      }
      double first_value = (raw_features[i+1].first - nn_model_info->nn_info.sample_means[i]) / nn_model_info->nn_info.sample_stds[i]; //raw_features[i+1].first
      features.set(fdeep::tensor_pos(i), first_value);
    }
    
    const auto result = nn_model->predict({features});
    

    
    double P = result[0].to_vector()[0];
    //Normalize P
    P = (P * nn_model_info->nn_info.target_std ) + nn_model_info->nn_info.target_mean;
    
    return P;
  }
}


double Graph::bfsPathLength( int s, bool weighted ){
  // returns number (weighted = 0) or total weight (weighted = 1) of edges on BFS-path 

  std::queue<int> q;
  
  q.push(s);
  auto path_lengths = std::vector<int>(n, INT_MAX);
  path_lengths[s] = 0.0;
  
  auto path_weights = std::vector<double>(n, DBL_MAX);
  path_weights[s] = 0.0;
  
  while (!q.empty()){
    auto u = q.front(); q.pop();
    if(target_set[u]){
      if (weighted){
        return (double)path_weights[u];
      } 
      else {
        return (double)path_lengths[u];
      }
    }
    for(auto e : outgoing[u]){
      if(path_lengths[u] + 1 < path_lengths[e.v]){
        q.push( e.v );
        path_lengths[e.v] = path_lengths[u] + 1;
        path_weights[e.v] = path_weights[u] + e.l;
      }
    }
  }
  return -1.0;
}


Path Graph::backtrack_path(const int s, const int t, const std::vector<double> &D){
  // extract actual shortest path by backtracking
  auto p = Path();
  int v = t;
  while(v != s){
    for(auto e: incoming[v]){
      if(D[e.u] + e.l == D[e.v]){ // if edge is tight we can use it as part of the shortest path 
        p.edges.push_back(e);
        v = e.u;
        break;
      }
    }
  }
  std::reverse(p.edges.begin(), p.edges.end());
  return p;
}


double Graph::dijkstra(const int s){
  // solve the single-source many-targets shortsest path problem using adapted Dijkstra
  
  reset_timers();
  if(target_set[s]) return 0; // nothing to do

  qin_c++;
  qin_t.tik();
  auto D = std::vector<double>(n, DBL_MAX);
  D[s] = 0;
  
  auto visited = std::vector<bool>(n, false);
  visited[s] = true;

  #ifdef FIBONACCI
    boost::heap::fibonacci_heap<std::pair<double, int>> q; // priority_queue with underlying fibonacci heap
    typedef typename boost::heap::fibonacci_heap<std::pair<double, int>>::handle_type handle_t; 
  #elif BINOMIAL
    boost::heap::binomial_heap<std::pair<double, int>> q; // priority_queue with underlying binomial heap
    typedef typename boost::heap::binomial_heap<std::pair<double, int>>::handle_type handle_t; 
  #endif 
  auto handles = std::vector<handle_t>(n);  // vector to store handles so that we can update them 

  qis_c++;
  handles[s] = q.push({0, s});
  int t = -1;
  qin_t.tok();
  
  while(!q.empty()) {

    qrm_c++;
    qrm_t.tik();
    auto [dist, u] = q.top(); q.pop();
    dist = -dist;
    qrm_t.tok();

    qcum_c += q.size();
  
    if(target_set[u]) { t = u; break; }

    for(auto e : outgoing[u]) {
      double tent = D[u] + e.l;

      if(tent >= D[e.v]) continue; 

      D[e.v] = tent;
      if(visited[e.v]){
        qdp_c++;
        qdp_t.tik();
        q.increase(handles[e.v], {-D[e.v], e.v});
        qdp_t.tok();
      }
      else{
        qis_c++;
        qis_t.tik();
        handles[e.v] = q.push({-D[e.v], e.v});
        visited[e.v] = true;
        qis_t.tok();
      }
    }
  }

  if (t != -1) return D[t]; 
  else return DBL_MAX;

}


double Graph::dijkstra_pruning(const int s, double B){
  // solve the single-source many-targets shortsest path problem using adapted Dijkstra with pruning
  
  reset_timers();
  if(target_set[s]) return 0; // nothing to do

  qin_c++;
  qin_t.tik();
  auto D = std::vector<double>(n, DBL_MAX);
  D[s] = 0;
  
  auto visited = std::vector<bool>(n, false);
  visited[s] = true;

  // switch (use -D when compiling) to use Fibonacci or binomial heap for priority queue 
  #ifdef FIBONACCI
    boost::heap::fibonacci_heap<std::pair<double, int>> q; // priority_queue with underlying fibonacci heap
    typedef typename boost::heap::fibonacci_heap<std::pair<double, int>>::handle_type handle_t; 
  #elif BINOMIAL
    boost::heap::binomial_heap<std::pair<double, int>> q; // priority_queue with underlying binomial heap
    typedef typename boost::heap::binomial_heap<std::pair<double, int>>::handle_type handle_t; 
  #endif 

  auto handles = std::vector<handle_t>(n);  // vector to store handles so that we can update them 

  qis_c++;
  handles[s] = q.push({0, s});
  int t = -1;
  qin_t.tok();
  
  while(!q.empty()) {

    qrm_c++;
    qrm_t.tik();
    auto [dist, u] = q.top(); q.pop();
    dist = -dist;
    qrm_t.tok();

      qcum_c += q.size();
    
    if(target_set[u]) { t = u; break; }

    for(auto e : outgoing[u]) {
      double tent = D[u] + e.l;

      if ( (tent >= D[e.v]) || (tent > B) ) continue;

      D[e.v] = tent;

      // update B if e.v in target set
      if(target_set[e.v]) B = tent;

      if(visited[e.v]) {
        qdp_c++;
        qdp_t.tik();
        q.increase(handles[e.v], {-D[e.v], e.v});
        qdp_t.tok();          
      }
      else {
        qis_c++;
        qis_t.tik();
        handles[e.v] = q.push({-D[e.v], e.v});
        visited[e.v] = true;
        qis_t.tok();          
      }
    }
  }

  if (t != -1) return D[t]; 
  else return DBL_MAX;

}


double Graph::dijkstra_oracle(const int s, double B) {
  // call dijkstra_pruning with B set to actual shortest path distance 

  return dijkstra_pruning(s, B);
}


double Graph::dijkstra_prediction_oracle(const int s, double alpha, double beta, double sp_value){
  // call dijkstra_prediction with prediction set to alpha times actual shortest path distance from beginning (prediction_iteration = 1)

  return dijkstra_prediction(s, 1, alpha, beta, sp_value*alpha, "");
}


double Graph::dijkstra_prediction_ML(const int s, const int prediction_iteration, double alpha, double beta, std::shared_ptr<fdeep::model> nn_model_,  std::shared_ptr<NN> nn_model_info_){
  // call dijkstra_prediction and use ML to obtain prediction in iteration prediction_iteration

  nn_model = nn_model_;
  nn_model_info = nn_model_info_;
  return dijkstra_prediction(s, prediction_iteration, alpha, beta, DBL_MAX, "ML");
}


double Graph::dijkstra_prediction_BFS(const int s, const int prediction_iteration, double alpha, double beta){
  // call dijkstra_prediction and use edge-BFS to obtain prediction in iteration prediction_iteration

  return dijkstra_prediction(s, prediction_iteration, alpha, beta, DBL_MAX, "bfs");
}


double Graph::dijkstra_prediction_BFS_w(const int s, const int prediction_iteration, double alpha, double beta){
  // call dijkstra_prediction and use weighted-BFS to obtain prediction in iteration prediction_iteration

  return dijkstra_prediction(s, prediction_iteration, alpha, beta, DBL_MAX, "bfs_weighted");
}


double Graph::dijkstra_prediction(const int s, const int prediction_iteration, double alpha, double beta, double P, std::string pred_method){
  // solves the single-source many targets problem by using an adapted Dijkstra algorithm with pruning and prediction
  // prediction P is either assumed to be given up front (prediction_iteratioj = 0) or determined by specific method in iteration prediction_iteration

  reset_timers();
  if (P == 0) return 0;
  if(target_set[s]) return 0; // nothing to do

  // initialize B to infinity  
  double B = DBL_MAX;
  
  // vector to keep track of features
  auto raw_features = std::vector<std::pair<double, double>>(prediction_iteration);

  // create reserve list (bucket list), actual init happens in iteration prediciton_iteration below
  bucket_list r;
  
  qin_c++;
  qin_t.tik();
  // vector to store (tentative) distances 
  auto D = std::vector<double>(n, DBL_MAX);
  D[s] = 0;

  // vector to keep track whether a node has been visited already
  auto visited = std::vector<bool>(n, false);
  visited[s] = true;
  
  // switch (use -D when compiling) to use Fibonacci or binomial heap for priority queue 
  #ifdef FIBONACCI
    boost::heap::fibonacci_heap<std::pair<double, int>> q; // priority_queue with underlying fibonacci heap
    typedef typename boost::heap::fibonacci_heap<std::pair<double, int>>::handle_type handle_t; 
  #elif BINOMIAL
    boost::heap::binomial_heap<std::pair<double, int>> q; // priority_queue with underlying binomial heap
    typedef typename boost::heap::binomial_heap<std::pair<double, int>>::handle_type handle_t; 
  #endif 

  // vector to store handle_t of each node to allow for direct access 
  auto handles = std::vector<handle_t>(n);  

  qis_c++;
  // insert source node 
  handles[s] = q.push({0, s});
  int t = -1;
  qin_t.tok();
  
  // iteration counter
  int iteration = 0;
  
  // std::cout << std::endl;
  // std::cout << "start prediction" << std::endl;

  do {

    // std::cout << "Begin trial with P equals " << P << std::endl;

    // stop if q runs empty or the minimum element is larger than P (needed because of reserve list)
    while(!q.empty() && -q.top().first <= P) {

      qrm_c++;
      qrm_t.tik();
      // extract minimum from q
      auto [dist, u] = q.top(); q.pop();
      dist = -dist;
      qrm_t.tok();
      
      qcum_c += q.size();

      // check whether a target node has been found
      if(target_set[u]) { t = u; return D[t]; }

      // std::cout << "Remove min: " << dist << " u: " << u << std::endl;

      iteration++;
      if(iteration <= prediction_iteration){
        // add current dist and B values to trace
        raw_features[iteration-1] = {dist, B};
      }
      if(iteration == prediction_iteration){

        // retrieve prediction based on trace (unless prediction_iteration = 1, from beginning)
        if (prediction_iteration != 1) {
          pre_t.tik();
          P = alpha * predict(0, pred_method, raw_features);
          // std::cout << " Prediction: " << P << std::endl;
          pre_t.tok();
        }

        // initialize reserve list (bucket list) based on P
        rin_c++;
        rin_t.tik();
        // determine max number of buckets needed based on bound sp_upper_bound
        int bn = ceil ( log (sp_upper_bound / P) / log( beta ) ) + 1;
        buck_c += bn;
        // std::cout << " (use " << bn << " buckets) " << std::flush;
        r.init(n, P, beta, bn);
        rin_t.tok(); 
      }

      // relax all outgoing edges of u
      for(auto e : outgoing[u]) {

        auto v = e.v;
        double tent = D[u] + e.l;

        // skip if new tentative distance tent is worse or can be pruned 
        if ( (tent >= D[v]) || (tent > B) ) continue;

        D[v] = tent;

        // update B if e.v in target set
        if(target_set[v]) B = tent;

        // std::cout << "   v: " << v << " tent: " << tent << " visited: " << visited[v] << "r.contains(v): " << r.contains(v) << "target_set[v]" << target_set[v] << std::endl;

        if(visited[v]) {
          if(!r.contains(v)) { 
            // already in q
            qdp_c++;
            qdp_t.tik();
            q.increase(handles[v], {-tent, v});   // increase priority
            qdp_t.tok();
          }
          else {
            // currently in r
            if(tent > P) {
              // stays in r
              rdp_c++;
              rdp_t.tik();
              // note: ensures that new priority tent > P
              r.move(v, tent); 
              rdp_t.tok();
            }
            else {
              // move from r to q (remove from r, insert in q)
              rrm_c++;
              rrm_t.tik();
              r.remove(v);            
              rrm_t.tok();
              qis_c++;
              qis_t.tik();
              handles[v] = q.push({-tent, v}); 
              qis_t.tok();
            }
          }
        }
        else { // e.v encountered for the first time
          visited[v] = true;
          if(tent <= P) { 
            // add v to PQ
            qis_c++;
            qis_t.tik();
            handles[v] = q.push({-tent, v});
            qis_t.tok();
          }
          else { 
            // add v to r
            ris_c++;
            ris_t.tik();
            // note: ensures that priority tent > P
            r.insert(v, tent);
            ris_t.tok();
          }
        }
      }
    }


    // std::cout << "Trial done. P equals " << P << "q.empty(): " << q.empty() << "r.size()" << r.size << std::endl;
    if(!q.empty()) assert(-q.top().first > P);

    // get smallest distance from q (if non-empty)
    double q_front = DBL_MAX;
    if(!q.empty()){
      q_front = -q.top().first;
      //std::cout<< "q not empty. -q.top().first: " <<  -q.top().first << " Prediction: " << P << " R empty: " << r.empty() << std::endl;
  	}
      
    // check for batch insertions
    if (!r.empty() or !q.empty() ) {

      // increases lowest bucket idx_0 until non-empty list, returns number of increments
      int i = r.advance(q_front); 
      // printf("\nentered batch insertion: P = %f, B = %f", P, B);
      // printf("\nadv = %d, bucket = %d", i, r.idx_0);
      // printf("\nmoving these nodes: ");      
  	  // std::cout<< "r.empty(r.idx_0)" << r.empty(r.idx_0) << std::endl;

      if (i == 0) continue;   // no non-empty bucket found for batch insertion
      else {
        // update P accordingly 
        P *= pow(beta, i);
        num_t += i;

        if (beta*B <= P) continue;   // all remaining items in used buckets exceed B; can stop

        // remove all items from lowest bucket and insert into q (unless > B)

        // r.bucket[r.idx_0]->print();

        while (!r.empty(r.idx_0)) {
          // remove first item
          rrm_c++;
          rrm_t.tik();
          int v = r.pop(r.idx_0);
          // std::cout << v << "-"<< std::flush;
          rrm_t.tok();

          if( D[v] <= B ) {
            qis_c++;
            qis_t.tik();
            handles[v] = q.push({-D[v], v});
            // std::cout << "+ " << std::flush;            
            qis_t.tok();
          }
          else {
            // std::cout << " " << std::flush;
            // mark node as unvisited because it is neither in q nor in r
            visited[v] = false;
          }
        }
      }
    }
  } while (!q.empty());

  return DBL_MAX;
}

