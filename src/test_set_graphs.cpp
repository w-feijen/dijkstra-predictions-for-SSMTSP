//Code by Ruben Brokkelkamp, Willem Feijen and Guido Sch\"afer
//March 2023

#include "graph.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>

#define ZERO std::chrono::steady_clock::steady_clock::duration::zero()

int main(int argc, char *argv[]) {

  // read command line arguments (needed for n and num_graphs)
  if (argc < 2) {
    std::cout << "\n\terror: please use " << std::string(argv[0]) << " t [alpha] [beta]\n\n";
    return 1; 
  }
  
  int num_graphs(std::stoi(argv[1]));

  // density and target probability default values
  double alpha = 1;
  double beta = 1.05;
  
  // overwrite by command line arguments (if specified)
  if (argc >= 3) alpha = std::stod(argv[2]);
  if (argc >= 4) beta = std::stod(argv[3]);

  mcpu_timer t1, t2, t3, t4, t5, t6, t7;
  time_stats s1, s2, s3, s4, s5, s6, s7;

  std::cout << std::endl;
  
  // read neural network and corresponding data
  std::string nn_model_filename;
  nn_model_filename = "based_on_27-01-2021_iterations_zeros_no_bool_no_norm2021_01_27_1725";
  std::shared_ptr<fdeep::model> nn_model;
  nn_model = std::make_shared<fdeep::model>(fdeep::load_model("../models/fdeep_model_" + nn_model_filename + ".json"));
  std::shared_ptr<NN> nn_model_info;
  nn_model_info = std::make_shared<NN>(NN(nn_model_filename));

  std::string prefix = "../Graphs/test_set/";
  std::string suffix = ".txt";
  int start_graph = 10000;

  for(int i = 0; i < num_graphs; i++){

    std::cout << "\r proceesing graph: " << std::setw(5) << (i + 1 ) << "/" << num_graphs << std::flush;

    std::string graph_name = prefix + std::to_string(i + start_graph) + suffix ;
    auto G = Graph(graph_name, "TestSet");
    int s = 0;

    // DIJKSTRA
    t1.tik();
    auto p_dijk = G.dijkstra(s);
    t1.tok();
    s1.add(G); 
	
    // DIJKSTRA PRUNING
    t2.tik();
    auto p_prun = G.dijkstra_pruning(s);
    t2.tok();
    s2.add(G);

    // extract shortest path distance for oracle
    double sp_value = p_dijk;
    if (sp_value == 0) { 
      std::cerr << "\n\n\t error: shortest path distance is 0\n\n";
    }
    else if (sp_value == DBL_MAX) {
      std::cerr << "\n\n\t error: shortest path distance is infinity\n\n";
    }

    // DIJKSTRA ORACLE
    t3.tik();
    auto p_orac = G.dijkstra_oracle(s, sp_value);
    t3.tok();
    s3.add(G); 

    // DIJKSTRA ORACLEPREDICTION
    t4.tik();
    auto p_orpr = G.dijkstra_prediction_oracle(s, alpha, beta, sp_value);
    t4.tok();
    s4.add(G); 

    // DIJKSTRA PREDICTION
    t5.tik();
    auto p_pred = G.dijkstra_prediction_ML(s, 10, alpha, beta, nn_model, nn_model_info);
    t5.tok();
    s5.add(G); 

    // DIJKSTRA PREDICTION BFS
    t6.tik();
    auto p_bfs = G.dijkstra_prediction_BFS(s, 10, alpha, beta);
    t6.tok();
    s6.add(G); 

    // DIJKSTRA PREDICTION BFS-w
    t7.tik();
    auto p_bfsw = G.dijkstra_prediction_BFS_w(s, 10, alpha, beta);
    t7.tok();
    s7.add(G); 

    // small check if all functions return the same value
    auto values = { p_dijk, p_prun, p_orac, p_orpr, p_pred, p_bfs, p_bfsw };
    for(auto lhs : values){
      for(auto rhs: values){
        if(lhs != rhs) {
          printf("\n\nError: encountered mismatch of computed distances:\n");
          printf("\tdijk: %6.2f \n", p_dijk); 
          printf("\tprun: %6.2f \n", p_prun);
          printf("\torac: %6.2f \n", p_orac);
          printf("\torpr: %6.2f \n", p_orpr);          
          printf("\tpred: %6.2f \n", p_pred);          
          printf("\tbfs:  %6.2f \n", p_bfs);          
          printf("\tbfsw: %6.2f \n", p_bfsw);          
          break;
        }
      }
    }
  }

std::cout << " done!" << std::endl;

// output statistics below (old school printf pretty)
printf("\nTotal runtime summary:\n\n");  
printf("\tdijk: %6.2f sec\n", t1.out_sec());
printf("\tprun: %6.2f sec\n", t2.out_sec());  
printf("\torac: %6.2f sec\n", t3.out_sec());
printf("\torpr: %6.2f sec\n", t4.out_sec());     
printf("\tpred: %6.2f sec\n", t5.out_sec());     
printf("\tbfs:  %6.2f sec\n", t6.out_sec());     
printf("\tbfsw: %6.2f sec\n", t7.out_sec());     


printf("\nMore detailed breakdown:\n\n");
printf("\tused parameters: n = %d, t = %d, c = %.0f, q = %.3f, alpha = %.2f, beta = %.2f\n\n", 1000, num_graphs, 8.0, 0.02, alpha, beta);

#ifndef OUTPUT_TIMING

  printf("\tNumber of operations:\n\n");
  printf("\t               qin      qrm      qis      qdp     qcum     numt     rin      rrm      ris      rdp\n");
  printf("\tdijk (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s1.q_c[0]), s1.q_c[1]/double(s1.q_c[0]), s1.q_c[2]/double(s1.q_c[0]), s1.q_c[3]/double(s1.q_c[0]), s1.cum_q/double(s1.q_c[0]), s1.num_t/double(s1.q_c[0]));    
  printf("\tprun (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s2.q_c[0]), s2.q_c[1]/double(s2.q_c[0]), s2.q_c[2]/double(s2.q_c[0]), s2.q_c[3]/double(s2.q_c[0]), s2.cum_q/double(s2.q_c[0]), s2.num_t/double(s2.q_c[0]));    
  printf("\torac (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s3.q_c[0]), s3.q_c[1]/double(s3.q_c[0]), s3.q_c[2]/double(s3.q_c[0]), s3.q_c[3]/double(s3.q_c[0]), s3.cum_q/double(s3.q_c[0]), s3.num_t/double(s3.q_c[0]));    
  printf("\torpr (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s4.q_c[0]), s4.q_c[1]/double(s4.q_c[0]), s4.q_c[2]/double(s4.q_c[0]), s4.q_c[3]/double(s4.q_c[0]), s4.cum_q/double(s4.q_c[0]), s4.num_t/double(s4.q_c[0]), double(s4.r_c[0]), s4.r_c[1]/double(s4.r_c[0]), s4.r_c[2]/double(s4.r_c[0]), s4.r_c[3]/double(s4.r_c[0])); 
  printf("\tpred (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s5.q_c[0]), s5.q_c[1]/double(s5.q_c[0]), s5.q_c[2]/double(s5.q_c[0]), s5.q_c[3]/double(s5.q_c[0]), s5.cum_q/double(s5.q_c[0]), s5.num_t/double(s5.q_c[0]), double(s5.r_c[0]), s5.r_c[1]/double(s5.r_c[0]), s5.r_c[2]/double(s5.r_c[0]), s5.r_c[3]/double(s5.r_c[0])); 
  printf("\tbfs  (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s6.q_c[0]), s6.q_c[1]/double(s6.q_c[0]), s6.q_c[2]/double(s6.q_c[0]), s6.q_c[3]/double(s6.q_c[0]), s6.cum_q/double(s6.q_c[0]), s6.num_t/double(s6.q_c[0]), double(s6.r_c[0]), s6.r_c[1]/double(s6.r_c[0]), s6.r_c[2]/double(s6.r_c[0]), s6.r_c[3]/double(s6.r_c[0])); 
  printf("\tbfsw (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n\n", double(s7.q_c[0]), s7.q_c[1]/double(s7.q_c[0]), s7.q_c[2]/double(s7.q_c[0]), s7.q_c[3]/double(s7.q_c[0]), s7.cum_q/double(s7.q_c[0]), s7.num_t/double(s7.q_c[0]), double(s7.r_c[0]), s7.r_c[1]/double(s7.r_c[0]), s7.r_c[2]/double(s7.r_c[0]), s7.r_c[3]/double(s7.r_c[0])); 
  printf("\taverage number of buckets orpr:     buck_c = %8.2f\n\n", s4.buck_c/double(s4.q_c[0]));


  printf("\tNote:\n");
  printf("\t- each row states the average number of operations executed\n");  
  printf("\t- qin = number of times shortest path computation was done \n\n");

#else

  printf("\tNumber of operations and time needed:\n\n");  
  printf("\t               qin      qrm      qis      qdp     qcum     numt     rin      rrm      ris      rdp\n"); 
  printf("\tdijk (t): %8.2f %8.2f %8.2f %8.2f\n", s1.q_t[0], s1.q_t[1], s1.q_t[2], s1.q_t[3]);  
  printf("\tdijk (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s1.q_c[0]), s1.q_c[1]/double(s1.q_c[0]), s1.q_c[2]/double(s1.q_c[0]), s1.q_c[3]/double(s1.q_c[0]), s1.cum_q/double(s1.q_c[0]), s1.num_t/double(s1.q_c[0]));    
  printf("\tprun (t): %8.2f %8.2f %8.2f %8.2f\n", s2.q_t[0], s2.q_t[1], s2.q_t[2], s2.q_t[3]);    
  printf("\tprun (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s2.q_c[0]), s2.q_c[1]/double(s2.q_c[0]), s2.q_c[2]/double(s2.q_c[0]), s2.q_c[3]/double(s2.q_c[0]), s2.cum_q/double(s2.q_c[0]), s2.num_t/double(s2.q_c[0]));    
  printf("\torac (t): %8.2f %8.2f %8.2f %8.2f\n", s3.q_t[0], s3.q_t[1], s3.q_t[2], s3.q_t[3]);      
  printf("\torac (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s3.q_c[0]), s3.q_c[1]/double(s3.q_c[0]), s3.q_c[2]/double(s3.q_c[0]), s3.q_c[3]/double(s3.q_c[0]), s3.cum_q/double(s3.q_c[0]), s3.num_t/double(s3.q_c[0]));    
  printf("\torpr (t): %8.2f %8.2f %8.2f %8.2f                   %8.2f %8.2f %8.2f %8.2f\n", s4.q_t[0], s4.q_t[1], s4.q_t[2], s4.q_t[3], s4.r_t[0], s4.r_t[1], s4.r_t[2], s4.r_t[3]);
  printf("\torpr (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s4.q_c[0]), s4.q_c[1]/double(s4.q_c[0]), s4.q_c[2]/double(s4.q_c[0]), s4.q_c[3]/double(s4.q_c[0]), s4.cum_q/double(s4.q_c[0]), s4.num_t/double(s4.q_c[0]), double(s4.r_c[0]), s4.r_c[1]/double(s4.r_c[0]), s4.r_c[2]/double(s4.r_c[0]), s4.r_c[3]/double(s4.r_c[0])); 
  printf("\tpred (t): %8.2f %8.2f %8.2f %8.2f                   %8.2f %8.2f %8.2f %8.2f\n", s5.q_t[0], s5.q_t[1], s5.q_t[2], s5.q_t[3], s5.r_t[0], s5.r_t[1], s5.r_t[2], s5.r_t[3]);
  printf("\tpred (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s5.q_c[0]), s5.q_c[1]/double(s5.q_c[0]), s5.q_c[2]/double(s5.q_c[0]), s5.q_c[3]/double(s5.q_c[0]), s5.cum_q/double(s5.q_c[0]), s5.num_t/double(s5.q_c[0]), double(s5.r_c[0]), s5.r_c[1]/double(s5.r_c[0]), s5.r_c[2]/double(s5.r_c[0]), s5.r_c[3]/double(s5.r_c[0])); 
  printf("\tbfs  (t): %8.2f %8.2f %8.2f %8.2f                   %8.2f %8.2f %8.2f %8.2f\n", s6.q_t[0], s6.q_t[1], s6.q_t[2], s6.q_t[3], s6.r_t[0], s6.r_t[1], s6.r_t[2], s6.r_t[3]);
  printf("\tbfs  (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s6.q_c[0]), s6.q_c[1]/double(s6.q_c[0]), s6.q_c[2]/double(s6.q_c[0]), s6.q_c[3]/double(s6.q_c[0]), s6.cum_q/double(s6.q_c[0]), s6.num_t/double(s6.q_c[0]), double(s6.r_c[0]), s6.r_c[1]/double(s6.r_c[0]), s6.r_c[2]/double(s6.r_c[0]), s6.r_c[3]/double(s6.r_c[0])); 
  printf("\tbfsw (t): %8.2f %8.2f %8.2f %8.2f                   %8.2f %8.2f %8.2f %8.2f\n", s7.q_t[0], s7.q_t[1], s7.q_t[2], s7.q_t[3], s7.r_t[0], s7.r_t[1], s7.r_t[2], s7.r_t[3]);
  printf("\tbfsw (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n\n", double(s7.q_c[0]), s7.q_c[1]/double(s7.q_c[0]), s7.q_c[2]/double(s7.q_c[0]), s7.q_c[3]/double(s7.q_c[0]), s7.cum_q/double(s7.q_c[0]), s7.num_t/double(s7.q_c[0]), double(s7.r_c[0]), s7.r_c[1]/double(s7.r_c[0]), s7.r_c[2]/double(s7.r_c[0]), s7.r_c[3]/double(s7.r_c[0])); 
  printf("\taverage number of buckets orpr:     buck_c = %8.2f\n", s4.buck_c/double(s4.q_c[0]));
  printf("\ttime pred used for prediction: pred   = %8.2f \n", s5.pre_t);
  printf("\ttime bfs  used for prediction: pred   = %8.2f \n", s6.pre_t);
  printf("\ttime bfsw used for prediction: pred   = %8.2f \n\n", s7.pre_t);

  printf("\tNote:\n");
  printf("\t- (t) row states CPU time (cumulative)\n");
  printf("\t- (c) row states operations (average)\n");  
  printf("\t- qin(c) = number of times shortest path computation was done\n");
  printf("\t- additional CPU time consumed outside of time-tracking\n\n");

#endif 

  return 0;
}
