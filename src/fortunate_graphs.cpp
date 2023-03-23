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
  if (argc < 3) {
    std::cout << "\n\terror: please use " << std::string(argv[0]) << " n t [r] [alpha] [beta]\n\n";
    return 1; 
  }
  
  int n(std::stoi(argv[1]));
  int num_graphs(std::stoi(argv[2]));

  // density and target probability default values
  double r = 0.25;
  double alpha = 1;
  double beta = 2.0;
  
  // overwrite by command line arguments (if specified)
  if (argc >= 4) r = std::stod(argv[3]);
  if (argc >= 5) alpha = std::stod(argv[4]);  
  if (argc >= 6) beta = std::stod(argv[5]);      

  mcpu_timer t1, t2, t3, t4;
  time_stats s1, s2, s3, s4;

  std::cout << std::endl;

  for(int i = 0; i < num_graphs; i++){

    std::cout << "\r proceesing graph: " << std::setw(5) << (i + 1 ) << "/" << num_graphs << std::flush;

    auto G = Graph(n, r);
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

    // small check if all functions return the same value
    auto values = { p_dijk, p_prun, p_orac, p_orpr };
    for(auto lhs : values){
      for(auto rhs: values){
        if(lhs != rhs) {
          printf("\n\nError: encountered mismatch of computed distances:\n");
          printf("\tdijk: %6.2f \n", p_dijk); 
          printf("\tprun: %6.2f \n", p_prun);
          printf("\torac: %6.2f \n", p_orac);
          printf("\torpr: %6.2f \n", p_orpr);          
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


printf("\nMore detailed breakdown:\n\n");
printf("\tused parameters: n = %d, t = %d, r = %.3f, alpha = %.2f, beta = %.2f\n\n", n, num_graphs, r, alpha, beta);
  
#ifndef OUTPUT_TIMING

  printf("\tNumber of operations:\n\n");
  printf("\t               qin      qrm      qis      qdp     qcum     numt     rin      rrm      ris      rdp\n");
  printf("\tdijk (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s1.q_c[0]), s1.q_c[1]/double(s1.q_c[0]), s1.q_c[2]/double(s1.q_c[0]), s1.q_c[3]/double(s1.q_c[0]), s1.cum_q/double(s1.q_c[0]), s1.num_t/double(s1.q_c[0]));    
  printf("\tprun (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s2.q_c[0]), s2.q_c[1]/double(s2.q_c[0]), s2.q_c[2]/double(s2.q_c[0]), s2.q_c[3]/double(s2.q_c[0]), s2.cum_q/double(s2.q_c[0]), s2.num_t/double(s2.q_c[0]));    
  printf("\torac (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s3.q_c[0]), s3.q_c[1]/double(s3.q_c[0]), s3.q_c[2]/double(s3.q_c[0]), s3.q_c[3]/double(s3.q_c[0]), s3.cum_q/double(s3.q_c[0]), s3.num_t/double(s3.q_c[0]));    
  printf("\torpr (c): %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", double(s4.q_c[0]), s4.q_c[1]/double(s4.q_c[0]), s4.q_c[2]/double(s4.q_c[0]), s4.q_c[3]/double(s4.q_c[0]), s4.cum_q/double(s4.q_c[0]), s4.num_t/double(s4.q_c[0]), double(s4.r_c[0]), s4.r_c[1]/double(s4.r_c[0]), s4.r_c[2]/double(s4.r_c[0]), s4.r_c[3]/double(s4.r_c[0])); 
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
  printf("\taverage number of buckets orpr:     buck_c = %8.2f\n", s4.buck_c/double(s4.q_c[0]));

  printf("\tNote:\n");
  printf("\t- (t) row states CPU time (cumulative)\n");
  printf("\t- (c) row states operations (average)\n");  
  printf("\t- qin(c) = number of times shortest path computation was done\n");
  printf("\t- additional CPU time consumed outside of time-tracking\n\n");

#endif 

  return 0;
}
