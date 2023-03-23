//#include <cmath>
//#include <list>
//#include <vector>

#include <fdeep/fdeep.hpp>

struct NNinfo {
  double sample_means[18];
  double sample_stds[18];
  double target_mean;
  double target_std;
  int prediction_iteration;
};

class NN{
  public:
  std::string nn_model_filename;
  NNinfo nn_info;
  // fdeep::model nn_model;
  
  void read_info_file( );
  NN(std::string nn_model_filename);
};