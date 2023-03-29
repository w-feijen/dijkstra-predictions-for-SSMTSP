#include "nn.hpp"

NN::NN(std::string nn_model_filename_ ){
  nn_model_filename = nn_model_filename_;
  // nn_model = fdeep::load_model("fdeep_model_" + nn_model_filename + ".json");
  read_info_file();
}

void NN::read_info_file( ){
  //Read the info for the network with specified filename.
  //Info like sample_means, sample_stds is read
  
  std::string infofilename =  "../models/info_" + nn_model_filename + ".txt";
  
  //std::cout << "Reading infofile" << infofilename;
  
  std::ifstream stream;
  stream.open(infofilename);
  
  int prediction_iteration = 10;
  int nrFeatures = 2 * ( prediction_iteration - 1);
  
  //NNinfo nn_info;
  double mean;
  double std;
  std::string s;
  //int prediction_iteration_check;
  
  for(int i = 0; i< nrFeatures; i ++ ){
    stream >> mean;
    nn_info.sample_means[i] = mean;
  }
  
  stream >> s; //----
  
  for(int i = 0; i< nrFeatures; i ++ ){
    stream >> std;
    nn_info.sample_stds[i] = std;
  }
  
  stream >> s; //-----
  stream >> nn_info.target_mean;
  stream >> s; //-----
  stream >> nn_info.target_std;
  stream >> s; //-----
  stream >> nn_info.prediction_iteration;
  
  // for(auto i : nn_info.sample_means ){
    // std::cout << i << std::endl;
  // }
  
  // for(auto i : nn_info.sample_stds ){
    // std::cout << i << std::endl;
  // }
  // std::cout << nn_info.target_mean << nn_info.target_std << std::endl;
}

