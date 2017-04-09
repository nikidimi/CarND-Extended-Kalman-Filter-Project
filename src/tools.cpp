#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    return rmse;
  }
  for(int i=0; i < estimations.size(); ++i) {
    VectorXd sum = estimations[i] - ground_truth[i];
    rmse += (sum.array() * sum.array()).matrix();
  }
  
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  
  Hj << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
  
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  

  
  if (px == 0 && py == 0) {
    px = 0.001;
    py = 0.001;
  }
  
  float sum = px * px + py * py;
  float sum_sqrt = sqrt(sum);
  float sum_3 = sum * sum_sqrt;
  
  Hj(0, 0) = px / sum_sqrt;
  Hj(0, 1) = py / sum_sqrt;
  Hj(1, 0) = -py / sum;
  Hj(1, 1) = px / sum;
  Hj(2, 0) = py * (vx * py - vy * px) / sum_3;
  Hj(2, 1) = px * (vy * px - vx * py) / sum_3;
  Hj(2, 2) = Hj(0, 0);
  Hj(2, 3) = Hj(0, 1);
  
  return Hj;
  
}
