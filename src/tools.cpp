#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
   VectorXd rmse(4);
   rmse << 0, 0, 0, 0;

   if (estimations.size() != ground_truth.size() || estimations.size() == 0)
   {
      std::cerr << "Invalid estimation or ground_truth data" << std::endl;
      return rmse;
   }

   // Calculate square sum of difference
   for (int i = 0; i < estimations.size(); i++)
   {
      VectorXd difference = estimations[i] - ground_truth[i];
      VectorXd difference_square = difference.array() * difference.array();
      rmse += difference_square;
   }

   rmse = rmse / estimations.size();
   rmse = rmse.array().sqrt();
   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
   MatrixXd jacobian(3, 4);
   double px = x_state(0);
   double py = x_state(1);
   double vx = x_state(2);
   double vy = x_state(3);

   double c1 = (px * px) + (py * py);
   double c2 = sqrt(c1);
   double c3 = (c1 * c2);

   if (fabs(c1) < 0.00000001)
   {
      std::cerr << "Error - Division by Zero!" << std::endl;
      return jacobian;
   }

   // compute the Jacobian matrix
   jacobian << (px / c2), (py / c2), 0, 0,
       -(py / c1), (px / c1), 0, 0,
       py * (vx * py - vy * px) / c3, px * (vy * px - vx * py) / c3, px / c2, py / c2;

   return jacobian;
}