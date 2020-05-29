#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>
#include <typeinfo>

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  /**
   * TODO: predict the state
   */

  x_ = F_ * x_;
  auto F_T = F_.transpose();
  auto F_P = F_ * P_;
  auto F_P_F_T = F_P * F_T;
  P_ = F_P_F_T * F_T + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  /**
   * TODO: update the state by using Kalman Filter equations
   */

  VectorXd y = z - H_ * x_;
  MatrixXd H_t = H_.transpose();
  auto H_P = H_ * P_;
  auto H_P_H_t = H_P * H_t;
  MatrixXd S = H_P_H_t + R_;
  MatrixXd S_i = S.inverse();
  // measurement covariance
  auto P_H_t = P_ * H_t;
  MatrixXd K = P_H_t * S_i;
  auto K_y = K * y;
  x_ = x_ + K_y;
  auto x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  auto K_H = K * H_;
  P_ = (I - K_H) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  double px = x_(0);
  double py = x_(1);

  double vx = x_(2);
  double vy = x_(3);

  double rho_pred = sqrt(pow(px, 2) + pow(py, 2));
  double phi_pred = 0.0;

  if (fabs(px) > 0.001)
  {
    phi_pred = atan2(py, px);
  }

  double rhodot_pred = 0.0;
  if (fabs(rho_pred) > 0.00000001)
  {
    rhodot_pred = (px * vx + py * vy) / rho_pred;
  }

  VectorXd z_pred(3);
  z_pred << rho_pred, phi_pred, rhodot_pred;

  //Now apply the udate equations again

  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  auto H_P = H_ * P_;
  auto H_P_Ht = H_P * Ht;
  MatrixXd S = H_P_Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimates
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd ::Identity(x_size, x_size);
  auto K_H = K * H_;
  P_ = (I - K_H) * P_;
}
