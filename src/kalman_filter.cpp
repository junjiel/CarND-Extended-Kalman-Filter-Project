#include "kalman_filter.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  std::cout << "x_ update starts " << std::endl;
  	float px = x_(0);
  	float py = x_(1);
  	float vx = x_(2);
    float vy = x_(3);
  	
  // pre-compute a set of terms to avoid repeated calculation
std::cout << "x_ conversion starts" << std::endl;
  	float rho = sqrt(px*px+py*py);
    float phi = atan2(py,px);
  	float rho_dot;
  	const float pi = 3.1415926;
  	//normalize phi
  std::cout << "normalization starts " << std::endl;
  	if (phi>pi){
		phi -=2*pi;  
  	}else if(phi < -pi){
    	phi+= 2*pi;
    }else{
    	phi = phi;
    }
  	//check if px*px+ py*py is close to 0
  std::cout << "check divide by 0 " << std::endl;
  	if (px*px+py*py < 0.0001){
    	std::cout<< "px**2 + py**2 is too close to 0, no valid EKF update" << std::endl;
      	rho_dot = (px*vx + py*vy)/(sqrt(px*px+py*py)+0.01);
    } else{
      	rho_dot = (px*vx + py*vy)/sqrt(px*px+py*py);
    }
    
  //map the predicted states x' from Cartesian coordinates to polar coordinates
  std::cout << "hx init " << std::endl;
  	VectorXd hx_(3) ;
  	hx_<< rho, phi, rho_dot;
  //EKF update
  std::cout << "EKF update " << std::endl;
  	VectorXd y = z - hx_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;
  
  	x_ = x_ + K*y;
    long x_size = x_.size();
  	MatrixXd I = MatrixXd::Identity(x_size, x_size);
  	P_ = (I - K*H_)*P_;
}
