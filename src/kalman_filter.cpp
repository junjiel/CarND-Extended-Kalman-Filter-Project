#include "kalman_filter.h"

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
  	float px = x_(0);
  	float py = x_(1);
  	float vx = x_(2);
    float vy = x_(3);
  	
  // pre-compute a set of terms to avoid repeated calculation

  	float rho = sqrt(px*px+py*py);
    float phi = atan2(py,px);
  	float rho_dot;
  	//normalize phi
  	if (phi>pi){
		phi -=2*pi;  
  	}else if(phi < -pi){
    	phi+= 2*pi;
    }else{
    	phi = phi;
    }
  	//check if px*px+ py*py is close to 0
  	if (px*px+py*py < 0.0001){
    	cout<< "px**2 + py**2 is too close to 0, no valid EKF update" << endl;
      	rho_dot = (px*vx + py*vy)/(sqrt(px*px+py*py)+0.01);
    } else{
      	rho_dot = (px*vx + py*vy)/sqrt(px*px+py*py);
    }
    
  //map the predicted states x' from Cartesian coordinates to polar coordinates
  	VectorXd hx_ << rho, phi, rho_dot;
  	VectorXd y = z - hx_;
  	MatrixXd Hj = CalculateJacobian(x_);
  	MatrixXd Hjt = Hj_.transpose();
    MatrixXd S = Hj * P_ * Hjt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Hjt * Si;
  
  	x_ = x_ + K*y;
    long x_size = x_.size();
  	MatrixXd I = MatrixXd::Identity(x_size, x_size);
  	P_ = (I - K*Hj)*P_;
}
