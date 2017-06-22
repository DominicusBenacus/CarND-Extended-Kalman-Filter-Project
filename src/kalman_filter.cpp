#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  TODO:
    * predict the state
  */
  // There is no external motion, so, we do not have to add "+u"
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	VectorXd y = z - H_ * x_;

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;

	// New state
	int x_size = x_.size();
	x_ = x_ + (K * y);
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	// Calculation of the inovation matrix. the correction of the system state
	// depends the result in y. In collelation with the Kalman gain it will be meaningful.
	VectorXd y = z - linearizeH(x_);
	
    // update the state by using Extended Kalman Filter equations
    // In C++, atan2() returns values between -pi and pi. When calculating phi in y = z - h(x)
    // for radar measurements, the resulting angle phi in the y vector should be adjusted so
    // that it is between -pi and pi. The Kalman filter is expecting small angle values between
    // the range -pi and pi. HINT: when working in radians, you can add 2? or subtract 2? until
    // the angle is within the desired range.
	y[1] -= (2 * M_PI) * floor((y[1] + M_PI) / (2 * M_PI));

	// Transpose the H matrix. Keep in mind the H matriy is replaced
	// by Jacobian Matrix before called in FusionEKF
	MatrixXd Ht = H_.transpose();
	// Calculatoe the absolute variance with bayes rule for further calculations
	MatrixXd S = H_ * P_ * Ht + R_;
	// calcluate the inverse of the absolute probability
	MatrixXd Si = S.inverse();
	// Define the Kalman gain to define in which values we trust more. Measurement or prediction
	MatrixXd K = P_ * Ht * Si;

	// Define new most probable state
	int x_size = x_.size();
	x_ = x_ + (K * y);
	// Define the new Covariance Matrix
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;


}
VectorXd KalmanFilter::linearizeH(const VectorXd &x) {
	// extract position and velocity
	float px = x(0);
	float py = x(1);
	float vx = x(2);
	float vy = x(3);
	float rho_dot;

	VectorXd linearizeH = VectorXd(3);

	float rho = sqrt(px*px + py*py);

	if (fabs(rho) < 0.0001) {
		rho_dot = 0;
	}
	else {
		rho_dot = (px*vx + py*vy) / rho;
	}

	float phi = atan2(py, px);

	linearizeH << rho, phi, rho_dot;


	return linearizeH;
}