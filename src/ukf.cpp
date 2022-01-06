#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  
  // time when the state is true, in us
  time_us_ = 0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3.0 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  for (int i = 0; i< 2*n_aug_+1;i++){
    if (i == 0)
      weights_(0) = lambda_/(lambda_ + n_aug_);
    else
      weights_(i) = 0.5/(lambda_ + n_aug_);
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    
    if (time_us_ == 0){
      time_us_ = meas_package.timestamp_;
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }
    
    else{
      UKF::Prediction(meas_package.timestamp_ - time_us_);
      time_us_ = meas_package.timestamp_;
      Zsig_ = MatrixXd(2, 2*n_aug_+1);
      S_ = MatrixXd(2,2);
      S_.fill(0.0);

      MatrixXd N(2,2);
      N << std_laspx_*std_laspx_, 0,
           0, std_laspy_*std_laspy_;

      z_pred_ = VectorXd(2);
      z_pred_ << VectorXd::Zero(2);

      for(int i = 0; i < 2*n_aug_; i++){
        double px = (Xsig_pred_(0,i));
        double py = (Xsig_pred_(1,i));
      
        Zsig_.col(i) << px, py;
      }
      for(int i = 0; i < 2*n_aug_+1; i++){
        z_pred_ = z_pred_ + weights_(i)*Zsig_.col(i);
      }
    
      for(int i = 0; i < 2*n_aug_+1; i++){
        S_ = S_ + weights_(i)*(Zsig_.col(i) - z_pred_)*(Zsig_.col(i) - z_pred_).transpose();
      }
      UKF::UpdateLidar(meas_package);
    }
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    
    if (time_us_ == 0){
      time_us_ = meas_package.timestamp_;
      double r = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rd = meas_package.raw_measurements_(2);

      double x = r*cos(phi);
      double y = r*sin(phi);
      double vx = rd*cos(phi);
      double vy = rd*sin(phi);
      double v = rd;
      x_ << x, y, v, r, rd;
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 0, std_radrd_*std_radrd_, 0, 0,
            0, 0, 0, std_radphi_*std_radphi_, 0,
            0, 0, 0, 0, std_radphi_*std_radphi_;
    }
    
    else{
      UKF::Prediction(meas_package.timestamp_ - time_us_);
      time_us_ = meas_package.timestamp_;
     
      Zsig_ = MatrixXd(3, 2*n_aug_+1);
      S_ = MatrixXd(3,3);
      S_.fill(0.0);

      MatrixXd N(3,3);
      N << std_radr_*std_radr_, 0, 0,
         0, std_radphi_*std_radphi_, 0,
         0, 0, std_radrd_*std_radrd_;

      z_pred_ = VectorXd(3);
      z_pred_ << VectorXd::Zero(3);

      for(int i = 0; i < 2*n_aug_; i++){
        double px = (Xsig_pred_(0,i));
        double py = (Xsig_pred_(1,i));
        double v  = (Xsig_pred_(2,i));
        double an = (Xsig_pred_(3,i));
        double vx = v*cos(an);
        double vy = v*sin(an);
        double a = sqrt(px*px + py*py);
        double b = atan2(py,px); 
        double c = (px*vx + py*vy)/a; 
      
        Zsig_.col(i) << a, b, c;
      }
    
      for(int i = 0; i < 2*n_aug_+1; i++){
        z_pred_ = z_pred_ + weights_(i)*Zsig_.col(i);
      }
    
      for(int i = 0; i < 2*n_aug_+1; i++){
        VectorXd zdiff = Zsig_.col(i) - z_pred_;
        while(zdiff(1) > M_PI) zdiff(1) -= 2.*M_PI;
        while(zdiff(1) <-M_PI) zdiff(1) += 2.*M_PI;
        S_ = S_ + weights_(i)*zdiff*zdiff.transpose();
      }
      UKF::UpdateRadar(meas_package);
    }
  
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  //std::cout << "Prediction state vector 1: " << x_(0) << std::endl;
  
  delta_t = delta_t/1e6;

  VectorXd x_aug_(n_aug_);
  MatrixXd P_aug_(n_aug_,n_aug_);
  MatrixXd Xsig_aug_(n_aug_, 2*n_aug_+1);
  
  x_aug_.head(n_x_) << x_;
  x_aug_.tail(n_aug_ - n_x_) << 0, 0;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) << P_;
  P_aug_.topRightCorner(n_x_, n_aug_ - n_x_) << MatrixXd::Zero(5,2);
  P_aug_.bottomLeftCorner(n_aug_ - n_x_, n_x_) << MatrixXd::Zero(2,5);
  P_aug_.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) << std_a_*std_a_, 0,
                                                            0, std_yawdd_*std_yawdd_;
  
  MatrixXd A = P_aug_.llt().matrixL();
  MatrixXd B = pow(lambda_+n_aug_, 0.5)*A;

  Xsig_aug_.col(0) = x_aug_;
  
  for (int i = 1; i < 2*n_aug_+1; i++){
    if (i <= n_aug_)
      Xsig_aug_.col(i) = x_aug_ + B.col(i - 1);
    
    else
      Xsig_aug_.col(i) = x_aug_ - B.col(i - 1 - n_aug_);
  }

  for(int i = 0; i < 2*n_aug_+1; i++){
    double a = Xsig_aug_(0,i);
    double b = Xsig_aug_(1,i);
    double c = Xsig_aug_(2,i);
    double d = Xsig_aug_(3,i);
    double e = Xsig_aug_(4,i);
    double f = Xsig_aug_(5,i);
    double g = Xsig_aug_(6,i);

    if (fabs(e) > 1e-10){
        Xsig_pred_.col(i) << c/e*(sin(d+e*delta_t)-sin(d)) + .5*delta_t*delta_t*cos(d)*f,
                             c/e*(-cos(d+e*delta_t)+cos(d)) + .5*delta_t*delta_t*sin(d)*f,
                             delta_t*f,
                             e*delta_t + .5*delta_t*delta_t*g,
                             delta_t*g;
    }
    else{
        Xsig_pred_.col(i) << c*cos(d) + .5*delta_t*delta_t*cos(d)*f,
                             c*sin(d) + .5*delta_t*delta_t*sin(d)*f,
                             delta_t*f,
                             e*delta_t + .5*delta_t*delta_t*g,
                             delta_t*g;
    }

    Xsig_pred_.col(i) = Xsig_pred_.col(i) + Xsig_aug_.col(i).head(n_x_);
  }

  x_ << VectorXd::Zero(n_x_);
  for (int i = 0; i < 2*n_aug_ + 1 ; i++){
      x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }

  P_ << MatrixXd::Zero(n_x_, n_x_);

  for(int i = 0; i < 2*n_aug_+1; i++){

        VectorXd xdiff = Xsig_pred_.col(i)-x_;
        while(xdiff(3) > M_PI) xdiff(3) -= 2.*M_PI;
        while(xdiff(3) <-M_PI) xdiff(3) += 2.*M_PI;
        P_ = P_ + weights_(i)*xdiff*xdiff.transpose();
  }
 //std::cout << "Prediction state vector 2: " << x_(0) << std::endl; 
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
 // std::cout << "Lidar state vector 1: " << x_(0) << std::endl;
  
  MatrixXd T(n_x_, 2);
  T.fill(0.0);

  for (int i = 0; i < 2*n_aug_+1; i++){

      VectorXd xdiff = Xsig_pred_.col(i)-x_;
      while(xdiff(3) > M_PI) xdiff(3) -= 2.*M_PI;
      while(xdiff(3) <-M_PI) xdiff(3) += 2.*M_PI;

      T = T + weights_(i)*xdiff*(Zsig_.col(i) - z_pred_).transpose();
  }
  MatrixXd K = T*S_.inverse();

  x_ = x_ + K*(meas_package.raw_measurements_ - z_pred_);
  P_ = P_ - K*S_*K.transpose();
  //std::cout << "Lidar state vector 2: " << x_(0) << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  //std::cout << "Radar state vector 1: " << x_(0) << std::endl;
  
  MatrixXd T(n_x_, 3);
  T.fill(0.0);

  for (int i = 0; i < 2*n_aug_+1; i++){

      VectorXd xdiff = Xsig_pred_.col(i)-x_;
      while(xdiff(3) > M_PI) xdiff(3) -= 2.*M_PI;
      while(xdiff(3) <-M_PI) xdiff(3) += 2.*M_PI;
      
      VectorXd zdiff = Zsig_.col(i) - z_pred_;
      while(zdiff(1) > M_PI) zdiff(1) -= 2.*M_PI;
      while(zdiff(1) <-M_PI) zdiff(1) += 2.*M_PI;

      T = T + weights_(i)*xdiff*zdiff.transpose();
  }

  MatrixXd K = T*S_.inverse();

  x_ = x_ + K*(meas_package.raw_measurements_ - z_pred_);
  P_ = P_ - K*S_*K.transpose();
  //std::cout << "Radar state vector 2: " << x_(0) << std::endl;
}