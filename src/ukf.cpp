#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // seperate initialization step for first measurement
  is_initialized_ = false;

  // for calculation of delta_t
  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // state dimensions
  n_x_ = 5;
  n_aug_ = 7;
  n_z_laser_ = 2;
  n_z_radar_ = 3;
  lambda_ = 3 - n_aug_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5; // Half of the maximum accelaration of a bike.

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.05;

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

  // initial predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // set weights for compensating lambda
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  // set laser measurement matrix
  H_ = MatrixXd(n_z_laser_, n_x_);
  H_.fill(0.0);
  H_(0, 0) = 1;
  H_(1, 1) = 1;

  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ << std_laspx_ *std_laspx_, 0, 0, std_laspy_ *std_laspy_;
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_ *std_radr_, 0, 0, 0, std_radphi_ *std_radphi_, 0, 0, 0,
      std_radrd_ *std_radrd_;
}

/**
TODO:

Complete the initialization. See ukf.h for other member properties.

Hint: one or more values initialized above might be wildly off...
*/

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    P_.fill(0.0);
    P_(0, 0) = 0.1;
    P_(1, 1) = 0.1;
    P_(2, 2) = 1;
    P_(3, 3) = 1;
    P_(4, 4) = 1;

    x_.fill(0.0);

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho, phi, rho_dot;
      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      rho_dot = meas_package.raw_measurements_[2];
      x_ << rho *cos(phi), rho * sin(phi), 0, 0, 0;
    }

    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float delta_t = (meas_package.timestamp_ - previous_timestamp_) /
                  1000000.0; // dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(delta_t);
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  createSigmaPoints(&Xsig_aug);
  predictSigmaPoints(Xsig_aug, delta_t);
  predictMeanCov();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_laser_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  x_ += K * y;
  P_ = (MatrixXd::Identity(n_x_, n_x_) - K * H_) * P_;

  NIS_laser_ = y.transpose()*S.inverse()*y;
  cout << "[DEBUG]: NIS_laser " << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  You'll also need to calculate the radar NIS.
  */
  MatrixXd Zsig_pred(n_z_radar_, 2 * n_aug_ + 1);
  Zsig_pred.fill(0.0);
  VectorXd z_pred(n_z_radar_);
  z_pred.fill(0.0);
  MatrixXd S_pred(n_z_radar_, n_z_radar_);
  S_pred.fill(0.0);
  VectorXd z = meas_package.raw_measurements_;

  predictMeasurementRadar(&Zsig_pred, &z_pred, &S_pred);
  innovateRadar(z, Zsig_pred, z_pred, S_pred);
  VectorXd y = z -z_pred;
  NIS_radar_ = y.transpose()*S_pred.inverse()*y;
  cout << "[DEBUG]: NIS_radar " << NIS_radar_ << endl;
}

void UKF::createSigmaPoints(MatrixXd *Xsig_aug) {

  // create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);

  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug->col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug->col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug->col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}

void UKF::predictSigmaPoints(MatrixXd Xsig_aug, double delta_t) {
  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yaw_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // avoid division by zero
    if (fabs(yaw_dot) < 0.0001) {
      Xsig_pred_(0, i) = px + v * cos(yaw) * delta_t +
                         0.5 * delta_t * delta_t * cos(yaw) * nu_a;
      Xsig_pred_(1, i) = py + v * sin(yaw) * delta_t +
                         0.5 * delta_t * delta_t * sin(yaw) * nu_a;
    } else {
      Xsig_pred_(0, i) =
          px + v / yaw_dot * (sin(yaw + yaw_dot * delta_t) - sin(yaw)) +
          0.5 * delta_t * delta_t * cos(yaw) * nu_a;
      Xsig_pred_(1, i) =
          py + v / yaw_dot * (-cos(yaw + yaw_dot * delta_t) + cos(yaw)) +
          0.5 * delta_t * delta_t * sin(yaw) * nu_a;
    }
    Xsig_pred_(2, i) = v + delta_t * nu_a;
    Xsig_pred_(3, i) =
        yaw + yaw_dot * delta_t + 0.5 * delta_t * delta_t * nu_yawdd;
    Xsig_pred_(4, i) = yaw_dot + delta_t * nu_yawdd;
  }
}

void UKF::predictMeanCov() {

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { // iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { // iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) {
      x_diff(3) -= 2. * M_PI;
      /*if (fabs(x_diff(3)) > 100) {
        x_diff(3) = 0;
      }*/
    }
    while (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
      /*if (fabs(x_diff(3)) > 100) {
        x_diff(3) = 0;
      }*/
    }
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::predictMeasurementRadar(MatrixXd *Zsig_pred, VectorXd *z_pred,
                                  MatrixXd *S_pred) {

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    // double yaw_dot = Xsig_pred_(4, i);

    double c1 = px * px + py * py;
    // avoid zero division
    if (fabs(c1) < 0.0001) {
      c1 = 0.0001;
    }
    (*Zsig_pred)(0, i) = sqrt(c1);
    (*Zsig_pred)(1, i) = atan2(py, px);
    (*Zsig_pred)(2, i) =
        (px * cos(yaw) * v + py * sin(yaw) * v) / (*Zsig_pred)(0, i);
  }

  // calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    *z_pred += weights_(i) * Zsig_pred->col(i);
  }

  // calculate measurement covariance matrix S
  MatrixXd Zsig_diff;
  Zsig_diff.fill(0.0);

  Zsig_diff = Zsig_pred->colwise() - *z_pred;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    while (Zsig_diff.col(i)(1) > M_PI)
      Zsig_diff.col(i)(1) -= 2. * M_PI;
    while (Zsig_diff.col(i)(1) < -M_PI)
      Zsig_diff.col(i)(1) += 2. * M_PI;
    *S_pred += weights_(i) * Zsig_diff.col(i) * Zsig_diff.col(i).transpose();
  }
  *S_pred += R_radar_;
}
void UKF::innovateRadar(VectorXd raw_measurements, MatrixXd Zsig_pred,
                        VectorXd z_pred, MatrixXd S_pred) {
  // calculate cross correlation matrix
  MatrixXd Tc(n_x_, n_z_radar_);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { // 2n+1 sigma points

    // residual
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI) {
      z_diff(1) -= 2. * M_PI;
      /*if (fabs(z_diff(1)) > 100) {
        z_diff(1) = 0;
      }*/
    }
    while (z_diff(1) < -M_PI) {
      z_diff(1) += 2. * M_PI;
      /*if (fabs(z_diff(3)) > 100) {
        z_diff(1) = 0;
      }*/
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) {
      x_diff(3) -= 2. * M_PI;
      /*if (fabs(x_diff(3)) > 100) {
        x_diff(3) = 0;
      }*/
    }
    while (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
      /*if (fabs(x_diff(3)) > 100) {
        x_diff(3) = 0;
      }*/
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S_pred.inverse();

  // residual
  VectorXd z = raw_measurements;
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI)
    z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI)
    z_diff(1) += 2. * M_PI;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += K * S_pred * K.transpose();
}
