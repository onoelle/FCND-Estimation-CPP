#  Writeup for the FCND-Estimation-CPP project

fixed output window size, using the advice from https://knowledge.udacity.com/questions/292314

# Step 1: Sensor Noise (Scenario 6)

The task is to process the logged files from the simulation to calculate the standard deviation of the GPS X signal and the IMU Accelerometer X signal.

To accomplish this I wrote a little python script parsing and evaluating the log files:

```
import numpy as np

# calculate standard deviation from a series of measurements
# 1. gps x coordinate
gps_x_data = np.genfromtxt('Graph1.txt', delimiter=',', skip_header=1)
x = gps_x_data[:,1]
stdDev = np.std(x, axis = 0)
print(f'StdDev GPS X Position: {stdDev}')

# 2. accelerometer x acceleration
acc_ax_data = np.genfromtxt('Graph2.txt', delimiter=',', skip_header=1)
ax = acc_ax_data[:,1]
stdDev = np.std(ax, axis = 0)
print(f'StdDev Accelerometer X Acceleration: {stdDev}')
```
The result for a simulation of more than 10 seconds was:
* StdDev GPS X Position: 0.7017824842601307 (SimulatedSensors.txt: 0.7)
* StdDev Accelerometer ax: 0.5093821200191182 (SimulatedSensors.txt: 0.5) 

With these values inserted, the dashed lines turned green to signal that approx. 68% of the measurements are within the interval defined by the standard deviation:

![Scenario 6](./scenario6_stddev.png)

# Step 2: Attitude Estimation (Scenario 7)

The task is to improve estimation with a better rate gyro attitude integration scheme in the complementary filter.

To do this, I implemented an integration based on a quaternion representation in the inertial frame.
Then, the body rate from the gyro can be integrated in that representation to get the updated pitch, roll, and yaw values.
The changed implementation part within QuadEstimatorEKF.cpp looks like this:
```
Quaternion<float> attitude = Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, ekfState(6));
attitude.IntegrateBodyRate(gyro, dtIMU);

float predictedPitch = attitude.Pitch();
float predictedRoll = attitude.Roll();
ekfState(6) = attitude.Yaw();
```

With this improved attitude estimation, the attitude error is within 0.1 rad for more than 3 seconds:

![Scenario 7](./scenario7_att_est.png)

# Step 3: Prediction Step (Scenarios 8-?)

The task is to implement the prediction step of the Kalman filter.

## Subtask: transition function, state prediction (Scenario 8)
To accomplish this, the transition function in PredictState() calculates the pose and velocities using a simple integration, based on the assumption that dt is small:
```
////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  predictedState(0) = predictedState(0) + curState(3) * dt;
  predictedState(1) = predictedState(1) + curState(4) * dt;
  predictedState(2) = predictedState(2) + curState(5) * dt;

  // calculate global acceleration in order to predict velocities
  V3F accel_global = attitude.Rotate_BtoI(accel);
  accel_global.z -= static_cast<float>(CONST_GRAVITY);

  predictedState(3) = predictedState(3) + accel_global.x * dt;
  predictedState(4) = predictedState(4) + accel_global.y * dt;
  predictedState(5) = predictedState(5) - accel_global.z * dt;

  // yaw integral already done in the IMU update => no update to yaw here
```

With this state prediction, the estimator keeps track of the current state with only a reasonable drift:

![Scenario 8](./scenario8_predict.png)

## Subtask: Covariance prediction (Scenario 9)
To calculate the partial derivative of the body-to-global rotation matrix in the function GetRbgPrime(), I implemented equation (52) from the estimation paper:
```
////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  // implements equation (52) from the estimation paper (theta = pitch, phi = roll, psi = yaw)
  RbgPrime(0, 0) = -cos(pitch) * sin(yaw);
  RbgPrime(0, 1) = -sin(roll) * sin(pitch) * sin(yaw) - cos(roll) * cos(yaw);
  RbgPrime(0, 2) = -cos(roll) * sin(pitch) * sin(yaw) + sin(roll) * cos(yaw);

  RbgPrime(1, 0) = cos(pitch) * cos(yaw);
  RbgPrime(1, 1) = sin(roll) * sin(pitch) * cos(yaw) - cos(roll) * sin(yaw);
  RbgPrime(1, 2) = cos(roll) * sin(pitch) * cos(yaw) + sin(roll) * sin(yaw);
/////////////////////////////// END STUDENT CODE ////////////////////////////
```

With that function implemented, I implemented the Predict() function following equation(51) to predict the state covariance forward:

```
////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  // equation (51): calculate the derivative of transition function: G_t = g'(u_t, x_t, \Delta t) -> eq. 51
  gPrime(0,3) = dt;
  gPrime(1,4) = dt;
  gPrime(2,5) = dt;

  gPrime(3, 6) = (RbgPrime(0) * accel).sum() * dt;
  gPrime(4, 6) = (RbgPrime(1) * accel).sum() * dt;
  gPrime(5, 6) = (RbgPrime(2) * accel).sum() * dt;

  // calculate the covariance: \bar{\Sigma}_t = G_t\Sigma_{t-1}G_t^T + Q_t
  MatrixXf term_1 = gPrime * ekfCov;
  gPrime.transposeInPlace();
  ekfCov = term_1 * gPrime + Q;

/////////////////////////////// END STUDENT CODE ////////////////////////////
```

Tuning the parameters to these values:
````
QPosXYStd = .0005
QVelXYStd = .3
````
results in a much better approximation of the errors regarding x and velocity of x, judging from these sample plots: 
![Scenario 9](./scenario9_covariance.png)

# Step 4: Magnetometer Update
## Tune the parameter QYawStd
In my simulation setting QYawStd to 0.01 resulted in this behavior, where the std deviation captures 81% of the values.

![Scenario 10](./scenario10_qyawstd_tuning.png)

However, I could not reproduce the drift of the yaw error over time as described in the project readme.
I'm suspecting a tuning issue, so after implementing the UpdateFromMag() method, I performed a lot of retuning of my so-far used parameters to pass this scenario successfully:

![Scenario 10](./scenario10_updateMag.png)

The resulting tuning parameters were:
```
QPosXYStd = .001
QPosZStd = .05
QVelXYStd = .3
QVelZStd = .1
QYawStd = .09

# GPS measurement std deviations
GPSPosXYStd = 1
# was: 3
GPSPosZStd = 300
GPSVelXYStd = .1
GPSVelZStd = .3
```

The implementation of the UpdateFromMag() method is
```
////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  // section 7.3.2 estimation paper, equations 56 - 58
  // equation (58)
  hPrime(0, 6) = 1.;

  zFromX(0) = ekfState(6);
  float diffYaw = magYaw - ekfState(6);
  if ( diffYaw > F_PI ) {
      zFromX(0) += 2.f * F_PI;
  } else if ( diffYaw < -F_PI ) {
      zFromX(0) -= 2.f * F_PI;
  }
  
/////////////////////////////// END STUDENT CODE ////////////////////////////

```

# Step 5: Closed Loop + GPS Update
