#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    // take length parameter (distance of rotors) into account
    float len = L / sqrtf(2.f);
    
    // total thrust
    float F_total = collThrustCmd;
    // calculate the three thrust components for x/y/z axis
    float F_roll_x = momentCmd.x / len; // moment on x axis/roll
    float F_pitch_y = momentCmd.y / len; // moment on y axis/pitch
    float F_yaw_z = - momentCmd.z / kappa; // moment on z axis/yaw

    // derive desired thrust for each rotor
    cmd.desiredThrustsN[0] = (F_total + F_roll_x + F_pitch_y + F_yaw_z)/4.f;  // front left / f1
    cmd.desiredThrustsN[1] = (F_total - F_roll_x + F_pitch_y - F_yaw_z)/4.f;  // front right / f2
    cmd.desiredThrustsN[2] = (F_total + F_roll_x - F_pitch_y - F_yaw_z)/4.f;  // rear left / f4
    cmd.desiredThrustsN[3] = (F_total - F_roll_x - F_pitch_y + F_yaw_z)/4.f;  // rear right / f3
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    // put parameters for moments of inertia into a vector
    V3F momentOfInertia(Ixx, Iyy, Izz);
    
    // calculate error between commanded/desired moment and actual one
    V3F bodyRateError = pqrCmd - pqr;
    
    // calculate vector of desired moments based on error and on inertia and gain parameters
    momentCmd = momentOfInertia * kpPQR * bodyRateError;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first


    
  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    // convert thrust force into acceleration (taking the mass into account)
    float accel = - collThrustCmd / mass;
    
    if ( collThrustCmd > 0 ) {
        // only required if total thrust is positive
        
        // limit roll command according to maxTiltAngle
        float body_roll_x_cmd = CONSTRAIN(accelCmd.x / accel, -maxTiltAngle, maxTiltAngle);
        
        // calculate error on x axis based on the given rotation matrix
        float body_roll_x_err = body_roll_x_cmd - R(0,2);
        // apply gain parameter kpBank
        float body_roll_x_p_term = kpBank * body_roll_x_err;
        
        // limit pitch command according to maxTiltAngle
        float body_pitch_y_cmd = CONSTRAIN(accelCmd.y / accel, -maxTiltAngle, maxTiltAngle);
        // calculate error on y axis based on the given rotation matrix
        float body_pitch_y_err = body_pitch_y_cmd - R(1,2);
        // apply gain parameter kpBank
        float body_pitch_y_p_term = kpBank * body_pitch_y_err;
        
        // derive rotational forces for x/y component of the pqr command
        pqrCmd.x = (R(1,0) * body_roll_x_p_term - R(0,0) * body_pitch_y_p_term) / R(2,2);
        pqrCmd.y = (R(1,1) * body_roll_x_p_term - R(0,1) * body_pitch_y_p_term) / R(2,2);
    } else {
        pqrCmd.x = 0.0;
        pqrCmd.y = 0.0;
    }
    
    // z component of pqr command is 0, as we only control roll and pitch here
    pqrCmd.z = 0;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    // calculate altitude error/error on z-axis
    float z_err = posZCmd - posZ;
    // calculate error in velocity in z direction
    float z_dot_err = velZCmd - velZ;
    // integrate the error over one time step
    integratedAltitudeError += z_err * dt;
    
    // calculate the components for a PID controller for altitude, using the corresponding parameters
    float p_term = kpPosZ * z_err;
    float d_term = kpVelZ * z_dot_err + velZ;
    float i_term = KiPosZ * integratedAltitudeError;
    float b_z = R(2,2);

    // determine u1 bar and acceleration
    float u_1_bar = p_term + d_term + i_term + accelZCmd;
    float accel = ( u_1_bar - CONST_GRAVITY ) / b_z;

    // derive collective thrust in [N], taking mass into accound and limiting to the range of [-maxDescentrate, maxAscentRate]
    thrust = -mass * CONSTRAIN(accel, -maxDescentRate / dt, maxAscentRate / dt);

  /////////////////////////////// END STUDENT CODE ////////////////////////////
  
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    // for easier handling: represent position and velocity parameters as vectors
    V3F kpPos(kpPosXY, kpPosXY, 0.f);
    V3F kpVel(kpVelXY, kpVelXY, 0.f);
    
    // limit velocity command to take limits of drone into account
    if ( velCmd.mag() > maxSpeedXY ) {
        velCmd = velCmd.norm() * maxSpeedXY;
    }
    
    // calculate position and velocity errors
    V3F posErr = posCmd - pos;
    V3F velErr = velCmd - vel;
    // derive which value to add to the feed-forward value
    // via proportionally taking both errors into account
    accelCmd += kpPos * posErr + kpVel * velErr;
    
    // limit resulting acceleration command to take physical imits of drone into account
    if ( accelCmd.mag() > maxAccelXY ) {
        accelCmd = accelCmd.norm() * maxAccelXY;
    }
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b]. 
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    // unwrap commanded yaw either to [-2Pi, 0] or [0, 2Pi]
    if (yawCmd > 0) {
        yawCmd = fmodf(yawCmd, 2.f * F_PI);
    } else {
        yawCmd = fmodf(yawCmd, -2.f * F_PI);
    }

    // normalize into range [-Pi, Pi]
    if (yawCmd <= -F_PI) {
        yawCmd += (2.0f * F_PI);
    } else {
        if (yawCmd > F_PI) {
            yawCmd -= (2.0f * F_PI);
        }
    }

    // calculate yaw error
    float yawErr = yawCmd - yaw;
    // calculate yaw rate proportionally on error, using kpYaw parameter
    yawRateCmd = kpYaw * yawErr;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
