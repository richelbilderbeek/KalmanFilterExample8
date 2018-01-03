#ifdef _WIN32
#undef __STRICT_ANSI__
#endif

///Kalman filter example
///Adapted from merge from www.adrianboeing.com and http://greg.czerniak.info/guides/kalman1
///following
/// * Simon, Dan. Kalman Filtering. Embedded Systems Programming. June 2001.

#include <iostream>
#include "maindialog.h"
#include "matrix.h"

///Context:
///A car that has a constant acceleration that has its position determined by GPS
///The car its speedometer is not used (as observation(1,1) is equal to 0.0),
///  and gives junk values (as x_real_noise(1) is equal to 10000000.0)
int main()
{
  const double t = 0.1;

  //The name of the states:
  //x: position
  //v: velocity
  const std::vector<std::string> state_names = { "x", "v" };

  //The real state vector
  //[ position ]
  //[ velocity ]
  const boost::numeric::ublas::vector<double> init_x_real = Matrix::CreateVector( { 0.0, 0.0 } );

  //Real measurement noise
  //[ standard deviation of noise in position ]   [ standard deviation of noise in GPS                       ]
  //[ standard deviation of noise in velocity ] = [ standard deviation of noise in defect/unused speedometer ]
  const boost::numeric::ublas::vector<double> x_real_measurement_noise = Matrix::CreateVector( { 10.0, 10000000.0 } );

  //Guess of the state matrix
  //Position and velocity guess is way off on purpose
  //[ position ]
  //[ velocity ]
  const boost::numeric::ublas::vector<double> x_first_guess = Matrix::CreateVector( { 100.0, 10.0 } );

  //Guess of the covariances
  //[ 1.0   0.0 ]
  //[ 0.0   1.0 ]
  const boost::numeric::ublas::matrix<double> p_first_guess = Matrix::CreateMatrix(2,2, { 1.0, 0.0, 0.0, 1.0 } );

  //Effect of inputs on state
  //Input = gas pedal, which gives acceleration 'a'
  //[ 1.0   0.5 * t * t ]   [teleportation (not used)                 x = 0.5 * a * t * t ]
  //[ 0.0   t           ] = [no effect of teleportation on velocity   v = a * t           ]
  const boost::numeric::ublas::matrix<double> control = Matrix::CreateMatrix(2,2, { 1.0, 0.0, 0.5 * t * t, t } );

  //Estimated measurement noise
  //[ 10.0          0.0 ]   [ Estimated noise in GPS   ?                                                     ]
  //[  0.0   10000000.0 ] = [ ?                        Estimated noise in speedometer (absent in this setup) ]
  const boost::numeric::ublas::matrix<double> measurement_noise = Matrix::CreateMatrix(2,2, { 10.0, 0.0, 0.0, 10000000.0 } );

  //Observational matrix
  //[ 1.0   0.0 ]   [GPS measurement   ?                                         ]
  //[ 0.0   0.0 ] = [?                 Speedometer (absent/unused in this setup) ]
  const boost::numeric::ublas::matrix<double> observation = Matrix::CreateMatrix(2,2, { 1.0, 0.0, 0.0, 0.0 } );

  //Real process noise
  //[ 0.001 ]   [ noise in position ]
  //[ 0.001 ] = [ noise in velocity ]
  const boost::numeric::ublas::vector<double> real_process_noise = Matrix::CreateVector( {0.01,0.01} );

  //Estimated process noise covariance
  //[ 0.01   0.01 ]
  //[ 0.01   0.01 ]
  const boost::numeric::ublas::matrix<double> process_noise = Matrix::CreateMatrix(2,2,{0.01,0.01,0.01,0.01});

  //State transition matrix, the effect of the current state on the next
  //[ 1.0     t ]   [ position keeps its value             a velocity changes the position ]
  //[ 0.0   1.0 ] = [ position has no effect on velocity   a velocity keeps its value      ]
  const boost::numeric::ublas::matrix<double> state_transition = Matrix::CreateMatrix(2,2,{1.0,0.0,t,1.0});

  //There is a constant push on the gas pedal. This has no direct effect on the position,
  //but it does increase velocity with accelation every state transition
  const double acceleration = 1.0;
  const boost::numeric::ublas::vector<double> input = Matrix::CreateVector( {0.0,acceleration} );

  const int time = 250;
  const MainDialog d(
    time,
    control,
    input,
    measurement_noise,
    observation,
    p_first_guess,
    process_noise,
    state_transition,
    init_x_real,
    real_process_noise,
    state_names,
    x_first_guess,
    x_real_measurement_noise);

  //Display header
  {
    const boost::numeric::ublas::vector<std::string> header = d.GetHeader(state_names);
    const std::size_t sz = header.size();
    for (std::size_t i=0; i!=sz; ++i)
    {
      std::cout << header(i);
      if (i != sz - 1) std::cout << ",";
    }
    std::cout << '\n';
  }
  //Display data
  {
    const boost::numeric::ublas::matrix<double>& data = d.GetData();
    const std::size_t n_cols = 6;
    assert(n_cols == d.GetHeader(state_names).size());
    for (int row=0; row!=time; ++row)
    {
      for (std::size_t col=0; col!=n_cols; ++col)
      {
        std::cout << data(row,col);
        if (col != 5) std::cout << ",";
      }
      std::cout << '\n';
    }
  }
}
