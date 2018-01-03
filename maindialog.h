#ifndef MAINDIALOG_H
#define MAINDIALOG_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

struct MainDialog
{
  MainDialog(
    const int time,
    const boost::numeric::ublas::matrix<double>& control,
    const boost::numeric::ublas::vector<double>& input,
    const boost::numeric::ublas::matrix<double>& measurement_noise,
    const boost::numeric::ublas::matrix<double>& observation,
    const boost::numeric::ublas::matrix<double>& p_first_guess,
    const boost::numeric::ublas::matrix<double>& process_noise,
    const boost::numeric::ublas::matrix<double>& state_transition,
    const boost::numeric::ublas::vector<double>& init_x_real,
    const boost::numeric::ublas::vector<double>& real_process_noise,
    const std::vector<std::string>& state_names,
    const boost::numeric::ublas::vector<double>& x_first_guess,
    const boost::numeric::ublas::vector<double>& x_real_measurement_noise);


  ///Obtain the simulation data
  const boost::numeric::ublas::matrix<double>& GetData() const { return m_data; }

  ///Obtain the header text from the states
  static const boost::numeric::ublas::vector<std::string> GetHeader(const std::vector<std::string>& state_names);

  private:

  ///The simulation data
  const boost::numeric::ublas::matrix<double> m_data;

  //Real initial state
  const boost::numeric::ublas::vector<double> m_init_x_real;

  //Real measurement noise
  const boost::numeric::ublas::vector<double> m_x_real_measurement_noise;

  //Guess of the state matrix
  const boost::numeric::ublas::vector<double> m_x_first_guess;

  //Guess of the covariances
  const boost::numeric::ublas::matrix<double> m_p_first_guess;

  //Effect of inputs on state
  const boost::numeric::ublas::matrix<double> m_control;

  //Estimated measurement noise
  const boost::numeric::ublas::matrix<double> m_measurement_noise;

  //Observational matrix
  const boost::numeric::ublas::matrix<double> m_observation;

  //Real process noise
  const boost::numeric::ublas::vector<double> m_real_process_noise;

  //Estimated process noise covariance
  const boost::numeric::ublas::matrix<double> m_process_noise;

  //State transition matrix, the effect of the current state on the next
  const boost::numeric::ublas::matrix<double> m_state_transition;


  ///Create the simulation data
  static const boost::numeric::ublas::matrix<double> CreateData(
    const int time,
    const boost::numeric::ublas::matrix<double>& control,
    const boost::numeric::ublas::vector<double>& input,
    const boost::numeric::ublas::matrix<double>& measurement_noise,
    const boost::numeric::ublas::matrix<double>& observation,
    const boost::numeric::ublas::matrix<double>& p_first_guess,
    const boost::numeric::ublas::matrix<double>& process_noise,
    const boost::numeric::ublas::matrix<double>& state_transition,
    const boost::numeric::ublas::vector<double>& init_x_real,
    const boost::numeric::ublas::vector<double>& real_process_noise,
    const std::vector<std::string>& state_names,
    const boost::numeric::ublas::vector<double>& x_first_guess,
    const boost::numeric::ublas::vector<double>& x_real_measurement_noise);
};

#endif // MAINDIALOG_H
