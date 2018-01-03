#ifndef WHITENOISESYSTEM_H
#define WHITENOISESYSTEM_H

#include <boost/numeric/ublas/matrix.hpp>

struct WhiteNoiseSystem
{
  WhiteNoiseSystem(
    const boost::numeric::ublas::matrix<double>& control,
    const boost::numeric::ublas::vector<double>& initial_state,
    const boost::numeric::ublas::vector<double>& real_measurement_noise,
    const boost::numeric::ublas::vector<double>& real_process_noise,
    const boost::numeric::ublas::matrix<double>& state_transition);

  ///Update reality, that is, let the real system (i.e. reality) go to its next state,
  ///without any input
  void GoToNextState();

  ///Update reality, that is, let the real system (i.e. reality) go to its next state
  void GoToNextState(const boost::numeric::ublas::vector<double>& input);

  ///Measure a value from this system with normally distributed noise
  const boost::numeric::ublas::vector<double> Measure() const;

  ///Peek what the real value is
  const boost::numeric::ublas::vector<double>& PeekAtRealState() const { return m_current_state; }

  private:
  ///The control matrix to determine the influence of the input (in GoToNextState)
  const boost::numeric::ublas::matrix<double> m_control;

  ///The real value of the system
  boost::numeric::ublas::vector<double> m_current_state;

  ///The standard deviation of the noise in the state transition (in GoToNextState)
  const boost::numeric::ublas::vector<double> m_process_noise;

  ///The amount of noise in the system
  ///A noise of zero indicates a system that can be measured accurately to infinite precision
  const boost::numeric::ublas::vector<double> m_real_measurement_noise;

  ///The state transitions in the system, used in GoToNextState
  const boost::numeric::ublas::matrix<double> m_state_transition;

  ///Obtain a random number from a normal distribution
  ///From http://www.richelbilderbeek.nl/CppGetRandomNormal.htm
  static double GetRandomNormal(const double mean = 0.0, const double sigma = 1.0);

};

#endif // WHITENOISESYSTEM_H
