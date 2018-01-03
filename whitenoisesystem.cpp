#ifdef _WIN32
#undef __STRICT_ANSI__
#endif

#include "whitenoisesystem.h"

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>

WhiteNoiseSystem::WhiteNoiseSystem(
  const boost::numeric::ublas::matrix<double>& control,
  const boost::numeric::ublas::vector<double>& initial_state,
  const boost::numeric::ublas::vector<double>& real_measurement_noise,
  const boost::numeric::ublas::vector<double>& real_process_noise,
  const boost::numeric::ublas::matrix<double>& state_transition)
  : m_control(control),
    m_current_state(initial_state),
    m_process_noise(real_process_noise),
    m_real_measurement_noise(real_measurement_noise),
    m_state_transition(state_transition)
{
  #ifndef NDEBUG
  //Check for correct dimensionality
  const auto sz = initial_state.size();
  assert(m_control.size1() == sz);
  assert(m_control.size2() == sz);
  assert(m_current_state.size() == sz);
  assert(m_process_noise.size() == sz);
  assert(m_real_measurement_noise.size() == sz);
  assert(m_state_transition.size1() == sz);
  assert(m_state_transition.size2() == sz);
  #endif
}

double WhiteNoiseSystem::GetRandomNormal(const double mean, const double sigma)
{
  boost::normal_distribution<double> norm_dist(mean, sigma);
  static boost::lagged_fibonacci19937 engine;
  const double value = norm_dist.operator () <boost::lagged_fibonacci19937>((engine));
  return value;
}

void WhiteNoiseSystem::GoToNextState()
{
  //Create a no-input vector
  const auto sz = m_current_state.size();
  boost::numeric::ublas::vector<double> input(sz,0.0);
  return GoToNextState(input);
}

void WhiteNoiseSystem::GoToNextState(const boost::numeric::ublas::vector<double>& input)
{
  //First do a perfect transition
  m_current_state
    = boost::numeric::ublas::prod(m_state_transition,m_current_state)
    + boost::numeric::ublas::prod(m_control,input);
  //Add process noise
  const auto sz = m_current_state.size();
  assert(m_current_state.size() == m_process_noise.size());
  for (std::size_t i = 0; i!=sz; ++i)
  {
    m_current_state(i) = GetRandomNormal(m_current_state(i),m_process_noise(i));
  }
}

const boost::numeric::ublas::vector<double> WhiteNoiseSystem::Measure() const
{
  const auto sz = m_current_state.size();
  assert(m_current_state.size() == m_real_measurement_noise.size());
  boost::numeric::ublas::vector<double> measured(sz);
  for (std::size_t i = 0; i!=sz; ++i)
  {
    measured(i) = GetRandomNormal(m_current_state(i),m_real_measurement_noise(i));
  }
  return measured;
}
