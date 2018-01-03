#ifdef _WIN32
#undef __STRICT_ANSI__
#endif

#include "kalmanfilter.h"
#include "matrix.h"

KalmanFilter::KalmanFilter(
  const boost::numeric::ublas::matrix<double>& control,
  const boost::numeric::ublas::vector<double>& first_x,
  const boost::numeric::ublas::matrix<double>& first_p,
  const boost::numeric::ublas::matrix<double>& measurement_noise,
  const boost::numeric::ublas::matrix<double>& observation,
  const boost::numeric::ublas::matrix<double>& process_noise_covariance,
  const boost::numeric::ublas::matrix<double>& state_transition)
  : m_control(control),
    m_measurement_noise(measurement_noise),
    m_observation(observation),
    m_p(first_p),
    m_process_noise_covariance(process_noise_covariance),
    m_state_transition(state_transition),
    m_x(first_x)
{
  #ifndef NDEBUG
  //Check for correct dimensionality
  const auto sz = m_x.size();
  assert(m_control.size1() == sz);
  assert(m_control.size2() == sz);
  assert(m_measurement_noise.size1() == sz);
  assert(m_measurement_noise.size2() == sz);
  assert(m_observation.size1() == sz);
  assert(m_observation.size2() == sz);
  assert(m_p.size1() == sz);
  assert(m_p.size2() == sz);
  assert(m_process_noise_covariance.size1() == sz);
  assert(m_process_noise_covariance.size2() == sz);
  assert(m_state_transition.size1() == sz);
  assert(m_state_transition.size2() == sz);
  assert(m_x.size() == sz);
  #endif
}

void KalmanFilter::SupplyMeasurement(
  const boost::numeric::ublas::vector<double>& x)
{
  const boost::numeric::ublas::vector<double> input
    = boost::numeric::ublas::vector<double>(x.size(),0.0);
  return SupplyMeasurementAndInput(x,input);
}

void KalmanFilter::SupplyMeasurementAndInput(
  const boost::numeric::ublas::vector<double>& x,
  const boost::numeric::ublas::vector<double>& input)
{
  using boost::numeric::ublas::identity_matrix;
  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::prod;
  using boost::numeric::ublas::trans;
  using boost::numeric::ublas::vector;
  ///Debug statements to keep code below clean
  /// 2/7) Covariance prediction
  assert(m_state_transition.size2() == m_p.size1());
  assert( prod(m_state_transition,m_p).size2()
    ==  trans(m_state_transition).size1() );
  assert(matrix<double>(prod(
      matrix<double>(prod(m_state_transition,m_p)),
      trans(m_state_transition))).size1()
    == m_process_noise_covariance.size1());
  assert(matrix<double>(prod(
      matrix<double>(prod(m_state_transition,m_p)),
      trans(m_state_transition))).size2()
    == m_process_noise_covariance.size2());

  /// 1/7) State prediction
  const vector<double> x_current
    = prod(m_state_transition,m_x)
    + prod(m_control,input);
  /// 2/7) Covariance prediction
  const matrix<double> p_current
    = matrix<double>(
      prod(
        matrix<double>(prod(m_state_transition,m_p)),
        trans(m_state_transition)
      )
    )
    + m_process_noise_covariance;
  /// 3/7) Innovation (y with a squiggle above it)
  const vector<double> z_measured = x; //x has noise itn it
  const vector<double> innovation = z_measured - prod(m_observation,x_current);
  /// 4/7) Innovation covariance (S)
  const matrix<double> innovation_covariance
    = matrix<double>(prod(
          matrix<double>(prod(m_observation,p_current)),
          trans(m_observation)
        ))
    + m_measurement_noise;
  /// 5/7) Kalman gain (K)
  if (Matrix::CalcDeterminant(innovation_covariance) == 0.0)
  {
    throw std::runtime_error("Innovation covariance became degenerate");
  }
  const matrix<double> kalman_gain
    = prod(
        matrix<double>(
        prod(
          p_current,
          trans(m_observation)
        )),
        Matrix::Inverse(innovation_covariance)
      );
  /// 6/7) Update state prediction
  m_x = x_current + prod(kalman_gain,innovation);
  /// 7/7) Update covariance prediction
  const identity_matrix<double> my_identity_matrix(kalman_gain.size1());
  m_p = prod(
    my_identity_matrix
      - matrix<double>(prod(kalman_gain,m_observation)),
    p_current
  );
}

