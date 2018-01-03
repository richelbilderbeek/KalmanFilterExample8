#ifndef QTMAINDIALOG_H
#define QTMAINDIALOG_H

#ifdef _WIN32
#undef __STRICT_ANSI__
#endif

#include <vector>

#include <QDialog>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace Ui {
  class QtMainDialog;
}

struct QwtPlot;
struct QwtPlotCurve;
struct QTableWidget;

class QtMainDialog : public QDialog
{
  Q_OBJECT
  
public:
  explicit QtMainDialog(QWidget *parent = 0);
  ~QtMainDialog();
  
private:
  Ui::QtMainDialog *ui;

  bool m_can_do_sim;
  std::vector<QwtPlotCurve *> m_curves;
  std::vector<QwtPlot *> m_plots;

  static const boost::numeric::ublas::matrix<double> ToMatrix(const QTableWidget * const table);
  static const boost::numeric::ublas::vector<double> ToVector(const QTableWidget * const table);
  void UpdateLegends();

private slots:
  void OnAnyChange();
  void on_box_n_states_valueChanged(int arg1);
  const std::vector<std::string> GetLegend() const;
  const std::vector<QTableWidget *> CollectMatrices() const;
  const std::vector<QTableWidget *> CollectVectors() const;
};

#endif // QTMAINDIALOG_H
