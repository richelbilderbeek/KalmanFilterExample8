#ifdef _WIN32
#undef __STRICT_ANSI__
#endif

#include <QtGui/QApplication>
#include "qtmaindialog.h"

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  QtMainDialog w;
  w.show();
  
  return a.exec();
}
