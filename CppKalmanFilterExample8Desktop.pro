QT       += core gui
TEMPLATE = app
QMAKE_CXXFLAGS += -Wall -Wextra -Werror -std=c++0x

#Prevents error:
#/my_boost_folder/boost/type_traits/detail/has_binary_operator.hp:50: Parse error at "BOOST_JOIN"
DEFINES += BOOST_TT_HAS_OPERATOR_HPP_INCLUDED

win32 {
  INCLUDEPATH += E:/boost_1_50_0 C:/qwt-6.0.1/include
  LIBS += \
    -LE:/boost_1_50_0/stage/lib  \
    -lboost_system-mgw44-mt-1_50 \
    -lboost_filesystem-mgw44-mt-1_50 \
    -lboost_regex-mgw44-mt-1_50 \
    -LC:/qwt-6.0.1/lib \
    -lqwtd  #Note: gives error 'QWidget: Must construct a QApplication before a QPaintDevice' when using '-lqwt'
    #-lQtSvg
}

unix {
  LIBS += -lqwt-qt4
  INCLUDEPATH += /usr/include/qwt-qt4/
}


SOURCES += \
    qtmain.cpp \
    qtmaindialog.cpp \
    kalmanfilter.cpp \
    whitenoisesystem.cpp \
    matrix.cpp \
    maindialog.cpp

HEADERS  += \
    qtmaindialog.h \
    whitenoisesystem.h \
    matrix.h \
    kalmanfilter.h \
    maindialog.h

FORMS    += qtmaindialog.ui
