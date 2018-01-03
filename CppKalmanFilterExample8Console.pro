TEMPLATE = app
CONFIG += console
CONFIG -= qt
QMAKE_CXXFLAGS += -Wall -Wextra -Werror -std=c++0x

win32 {
  INCLUDEPATH += E:/boost_1_50_0

  LIBS += \
    -LE:/boost_1_50_0/stage/lib  \
    -lboost_system-mgw44-mt-1_50 \
    -lboost_filesystem-mgw44-mt-1_50 \
    -lboost_regex-mgw44-mt-1_50
}


SOURCES += main.cpp \
    kalmanfilter.cpp \
    whitenoisesystem.cpp \
    matrix.cpp \
    maindialog.cpp

HEADERS += \
    kalmanfilter.h \
    whitenoisesystem.h \
    matrix.h \
    maindialog.h

