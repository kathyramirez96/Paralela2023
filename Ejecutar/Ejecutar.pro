QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        main.cpp



win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../build-IMAGENES-Desktop-Debug/release/ -lIMAGENES
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../build-IMAGENES-Desktop-Debug/debug/ -lIMAGENES
else:unix: LIBS += -L$$PWD/../../build-IMAGENES-Desktop-Debug/ -lIMAGENES

INCLUDEPATH += $$PWD/../IMAGENES
DEPENDPATH += $$PWD/../IMAGENES

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../build-IMAGENES-Desktop-Debug/release/libIMAGENES.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../build-IMAGENES-Desktop-Debug/debug/libIMAGENES.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../build-IMAGENES-Desktop-Debug/release/IMAGENES.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../build-IMAGENES-Desktop-Debug/debug/IMAGENES.lib
else:unix: PRE_TARGETDEPS += $$PWD/../../build-IMAGENES-Desktop-Debug/libIMAGENES.a
