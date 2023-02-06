QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        main.cpp

INCLUDEPATH += /home/user/PARALELA/PADRE/classification \ # Aqui es el path del proyecto con las cabeceras


LIBS += -L/home/user/PARALELA/PADRE/build \  # Directorio con los archivos compilados
    -lclassification \
    -ldcmdata \
    -ldcmimgle \
    -ldcmimage \
    -ldcmjpeg \
    -lpthread




unix {
    target.path = $$[QT_INSTALL_PLUGINS]/generic
}
!isEmpty(target.path): INSTALLS += target

