# Dicom Classifier

Dicom Classifier es una libreria para manejar imagenes médicas dicom y con diferentes técnicas de clasificación.

## Dependencias
La libreria usa:
- **libdcmtk-dev**: para el manejo de imágenes dicom.

Para instalarlos puedes ejecutar el script `checker.sh`

```bash
chmod +x ./checker.sh
./checker.sh
```

Tambien se pueden instalar de manera manual:
```bash
sudo apt install libdcmtk-dev build-essential cmake
```

## Instalación

Antes de usar la libreria se necesita realizar el proceso de compilación para generar los archivos necesarios e importarlos directamente en su proyecto.

### QT

Para compilar la libreria se puede usar QT creator.

1. Primero se clona el proyecto en su computadora.
2. Abrir el archivo .pro con QT Creator.
3. Configurar la carpeta donde se compilará
4. Compilar el proyecto

### CMake

En caso de compilar con cmake asegurese de tener instalado la version `>3.0`

```
mkdir build
cd build
cmake ..
make
```

## Uso

Para usar esta libreria asegurese de haber compilado este proyecto y haber generado los archivos `.o` y `.a`.

La libreria se ha usado con los IDEs QT Creator y CodeBlocks.

> **Nota:** Algunos algoritmos hacen uso de procesamiento multi-hilo, para ellos asegurese de usar el flag `-pthread` e incluir la libreria `-lpthread`.

### QT
Para usarlo con QT inicialice un nuevo proyecto y asegurece que su archivo de configuracion `.pro` tenga lo siguiente:

```
INCLUDEPATH += /dip/classification \ # Aqui es el path del proyecto con las cabeceras.

LIBS += -L/dip/build \  # Directorio con los archivos compilados
    -lclassification \
    -ldcmdata \
    -ldcmimgle \
    -ldcmimage \
    -ldcmjpeg \
    -lpthread

```

### Codeblocks

Para usarlo en codeblocks inicie un nuevo proyecto y asegurece de ejecutar el script `checker.sh`.

Dirigirse a configuraciones y compilador, en la pestaña __linker settings__ agregue el archivo `libDicomClasifier.a` que se encuentra en la carpeta __build__ despues del proceso de compilacion. Tambien incluya las siguientes opciones:

```
-lclassification
-ldcmdata
-ldcmimgle
-ldcmimage
-ldcmjpeg
-lpthread
```

Ahora en la pestaña `Search directories > compiler` agregue el directorio donde se encuentran los archivos `.h` del proyecto.

```
/home/user/projects/DPI/DPI
```

En la pestaña `Search directories > linker` agregue el directorio donde se se encuentran los archivos compilados.

```
/home/usuario/dip/build
```

### CMake

Para IDEs que hace uso de CMAKE para la compilacion se puede ejecutar haciendo uso de la siguiente configuracion en el archivo `CMakeLists.txt`

```
cmake_minimum_required(VERSION 3.5)

project(knn_cmake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(classification STATIC IMPORTED)
set_target_properties(classification PROPERTIES
    IMPORTED_LOCATION "/home/usuario/dip/build/libclassification.a"
    INTERFACE_INCLUDE_DIRECTORIES "/home/usuario/dip/classification"
)

add_executable(knn_cmake main.cpp)
target_link_libraries(knn_cmake classification -ldcmdata -ldcmimgle -ldcmimage -ldcmjpeg -lpthread)


```


## Ejemplos

Puede encontrar ejemplos de uso [aquí](examples)

## Autores

- [Wilson Aguilar](https://www.github.com/waguilars)
- [Gabriel Cacuango]()
