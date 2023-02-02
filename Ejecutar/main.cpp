#include <QCoreApplication>
#include <iostream>
//#include <libimg.h>
#include <imagenes.h>

using namespace std;
int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    cout << "------COMPUTACIÓN PARALELA-------\n";
    IMAGENES img;
    img.inicio();
    //LIBIMG img2;
    //img2.ExecHilo1();
    return 0;
}
