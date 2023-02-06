#ifndef IMAGENES_H
#define IMAGENES_H

#include <iostream>
#include <vector>

using namespace std;


class IMAGENES
{
public:
    IMAGENES();
    void inicio() const;
    void RUN(vector<vector<int>>  image_path, int min, int max) const;
};

#endif // IMAGENES_H
