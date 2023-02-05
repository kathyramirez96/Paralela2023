#ifndef DICOMREADER_H
#define DICOMREADER_H

#include<vector>
#include <string>

class DicomImage;

class DicomReader
{
private:
    DicomImage *image;
    int width;
    int height;

public:
/**
 * @brief Construct a new Dicom Reader object
 * 
 * @param path Dirreccion de la imagen DICOM 
 */
    DicomReader(char const *path);
    /**
     * @brief Get the Width object
     * 
     * @return int ancho de la imagen 
     */
    int getWidth();
    /**
     * @brief Set the Width object
     * 
     * @param width envia el ancho de la imagen
     */
    void setWidth(int width);

    /**
     * @brief Get the Height object
     * 
     * @return int alto de la imagen
     */
    int getHeight();
    void setHeight(int height);
    int **getImageArray(int depth);
    /**
     * @brief Guardar los resultados en un csv
     * 
     * @param filename nombre del archivo para almacenar el csv
     * @param data informacion que se escribe en el csv
     */
    void saveCSV(const char *filename, int **data);
    /**
     * @brief  elimina la imagen actual para limpiar memoria
     * 
     */
    void clear();

    std::vector<std::vector<int>> getIntImageMatrix(int depth);

    std::vector<std::vector<float>> getFloatImageMatrix(int depth);

    std::vector<std::vector<double>> getDoubleImageMatrix(int depth);

    static void saveData(std::vector<std::vector<int>> data, const char* filename, std::string delimiter, bool header);


};



class Test
{
public:
    void run();
};

#endif // DICOMREADER_H
