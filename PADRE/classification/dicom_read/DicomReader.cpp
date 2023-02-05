#include "DicomReader.h"
#include <iostream>
#include <fstream>

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h" /* for DicomImage */
#include "dcmtk/dcmjpeg/djdecode.h"

using namespace std;
using namespace dcmtk;

void Test::run() {
    cout << "Hello" << endl;
}

DicomReader::DicomReader(const char *path)
{
    DJDecoderRegistration::registerCodecs();

    DicomImage *image = new DicomImage(path);
    this->image = image;
    this->width = this->image->getWidth();
    this->height = this->image->getHeight();
}

void DicomReader::saveCSV(const char *filename, int **data) {
    std::ofstream myfile;
    myfile.open (filename);
    for (int i = 0; i < this->getHeight(); ++i) {
        for (int j = 0; j < this->getWidth(); ++j) {
            myfile << data[i][j] << ";";
        }
        myfile << endl;
    }

    myfile.close();
    cout << filename << " has been created!" <<endl;
}

void DicomReader::clear()
{
    delete this->image;
}

std::vector<std::vector<float> > DicomReader::getFloatImageMatrix(int depth)
{
    int **imageArr = this->getImageArray(depth);

    vector<vector<float>> imageMatrix;

    for (int i = 0; i < this->height; ++i) {
        vector<float> row;
        for (int j = 0; j < this->width; ++j) {
            row.push_back(imageArr[i][j]);
        }
        imageMatrix.push_back(row);
    }

    for (int i = 0; i < this->height; ++i) {
        delete[] imageArr[i];
    }
    delete[] imageArr;

    return imageMatrix;

}

std::vector<std::vector<double> > DicomReader::getDoubleImageMatrix(int depth)
{

    int **imageArr = this->getImageArray(depth);

    vector<vector<double>> imageMatrix;

    for (int i = 0; i < this->height; ++i) {
        vector<double> row;
        for (int j = 0; j < this->width; ++j) {
            row.push_back(imageArr[i][j]);
        }
        imageMatrix.push_back(row);
    }

    for (int i = 0; i < this->height; ++i) {
        delete[] imageArr[i];
    }
    delete[] imageArr;

    return imageMatrix;

}

void DicomReader::saveData(std::vector<std::vector<int> > data, const char *filename, string delimiter, bool header)
{
    std::ofstream output_file(filename);

    // Write header
    if (header) {
        for (size_t i = 0; i < data[0].size(); ++i) {
            if (i == data[0].size() - 1) {
                output_file << "LABEL" << endl;
            } else {
                output_file << "X" << i << delimiter;
            }
        }

    }

    // Write data
    for (const auto &row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            if (i == row.size() - 1) {
                output_file << row[i];
            } else {
                output_file << row [i] << delimiter;
            }
        }
        output_file << "\n";
    }
}

int DicomReader::getWidth() { return this->width; }

void DicomReader::setWidth(int width) { this->width = width; }

int DicomReader::getHeight() { return this->height; }

void DicomReader::setHeight(int height) { this->height = height;}

int **DicomReader::getImageArray(int depth = 16) {
    int width = image->getWidth();
    int height = image->getHeight();

    int **imageArray;
    imageArray = new int*[height];

    for (int i = 0; i < height; ++i) {
        imageArray[i] = new int[width];
    }

    if(this->image == NULL) {
        cerr << "Error: no image has been provided" << endl;
        return NULL;
    }

    if (this->image->getStatus() == EIS_Normal) {
        if (image->isMonochrome())
        {

            image->setMinMaxWindow();
            const DiPixel* pixel = image->getInterData();
            for (int i = 0; i < height; ++i) {
                for (int j = 0; j < width; ++j) {
                    int pixel_value = *((short*)(pixel->getData()) + (i * image->getWidth() + j));
                    imageArray[i][j] = pixel_value;
                }
            }
        }
    } else {
        cerr << "Error: cannot load DICOM image (" << DicomImage::getString(image->getStatus()) << ")" << endl;
    }

    return imageArray;
}

std::vector<std::vector<int>> DicomReader::getIntImageMatrix(int depth)
{
    int **imageArr = this->getImageArray(depth);

    vector<vector<int>> imageMatrix;

    for (int i = 0; i < this->height; ++i) {
        vector<int> row;
        for (int j = 0; j < this->width; ++j) {
            row.push_back(imageArr[i][j]);
        }
        imageMatrix.push_back(row);
    }

    for (int i = 0; i < this->height; ++i) {
        delete[] imageArr[i];
    }
    delete[] imageArr;

    return imageMatrix;

}
