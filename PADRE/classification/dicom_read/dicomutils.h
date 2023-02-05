#ifndef DICOMUTILS_H
#define DICOMUTILS_H

#include <vector>
#include <string>
#include <iostream>

#include "kmeans/kmeans.h"

using namespace std;

class DicomUtils
{
public:
    DicomUtils();
    // Return list with path of dicom files from and directory
    static vector<string> getDicomFilesPath(const char *path);

    // return a list of dicom files path based in file list
    // Format:
    // filename.dcm, label
    static vector<string> makeDicomFilesPath(string dicomDirectory, string summaryFile);

    // Return a list with the labels of each dcmfile
    static vector<int> makeDicomLabels(string filepath);

    // Read labels generated based on the roi ([filename].target)
    // path -> directory of target values
    // DicomFile path -> Path of dicom file to get target values
    static vector<int> getDicomTargetRoi(const char *path, string dicomFilePath);

    static int getDataWidth(vector<string> dicomFiles);
    static int getDataHeight(vector<string> dicomFiles);
    static string base_name(string const & path, string const &delims);


    static double **asFCMPointsData(int **arr, int numPoints, int numDims);
    static void writeMetrics(ofstream &file, int i, string filename, double cpu, int mem, double time, double param);

    static vector<Point> getKMeansFilePoints(string filename);
    static vector<Point> getKMeansPoints(vector<vector<double>> data);

    static vector<int> readLabels(string filename, string delimiter);

    static vector<int> genTargetValues(vector<vector<double>> data, int numClass);
    static vector<int> genTargetValuesInt(vector<vector<int>> data, int numClass);


    static vector<int> genTestDataIdx(vector<vector<double>> data, int size);
    static vector<int> genTestDataIdxInt(vector<vector<int>> data, int size);

    // Generate values and labels for testing using index
    static vector<int> getTestingLabels(vector<int>labels, vector<int> indexes);
    static vector<vector<double>> getTestingValues(vector<vector<double>> data, vector<int> indexes);
    static vector<vector<int>> getTestingValuesInt(vector<vector<int>> data, vector<int> indexes);

    // Parse to KNN input values
    static double *parseKNNData(vector<vector<double>>data);
    static int *parseKNNLabels(vector<int> labels);

    // Save vector label
    static void saveFile(vector<int> dataList, string delimiter, string filename);
    // Save matrix data
    static void saveData(std::vector<std::vector<double>> data, string delimiter, string filename, bool header);
    static void saveDataInt(std::vector<std::vector<int>> data, string delimiter, string filename, bool header);

};

#endif // DICOMUTILS_H
