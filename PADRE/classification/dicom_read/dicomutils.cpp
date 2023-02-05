#include "dicomutils.h"
#include "DicomReader.h"
#include <dirent.h>
#include <fstream>
#include <sstream>


#include "kmeans/kmeans.h"

using namespace std;

vector<string> DicomUtils::getDicomFilesPath(const char *path) {
    vector<string> files;
    struct dirent *d;
    DIR *dr;
    dr = opendir(path);
    if(dr!=NULL)
    {

        for(d=readdir(dr); d!=NULL; d=readdir(dr))
        {
            string filename = d->d_name;
            if(filename.substr(filename.find_last_of(".") + 1) == "dcm") {
                string fullpath = path;
                fullpath += "/" + filename;
                files.push_back(fullpath);
            }
        }
        closedir(dr);
    }
    else {
        cout<<"\nError Occurred when reading dir!" << endl;
        throw std::runtime_error("Error Occurred when reading dir");
        return files;
    }
    return files;
}

vector<string> DicomUtils::makeDicomFilesPath(string dicomDirectory, string summaryFile)
{
    vector<string> dcmFiles;
    string line;

    fstream file (summaryFile, ios::in);
    if(file.is_open())
    {
        while(getline(file, line))
        {
            vector<string> content;
            stringstream str(line);
            string word;

            while(getline(str, word, ','))
                content.push_back(word);

            string fullpath = dicomDirectory + "/" + content[0];
            dcmFiles.push_back(fullpath);

        }
    }
    else{
        cout<<"Could not open the file\n";
        throw std::runtime_error("Error Occurred when reading dir");
    }


    file.close();

    return dcmFiles;
}

vector<int> DicomUtils::makeDicomLabels(string filepath)
{
    vector<int> dcmLabels;
    string line;

    fstream file (filepath, ios::in);
    if(file.is_open())
    {
        while(getline(file, line))
        {
            vector<string> content;
            stringstream str(line);
            string word;

            while(getline(str, word, ','))
                content.push_back(word);

            dcmLabels.push_back(stoi(content[1]));

        }
    }
    else
        cout<<"Could not open the file\n";

    file.close();

    return dcmLabels;
}

vector<int> DicomUtils::getDicomTargetRoi(const char *path, string dicomFilePath)
{
    vector<int> target;
    struct dirent *d;
    DIR *dr;
    dr = opendir(path);
    if(dr!=NULL)
    {

        for(d=readdir(dr); d!=NULL; d=readdir(dr))
        {
            string filename = d->d_name;

            if(filename.substr(filename.find_last_of(".") + 1) == "target") {
                // remove ext from .target
                size_t lastindex = filename.find_last_of(".");
                string targetraw = filename.substr(0, lastindex);


                // remove ext from dicom file
                size_t lastindex2 = dicomFilePath.find_last_of(".");
                string dicomraw = dicomFilePath.substr(0, lastindex2);

                // if exist target file load it else load test values
                if(dicomraw == targetraw) {
                    string targetpath = string(path)+"/"+filename.c_str();
                    target = DicomUtils::readLabels(targetpath, ",");

                }
            }
        }
        closedir(dr);
    }
    else {
        cout<<"\nError Occurred when reading dir!" << endl;
        return target;
    }

    return target;



}

// return width for the data set
int DicomUtils::getDataWidth(vector<string> dicomFiles)
{
    int maxW = 0;
    for (size_t i = 0; i < dicomFiles.size(); ++i) {
        DicomReader dcm(dicomFiles[i].c_str());

        if (dcm.getWidth() > maxW) {
            maxW = dcm.getWidth();
        }
    }

    return maxW;
}

// return height for the data set
int DicomUtils::getDataHeight(vector<string> dicomFiles)
{
    int maxH = 0;
    for (size_t i = 0; i < dicomFiles.size(); ++i) {
        DicomReader dcm(dicomFiles[i].c_str());


        if (dcm.getHeight() > maxH){
            maxH = dcm.getHeight();
        }
    }

    return maxH;

}

// Return basename from path file
string DicomUtils::base_name(const string &path, const string &delims)
{      return path.substr(path.find_last_of(delims) + 1);
}

// return dcm data as array pointer (for FCM)
double **DicomUtils::asFCMPointsData(int **arr, int numPoints, int numDims)
{
    double **data;
    data= new double*[numPoints];

    for (int i = 0; i < numPoints; ++i) {
        data[i] = new double[numDims];
        for (int j = 0; j < numDims; ++j) {
            data[i][j] = arr[i][j];
        }
    }
    return data;
}


// Write metrics to file instance
// i, nombre, CPU, memoria, tiempo, alg_param
// alg_param: media del parametro a evaluar
void DicomUtils::writeMetrics(ofstream &file, int i, string filename, double cpu, int mem, double time, double param)
{
    file << i << "," << filename << "," << cpu << "," << mem << "," << time << "," << param << endl;
}

vector<Point> DicomUtils::getKMeansPoints(vector<vector<double> > data)
{
    vector<Point> all_points;
    for (int i = 0; i < data.size(); ++i) {
        Point point(i+1, data[i]);
        all_points.push_back(point);
    }
    return all_points;
}

vector<int> DicomUtils::readLabels(string filename, string delimiter)
{
    ifstream file(filename);
    string line = "";
    vector<int> labels;

    // Read label values and push to vector
    while (getline(file, line))
    {
        size_t pos = 0;
        string token;
        while ((pos = line.find(delimiter)) != string::npos) {
            token = line.substr(0, pos);
            line.erase(0, pos + delimiter.length());
            labels.push_back(stoi(token) );
        }
        labels.push_back(stoi(line));
    }
    // Close the File
    file.close();

    return labels;
}

vector<int> DicomUtils::genTargetValues(vector<vector<double> > data, int numClass)
{
    vector<int> values;
    for (size_t i = 0; i < data.size(); ++i) {
        double n = 0;
        for (size_t j = 0; j < data[0].size(); ++j) {
            n += data[i][j];
        }

        int m =(int)(n/data.size());
        values.push_back(m%numClass);
    }

    return values;
}

vector<int> DicomUtils::genTargetValuesInt(vector<vector<int> > data, int numClass)
{
    vector<int> values;
    for (size_t i = 0; i < data.size(); ++i) {
        double n = 0;
        for (size_t j = 0; j < data[0].size(); ++j) {
            n += data[i][j];
        }

        int m =(int)(n/data.size());
        values.push_back(m%numClass);
    }

    return values;
}

vector<int> DicomUtils::genTestDataIdx(vector<vector<double>> data, int size)
{
    vector<int> idx;
    int max = data.size()-1;
    int min = 0;
    int range = max - min + 1;

    for (int i = 0; i < size; ++i) {
        int num = rand() % range + min;
        idx.push_back(num);
    }

    return idx;

}

vector<int> DicomUtils::genTestDataIdxInt(vector<vector<int> > data, int size)
{
    vector<int> idx;
    int max = data.size()-1;
    int min = 0;
    int range = max - min + 1;

    for (int i = 0; i < size; ++i) {
        int num = rand() % range + min;
        idx.push_back(num);
    }

    return idx;
}

vector<int> DicomUtils::getTestingLabels(vector<int> labels, vector<int> indexes)
{

    vector<int> testingLabels;

    for (size_t i = 0; i < indexes.size(); ++i) {
        testingLabels.push_back(labels[indexes[i]]);
    }

    return testingLabels;

}

vector<vector<double> > DicomUtils::getTestingValues(vector<vector<double>> data, vector<int> indexes)
{
    vector<vector<double>> testingValues;
    for (size_t i = 0; i < indexes.size(); ++i) {
        testingValues.push_back(data[indexes[i]]);
    }

    return testingValues;

}

vector<vector<int> > DicomUtils::getTestingValuesInt(vector<vector<int>> data, vector<int> indexes)
{
    vector<vector<int>> testingValues;
    for (size_t i = 0; i < indexes.size(); ++i) {
        testingValues.push_back(data[indexes[i]]);
    }

    return testingValues;
}

double *DicomUtils::parseKNNData(vector<vector<double> > data)
{
    int rows = data.size();
    int cols = data[0].size();
    double *trainData = new double[rows*cols];

    int k = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            trainData[k] = data[i][j];
            k++;
        }
    }

    return trainData;
}

int *DicomUtils::parseKNNLabels(vector<int> labels)
{
    int rows = labels.size();
    int *trainLabels = new int[rows];
    for (int i = 0; i < labels.size(); ++i) {
        trainLabels[i] = labels[i];
    }

    return trainLabels;
}

void DicomUtils::saveFile(vector<int> dataList, string delimiter, string filename)
{
    ofstream myfile;
    myfile.open (filename);

    for (int i = 0; i < dataList.size(); ++i) {
        if (i == dataList.size()-1) {
            myfile << dataList[i];
        } else {
            myfile << dataList[i] << delimiter;
        }

    }
    myfile.close();

}

void DicomUtils::saveData(std::vector<std::vector<double> > data, string delimiter, string filename, bool header)
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
    output_file.close();
}

void DicomUtils::saveDataInt(std::vector<std::vector<int> > data, string delimiter, string filename, bool header)
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
    output_file.close();

}


vector<Point> DicomUtils::getKMeansFilePoints(string filename)
{
    ifstream infile(filename.c_str());
    if (!infile.is_open())
    {
        cout << "Error: Failed to open file: " << filename << endl;
        return {};
    }
    int pointId = 1;
    vector<Point> all_points;
    string line;

    while (getline(infile, line))
    {
        Point point(pointId, line);
        all_points.push_back(point);
        pointId++;
    }

    infile.close();
    cout << "\nData fetched successfully!" << endl;

    return all_points;
}
