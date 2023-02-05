#include <iostream>
#include "kmeans/kmeans.h"

#include "dicom/dicomutils.h"
#include "dicom/DicomReader.h"

using namespace std;

int main() {
    DicomReader dicomObj("/home/will/Projects/dicom-classifier/data/20586908_6c613a14b80a8591_MG_R_CC_ANON.dcm");
    vector<vector<double>> data = dicomObj.getDoubleImageMatrix(12);

    vector<Point> points = DicomUtils::getKMeansPoints(data);

    int k = 3;

    int iters = 5;

    KMeans kmeans(k, iters, ".");
    kmeans.run(points); // last

    string pref1 = "cluster";
    string pref2 = "points";
    kmeans.saveClusters(pref1);
    kmeans.savePoints(points, pref2);

    return 0;
}

