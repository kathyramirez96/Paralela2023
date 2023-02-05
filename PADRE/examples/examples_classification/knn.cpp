#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "dicom/DicomReader.h"
#include "dicom/dicomutils.h"

#include "knn/dataset.h"
#include "knn/knn.h"
#include "knn/knnUtils.h"
#include "knn/Preprocessing.h"

using namespace std;

int main()
{
    DicomReader dicomObj("/home/will/Projects/dicom-classifier/data/20586908_6c613a14b80a8591_MG_R_CC_ANON.dcm");

    vector<vector<double>> data= dicomObj.getDoubleImageMatrix(12);
    vector<int> dataLabels = DicomUtils::genTargetValues(data, 2);

    vector<int> testIdx = DicomUtils::genTestDataIdx(data, 50);
    vector<int> testLabels = DicomUtils::getTestingLabels(dataLabels, testIdx);
    vector<vector<double>> testData = DicomUtils::getTestingValues(data, testIdx);

    double *parsedTrainData = DicomUtils::parseKNNData(data);
    int *parsedTrainLabel = DicomUtils::parseKNNLabels(dataLabels);

    double *parsedTestData = DicomUtils::parseKNNData(testData);
    int *parsedTestLabel = DicomUtils::parseKNNLabels(testLabels);


    int rows = data.size();
    int cols = data[0].size();

    int numLabels = 2;
    DatasetPointer dataset = makeDataset(rows, cols, numLabels, parsedTrainData, parsedTrainLabel);
    DatasetPointer test = makeDataset(testData.size(), cols, numLabels, parsedTestData, parsedTestLabel);

    MatrixPointer meanData = MeanNormalize(dataset);

    KNN knn(dataset);

    ApplyMeanNormalization(test, meanData);


    int bestK = 2;
    KNNResults results = knn.run(bestK, test);
    cout << "Consolidating results" << endl;

    SingleExecutionResults top1 = results.top1Result();
    SingleExecutionResults top2 = results.topXResult(2);
    SingleExecutionResults top3 = results.topXResult(3);

    double successRate = top1.successRate();


//    MatrixPointer predictions = results.getPredictions();


   // Confusion matrix
//    MatrixPointer confusionMatrix = results.getConfusionMatrix();
//    for(size_t i = 0; i< confusionMatrix->rows; i++) {
//        for(size_t j = 0; j< confusionMatrix->cols; j++) {
//            if (j!=0) printf(",");
//            printf("%d", (int)confusionMatrix->pos(i,j));
//        }
//        printf("\n");
//    }


    // Read Predicted values for test
//    for (size_t i = 0; i < predictions->rows; ++i) {
//        for (size_t j = 0; j < predictions->cols; ++j) {
//            cout<< "Pred" << predictions->pos(i, j) << " ";
//            cout<< "Pred" << predictions->pos(i, j) << " ";
//        }

//        cout << endl;
//    }

}
