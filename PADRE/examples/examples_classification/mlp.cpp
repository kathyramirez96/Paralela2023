#include <iostream>
#include <fstream>

#include "tiny_dnn/tiny_dnn.h"
#include "tiny_dnn/config.h"
#include "dicom/DicomReader.h"
#include "opencv2/opencv.hpp"
#include "opencv2/core/hal/interface.h"

using namespace tiny_dnn;
using namespace std;

template<typename _Tp> static  cv::Mat toMat(const vector<vector<_Tp> > vecIn) {
    cv::Mat_<_Tp> matOut(vecIn.size(), vecIn.at(0).size());
    for (int i = 0; i < matOut.rows; ++i) {
        for (int j = 0; j < matOut.cols; ++j) {
            matOut(i, j) = vecIn.at(i).at(j);
        }
    }
    return matOut;
}

vector<vector<double> > toVec(const cv::Mat matIn) {
    vector<vector<double> > vecOut(matIn.rows);
    for (int i = 0; i < matIn.rows; ++i) {
        vecOut[i].resize(matIn.cols);
        for (int j = 0; j < matIn.cols; ++j) {
            vecOut[i][j] =  matIn.at<double>(i, j);
        }
    }
    return vecOut;
}

void saveWeightValues(string filename, std::vector<vec_t*> weights) {
      ofstream file;
      file.open(filename);

    for (int i = 0; i < weights[0]->size(); ++i) {
        file<< weights[0]->at(i) << endl;
    }
      file.close();
      return;
}


vec_t flatImage(vector<vector<double>> image) {

    vec_t flattedImage;
    for (auto row : image) {
        for (auto val : row) {
            flattedImage.push_back(val);
        }
    }

    return flattedImage;
}


int main()
{
    DicomReader dicomObj("/home/will/Projects/dicom-classifier/data/20586908_6c613a14b80a8591_MG_R_CC_ANON.dcm");

    vector<vector<double>> imageData = dicomObj.getDoubleImageMatrix(12);

    // Declare resize Mat
    cv::Mat img_resize;
    cv::Mat dicomMat = toMat(imageData); // Parse dicom vec to mat

    // Resize
    cv::resize(dicomMat, img_resize,cv::Size(512,512), 0, 0, cv::INTER_CUBIC);

    // Parse mat to vec
    vector<vector<double>> dicomResized = toVec(img_resize);
    vec_t flatImg = flatImage(dicomResized);

    vector<vec_t> input_data;
    vector<label_t> labels {0};

    input_data.push_back(flatImg);

    network<sequential> net;
    adam opt;

//    opt.alpha // Learning rate

    cout << "Adding layers"<<endl;
    net << layers::fc(512*512, 300) << activation::tanh()
        << layers::fc(300, 1) << activation::tanh();

    size_t batch_size = 1;
    size_t epochs = 1 ;

    cout << "Begin training"<<endl;
    net.train<mse>(opt, input_data, labels, batch_size, epochs);

//    net.save("model", content_type::model, file_format::json);

    for (int i = 0; i < net.depth(); i++) {
        cout << "#layer:" << i << "\n";
        cout << "layer type:" << net[i]->layer_type() << "\n";
        cout << "input:" << net[i]->in_data_size() << "(" << net[i]->in_data_shape() << ")\n";
        cout << "output:" << net[i]->out_data_size() << "(" << net[i]->out_data_shape() << ")\n";
    }


    // Print weight size and bias size for layer 1 (fc)
    //    cout << "fully connected layer 1" << endl;
//    std::vector<vec_t*> weights1 = net[0]->weights();
//    cout << "weight 1 size: " << weights1[0]->size() << endl;
//    cout << "bias size: " << weights1[1]->size() << endl;


    // Print weight size and bias size for layer 2 (fc)
//    std::vector<vec_t*> weights2 = net[2]->weights();
//    cout << "fully connected layer 2" << endl;
//    cout << "weight size: " << weights2[0]->size() << endl;
//    cout << "bias size: " << weights2[1]->size() << endl;

    // Print net in json
//    string json = net.to_json();
//    cout << json << endl;

    saveWeightValues("fc2-weights.txt", weights2);


    cout << "finished" << endl;

    return 0;
}
