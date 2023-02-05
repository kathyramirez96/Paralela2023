#include <vector>
#include <iostream>
#include <fstream>


#include "dicom/DicomReader.h"
#include "svm/svm.h"



#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;

using namespace std;

vector<vector<int> > parseData(int **img, int size, int feat)
{
    std::vector<std::vector<int> > data;

    for (int i=0; i< size; ++i)
    {
        std::vector<int> featureSet;
        for (int j=0; j< feat; ++j)
        {
            featureSet.push_back(img[i][j]);
        }

        data.push_back(featureSet);
    }
    return data;
}

vector<int> generateLabels(int labelsSize)
{
    std::vector<int> labels;
    for (int i=0; i<labelsSize; ++i)
    {
        int lbl = i%4;
        labels.push_back(lbl);
    }
    return labels;
}



void saveData (vector<vector<int>> data, const char* filename) {
    std::ofstream output_file(filename);

    for (const auto &row : data) {
        for (const auto &val: row) {
         output_file << val << ",";
        }
        output_file << "\n";
    }
}

void saveLabels (vector<int> labels) {
    std::ofstream output_file("labels.csv");

    for (const auto &v: labels) {
        output_file << v << ",";
    }
}

void initProblemNodes(svm_node *x_space, svm_problem prob, vector<vector<int>> data) {
    //initialize the svm_node vector with input data array as follows:
    int j=0; //counter to traverse x_space[i];
    for (int i=0;i < prob.l; ++i)
    {
        //set i-th element of prob.x to the address of x_space[j].
        //elements from x_space[j] to x_space[j+data[i].size] get filled right after next line
        prob.x[i] = &x_space[j];
        for (int k=0; k<data[i].size(); ++k, ++j)
        {
            x_space[j].index=k+1; //index of value
            x_space[j].value=data[i][k]; //value
        }
        x_space[j].index=-1;//state the end of data vector
        x_space[j].value=0;
        j++;
    }
}

struct svm_node *createSVMNode(vector<int> points) {
    struct svm_node *node;

    node = Malloc(struct svm_node, (points.size()+1) * 1);
    int j = 0;

    for (int i = 0; i < points.size(); ++i) {
        node[j].index =  i+1;
        node[j].value = points[i];
//        cout << node[j].index  << ": " << node[j].value << endl;
        j++;
    }
    node[j].index =  -1;
    node[j].value = 0;

    return node;
}



int main(int argc, char *argv[])
{
    DicomReader dicomObj("/home/will/Projects/dicom-classifier/data/20586908_6c613a14b80a8591_MG_R_CC_ANON.dcm");
    int size = dicomObj.getHeight();
    int elements = dicomObj.getWidth();

    vector<vector<int>> data = parseData(dicomObj.getImageArray(12), size, elements);
    vector<int> labels = generateLabels(size);

    vector<int> testIdx = {432, 2240, 860, 678, 1567};
    vector<int> train_labels;

    vector<vector<int>> test;
    for (auto id : testIdx) {
        test.push_back(data[id]);
        train_labels.push_back(labels[id]);
    }


    cout<<"data size = "<<data.size()<<endl;
    cout<<"labels size = "<<labels.size()<<endl;



//     init problem
    prob.l = size;

    // @param prob.l = number of labels
    // @param elements = number of features for each label
    prob.y = Malloc(double,prob.l); //space for prob.l doubles
    prob.x = Malloc(struct svm_node *, prob.l); //space for prob.l pointers to struct svm_node
    x_space = Malloc(struct svm_node, (elements+1) * prob.l); //memory for pairs of index/value//here we need to give some memory to our structures

    //initialize the different lables with an array of labels
    for (int i=0; i < prob.l; ++i)
        prob.y[i] = labels[i];

    // Initialize problem nodes with data values
    initProblemNodes(x_space, prob, data);

    //set all default parameters for param struct
    param.svm_type = NU_SVC;
    param.kernel_type = LINEAR;

    param.nu = 0.5;
    param.cache_size = 1000;
    param.C = 1;
    param.eps = 1e-3;
    param.p = 0.1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;

    //try to actually execute it
    model = svm_train(&prob, &param);

//     svm_set_print_string_function(NULL);
     cout << "Total classes: " << model->nr_class << endl;
     cout << model->label[0] << endl;
     cout << model->label[1] << endl;

     // Predict Labels
     vector<int> test_result;

     for (int i = 0; i < test.size(); ++i) {
         struct svm_node *point = createSVMNode(test[i]);
         int lbl = svm_predict(model, point);
         test_result.push_back(lbl);

         cout <<  "train: " << train_labels[i] << "| predict: " << lbl << endl;
     }


    return 0;
}
