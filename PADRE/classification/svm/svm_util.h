#ifndef SVM_UTIL_H
#define SVM_UTIL_H
#include <vector>
#include <iostream>

#include "svm.h"


#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;

class svm_util
{
public:
    svm_util();

    static struct svm_node *createSVMNode(vector<int> points);

    static void initProblemNodes(svm_node *x_space, svm_problem prob, vector<vector<int>> data);

};

#endif // SVM_UTIL_H
