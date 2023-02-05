#include "svm_util.h"

svm_util::svm_util()
{

}


svm_node *svm_util::createSVMNode(vector<int> points)
{

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

void svm_util::initProblemNodes(svm_node *x_space, svm_problem prob, vector<vector<int> > data)
{

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

