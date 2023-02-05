#ifndef RF_H
#define RF_H

#include <vector>
#include <string>

#include "random_forest/globals.h"
#include "random_forest/ForestClassification.h"
#include "random_forest/utility.h"
#include "random_forest/Data.h"

using namespace std;
using namespace ranger;

class RF
{
public:
    std::unique_ptr<Forest> forest;
    int totalTrees;
/**
 * @brief Construct a new RF object
 * 
 * @param trees numero de arboles a crear
 */
    RF(int trees);
    ~RF();
/**
 * @brief Set the Train Data object
 * 
 * @param data matriz inicial para entrenamiento
 * @param target etiquetas del entrenamiento
 */
    void setTrainData(vector<vector<double>> data, vector<int> target);
/**
 * @brief metodo para ver las predicciones de la tecnica
 * 
 * @param data vector con los valores de inciales
 * @param target etiquetas de entrenamiento
 * @param showOutput resultados obtenidos despues del proceso
 */
    void predict(vector<vector<double>> data, vector<int> target, bool showOutput);
    void init(bool showOutput);
/**
 * @brief Set the Train File object
 * 
 * @param filename nombre del archivo 
 */
    void setFile(string filename);
    /**
     * @brief Set the Predict File object
     * 
     * @param filename nombre del archivo
     */
    void setPredictFile(string filename);
/**
 * @brief Set the Dep Var Name object
 * 
 * @param name nombre del archivo
 */
    void setDepVarName(string name);

private:
        // Config
        uint mtry;
        MemoryMode mode;
        string  filename;
        string outprefix;
        string depvarname;
//        std::ofstream logfile { arg_handler.outprefix + ".log" };
        uint default_seed;
        string predict_file;
        string split_weights_file;
        vector<std::string> split_vars;
        string status_var_name;
        bool replacement;
        vector<std::string> cat_vars;
        bool save_memory;
        string weights_file;
        bool predall;
        double samplefraction;
        bool holdout;
        vector<double> reg_factor;
        bool reg_usedepth;

    string generateDataFormat(vector<vector<double>> data, vector<int> target);
};

#endif // RF_H
