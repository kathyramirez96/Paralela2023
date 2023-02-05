#ifndef FCM_H
#define FCM_H

#define MAX_DATA_POINTS 10000
#define MAX_CLUSTER 100
#define MAX_DATA_DIMENSION 5000

class FCM
{
public:
    FCM(double fuzziness, double epsilon);
    /**
     * @brief Inicializa con la matriz de la imagen, el numero de clusters, el numero de puntos, numero de dimensiones
     * 
     * @param data matriz de la imagen 
     * @param clusters numero de clusters a calcular
     * @param num_points numero de puntos deseados
     * @param num_dimensions dimensiones de alcance de la tecnica
     */
    void init(double **data, int clusters, int num_points, int num_dimensions);
    /**
     * @brief Evalua el valor para verificar cambios en la tecnica
     * 
     */
    void eval();
    /**
     * @brief Get the Centers object obtiene los centroides en base a los clusters
     * 
     * @return double** 
     */
    double **getCenters();
    double **getMembershipMatrix();
    /**
     * @brief Muestra los puntos del valor y los muestra en consola
     * 
     */
    void printPoints();
    /**
     * @brief almacena los nuevos clusters
     * 
     * @param prefix  
     */
    void saveClusters(const char *prefix);
    void saveMembershipMatrixU(const char *name);
    /**
     * @brief Almacena los centroides calculados
     * 
     * @param name almacena con los nombres enviados
     */
    void saveCenters(const char *name);
    /**
     * @brief Get the Center A V G object calcula el averague de los valores del centro 
     * 
     * @return double 
     */
    double getCenterAVG();


private:
    int num_data_points;
    int num_clusters;
    int num_dimensions;
    double low_high[MAX_DATA_DIMENSION][2];
    double degree_of_memb[MAX_DATA_POINTS][MAX_CLUSTER];
    double epsilon;
    double fuzziness;
    double data_point[MAX_DATA_POINTS][MAX_DATA_DIMENSION];
    double cluster_centre[MAX_CLUSTER][MAX_DATA_DIMENSION];


    void calculate_centre_vectors();
    double update_degree_of_membership();
    double get_new_value(int i, int j);
    double get_norm(int i, int j);
};

#endif // FCM_H
