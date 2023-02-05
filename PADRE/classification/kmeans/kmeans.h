#ifndef KMEANS_H
#define KMEANS_H


#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Point
{
private:
    int pointId, clusterId;
    int dimensions;
    vector<double> values;

    vector<double> lineToVec(string &line)
    {
        vector<double> values;
        string tmp = "";

        for (int i = 0; i < (int)line.length(); i++)
        {
            if ((48 <= int(line[i]) && int(line[i])  <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
            {
                tmp += line[i];
            }
            else if (tmp.length() > 0)
            {

                values.push_back(stod(tmp));
                tmp = "";
            }
        }
        if (tmp.length() > 0)
        {
            values.push_back(stod(tmp));
            tmp = "";
        }

        return values;
    }

public:
/**
 * @brief Construct a new Point object contruye un nuevo punto para calcular en cada iteracion
 * 
 * @param id identificador de nuevo punto.
 * @param line tamaño del vector para almacenar los puntos.
 */
    Point(int id, string line)
    {
        pointId = id;
        values = lineToVec(line);
        dimensions = values.size();
        clusterId = 0; // Initially not assigned to any cluster
    }

/**
 * @brief Construct a new Point object 
 * 
 * @param id identificador del nuevo punto
 * @param dataPoints datos de los puntos almacenados en un vector
 */

    Point(int id, vector<double> dataPoints) {
        pointId = id;
        values = dataPoints;
        dimensions = dataPoints.size();
        clusterId = 0;
    }
/**
 * @brief Get the Dimensions object devuelve las dimensiones de la matriz
 * 
 * @return int regresa el numero de dimensiones
 */

    int getDimensions() { return dimensions; }

/**
 * @brief Get the Cluster object regresa el numero de cluster 
 * 
 * @return int regresa el numero de clusters
 */

    int getCluster() { return clusterId; }
/**
 * @brief devuelve el ID del punto calculado
 * 
 * @return int valor del id del punto
 */

    int getID() { return pointId; }
/**
 * @brief Set the Cluster object envia el nuevo valor del cluster
 * 
 * @param val valor del clusterid a actualizar
 */

    void setCluster(int val) { clusterId = val; }

/**
 * @brief Get the Val object obtiene el valor del punto
 * 
 * @param pos la posicion del punto
 * @return double el valor del punto
 */

    double getVal(int pos) { return values[pos]; }
};

class Cluster
{
private:
    int clusterId;
    vector<double> centroid;
    vector<Point> points;

public:
/**
 * @brief Construct a new Cluster object
 * 
 * @param clusterId identificador del cluster
 * @param centroid identificador del centroide
 */
    Cluster(int clusterId, Point centroid)
    {
        this->clusterId = clusterId;
        for (int i = 0; i < centroid.getDimensions(); i++)
        {
            this->centroid.push_back(centroid.getVal(i));
        }
        this->addPoint(centroid);
    }

    /**
     * @brief Agrega un nuevo punto
     * 
     * @param p punto a agregar
     */

    void addPoint(Point p)
    {
        p.setCluster(this->clusterId);
        points.push_back(p);
    }
    /**
     * @brief Elimina un punto
     * 
     * @param pointId  identificador del punto
     * @return true condicion para eliminar o mantener el punto
     * @return false condicion para eliminar o mantener el punto
     */
    bool removePoint(int pointId)
    {
        int size = points.size();

        for (int i = 0; i < size; i++)
        {
            if (points[i].getID() == pointId)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Remuev todos los puntos
     * 
     */
    void removeAllPoints() { points.clear(); }

    /**
     * @brief Get the Id object
     * 
     * @return int valor del identificador del cluster
     */

    int getId() { return clusterId; }

    /**
     * @brief Get the Point object
     * 
     * @param pos valor de la posicion del punto
     * @return Point regresa el valor del punto en la posicion enviada
     */

    Point getPoint(int pos) { return points[pos]; }
    /**
     * @brief Get the Size object
     * 
     * @return int numero de puntos del vector
     */

    int getSize() { return points.size(); }
    /**
     * @brief Get the Centroid By Pos object
     * 
     * @param pos numero de la posicion del centroide
     * @return double valor del centroide en base a la posicion
     */

    double getCentroidByPos(int pos) { return centroid[pos]; }

    /**
     * @brief Set the Centroid By Pos object
     * 
     * @param pos valor de la posicion del vector a enviar
     * @param val valor del centroide a enviar
     */

    void setCentroidByPos(int pos, double val) { this->centroid[pos] = val; }
};

class KMeans
{
private:
    int K, iters, dimensions, total_points;
    vector<Cluster> clusters;
    string output_dir;

    void clearClusters()
    {
        for (int i = 0; i < K; i++)
        {
            clusters[i].removeAllPoints();
        }
    }

    int getNearestClusterId(Point point)
    {
        double sum = 0.0, min_dist;
        int NearestClusterId;
        if(dimensions==1) {
            min_dist = abs(clusters[0].getCentroidByPos(0) - point.getVal(0));
        }
        else
        {
            for (int i = 0; i < dimensions; i++)
            {
                sum += pow(clusters[0].getCentroidByPos(i) - point.getVal(i), 2.0);
                // sum += abs(clusters[0].getCentroidByPos(i) - point.getVal(i));
            }
            min_dist = sqrt(sum);
        }
        NearestClusterId = clusters[0].getId();

        for (int i = 1; i < K; i++)
        {
            double dist;
            sum = 0.0;

            if(dimensions==1) {
                dist = abs(clusters[i].getCentroidByPos(0) - point.getVal(0));
            }
            else {
                for (int j = 0; j < dimensions; j++)
                {
                    sum += pow(clusters[i].getCentroidByPos(j) - point.getVal(j), 2.0);
                    // sum += abs(clusters[i].getCentroidByPos(j) - point.getVal(j));
                }

                dist = sqrt(sum);
                // dist = sum;
            }
            if (dist < min_dist)
            {
                min_dist = dist;
                NearestClusterId = clusters[i].getId();
            }
        }

        return NearestClusterId;
    }

public:
/**
 * @brief Construct a new KMeans object
 * 
 * @param K valor para k
 * @param iterations numero de iteraicones a ejecutar
 * @param output_dir directorio de salida del resultado
 */
    KMeans(int K, int iterations, string output_dir)
    {
        this->K = K;
        this->iters = iterations;
        this->output_dir = output_dir;
    }

    void run(vector<Point> &all_points)

    {
        total_points = all_points.size();
        dimensions = all_points[0].getDimensions();

        // Initializing Clusters
        vector<int> used_pointIds;

        for (int i = 1; i <= K; i++)
        {
            while (true)
            {
                int index = rand() % total_points;

                if (find(used_pointIds.begin(), used_pointIds.end(), index) ==
                        used_pointIds.end())
                {
                    used_pointIds.push_back(index);
                    all_points[index].setCluster(i);
                    Cluster cluster(i, all_points[index]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }
        //        cout << "Clusters initialized = " << clusters.size() << endl << endl;

        cout << "Running K-Means Clustering.." << endl;

        int iter = 1;
        while (true)
        {
            cout << "Iter - " << iter << "/" << iters << endl;
            bool done = true;

            // Add all points to their nearest cluster
#pragma omp parallel for reduction(&&: done) num_threads(16)
            for (int i = 0; i < total_points; i++)
            {
                int currentClusterId = all_points[i].getCluster();
                int nearestClusterId = getNearestClusterId(all_points[i]);

                if (currentClusterId != nearestClusterId)
                {
                    all_points[i].setCluster(nearestClusterId);
                    done = false;
                }
            }

            // clear all existing clusters
            clearClusters();

            // reassign points to their new clusters
            for (int i = 0; i < total_points; i++)
            {
                // cluster index is ID-1
                clusters[all_points[i].getCluster() - 1].addPoint(all_points[i]);
            }

            // Recalculating the center of each cluster
            for (int i = 0; i < K; i++)
            {
                int ClusterSize = clusters[i].getSize();

                for (int j = 0; j < dimensions; j++)
                {
                    double sum = 0.0;
                    if (ClusterSize > 0)
                    {
#pragma omp parallel for reduction(+: sum) num_threads(16)
                        for (int p = 0; p < ClusterSize; p++)
                        {
                            sum += clusters[i].getPoint(p).getVal(j);
                        }
                        clusters[i].setCentroidByPos(j, sum / ClusterSize);
                    }
                }
            }

            if (done || iter >= iters)
            {
                cout << "Clustering completed in iteration : " << iter << endl
                     << endl;
                break;
            }
            iter++;
        }


    }

    vector<Cluster> getClustersValues() {
        return clusters;
    }
    /**
     * @brief almacena los nuevos valores de los puntos en el vector
     * 
     * @param all_points un vector que almacena todos los puntos 
     * @param prefix cadena de texto para almacenar junto con los puntos
     */

    void savePoints(vector<Point> &all_points ,string &prefix) {

        // Write Points
        ofstream pointsFile;
        pointsFile.open(output_dir + "/" + prefix + to_string(K) + "-points.txt", ios::out);

        for (int i = 0; i < total_points; i++)
        {
            pointsFile << all_points[i].getCluster() << endl;
        }

        pointsFile.close();

        return;
    }

    /**
     * @brief Almacena todos los valores de los clusters calculados 
     * 
     * @param prefix extension de los cluster almacenados
     */

    void saveClusters(string &prefix) {
        // Write cluster centers to file
        ofstream outfile;

        outfile.open(output_dir + "/" + prefix +to_string(K) + "-clusters.txt");
        if (outfile.is_open())
        {
            for (int i = 0; i < K; i++)
            {
//                cout << "Cluster " << clusters[i].getId() << " centroid : ";
                for (int j = 0; j < dimensions; j++)
                {
//                    cout << clusters[i].getCentroidByPos(j) << " ";    // Output to console
                    outfile << clusters[i].getCentroidByPos(j) << " "; // Output to file
                }
//                cout << endl;
                outfile << endl;
            }
            outfile.close();
        }
        else
        {
            cout << "Error: Unable to write to clusters.txt";
        }

        return;
    }
    /**
     * @brief Get the Avg Clusters object
     * 
     * @return double obtiene el valor de los centroides y calcula el porcentaje de average
     */

    double getAvgClusters() {

        double accum = 0;
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < dimensions; j++) {
                accum = accum + clusters[i].getCentroidByPos(j);
            }
        }
        return accum / (K*dimensions);
    }
};


#endif // KMEANS_H
