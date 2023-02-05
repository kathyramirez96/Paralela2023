#ifndef KNN_H
#define KNN_H

#include "dataset.h"

class SingleExecutionResults {
public:
	/**
	 * @brief Construct a new Single Execution Results object
	 * 
	 * @param nExamples numero de ejemplos 
	 * @param nSuccess numero de sucesos
	 * @param nRejected numero de rechazados
	 */
	SingleExecutionResults(int nExamples, int nSuccess, int nRejected) {
		this->nExamples = nExamples;
		this->nSuccess = nSuccess;
		this->nRejected = nRejected;
	}

	double successRate() { return ((double)nSuccess) / (nExamples - nRejected);}
	double rejectionRate() { return ((double)nRejected) / nExamples; }

	int nExamples;
	int nSuccess;
	int nRejected;
};


class KNNResults {
public:
	SingleExecutionResults top1Result();
	SingleExecutionResults topXResult(int n);
	MatrixPointer getConfusionMatrix();
	MatrixPointer getPredictions();

	
	DatasetPointer getRawResults() {return results;}

	/**
	 * @brief Construct a new KNNResults object
	 * 
	 * @param results almacena los resultados con los puntos trabajados
	 */

	KNNResults(DatasetPointer results) { this->results = results;}

private:
	DatasetPointer results;

};


class KNN {
public:
/**
 * @brief Construct a new KNN object
 * 
 * @param train conjunto de datos de los puntos con los cuales se ejecuta la tecnica
 */
	KNN (DatasetPointer train) { this->data = train; }

	KNNResults run(int k, DatasetPointer target);

private:
	DatasetPointer data;
};

#endif /* KNN_H */
