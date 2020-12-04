#ifndef IDWOPERATOR_H
#define IDWOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class ReclassifyOperator : public RasterOperator<double> 
{
public:
	ReclassifyOperator()
		:RasterOperator<double>(){};


    void setInputLayer(RasterLayer<double>& inputLayer);
    void setOutputLayer(RasterLayer<double>& outputLayer){_outputLayer=&outputLayer;}
private:
    bool Operator(const CellCoord &coord, bool operFlag) OVERRIDE;

	RasterLayer<double> *_inputLayer;
	RasterLayer<double> *_outputLayer;
    


};

#endif