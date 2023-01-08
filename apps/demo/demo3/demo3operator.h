#ifndef DEMO3_H
#define DEMO3_H

#include "utility.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;
#define Eps 0.0000001

class Operator : public RasterOperator<double> 
{
  public:
    Operator()
      :RasterOperator<double>(),
       cellSize(0), noData(-9999.), num(0),
       _pDEMLayer(nullptr), _pssLayer(nullptr),_psmoLayer(nullptr){}
   
    ~Operator() DEFAULT;
  
    void demLayer(RasterLayer<double> &layerD);
	void ssLayer(RasterLayer<double> &layerD);
	void smoLayer(RasterLayer<double> &layerD);

	bool isTermination() OVERRIDE;
    bool Operator(const CellCoord &coord, bool operFlag) OVERRIDE;

  protected:
	double cellSize;
	double noData;
	int num;
	//float p;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pssLayer;
	RasterLayer<double> *_psmoLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif