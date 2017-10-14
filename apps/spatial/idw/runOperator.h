#ifndef RUNOPERATOR_H
#define RUNOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.0000001
class RUNOperator : public RasterOperator<double> 
{
public:
	RUNOperator()
		:RasterOperator<double>(),
		_pDEMLayer(0), _pPitLayer(0), num(0), flag(true) {}
	//�����ĳ�ʼ��

	~RUNOperator() {}

	void demLayer(RasterLayer<double> &layerD);
	void pitLayer(RasterLayer<double> &layerD);
	
	virtual bool isTermination();

	virtual bool Operator(const CellCoord &coord, bool operFlag);


private:
	int _cellSize;
	int _xSize, _ySize;	//��ǰ������DEM�������
	int _nRows, _nCols; //����ͼ����������������
	double _noData;
	int _rank;
	int num;//��������
	bool flag;

	RasterLayer<double> * _pDEMLayer;
	RasterLayer<double> * _pPitLayer;

	Neighborhood<double> *_pDEMNbrhood;

};

#endif