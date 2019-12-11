#ifndef IDWOPERATOR_H
#define IDWOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.0000001
#define  NODATA_DEFINE -9999

struct Sample_Point{
	int x;
	int y;
	double value;
};

struct extent_info {
	double minX;
	double maxX;
	double minY;
	double maxY;
};

struct Sample_block{
	//CoordBR _MBR;	//�Ƿ���Ҫ����
	//extent_info blockExtent;
	vector<Sample_Point> sample_Points;
};

class IDWOperator : public RasterOperator<double> 
{
public:
	IDWOperator()
		:RasterOperator<double>(),
		 _iterNum(0), flag(true), _sample_nums(0)
		{}

	//�����ĳ�ʼ��
	IDWOperator( float cellsize, int nbrPoints, int idw_power, int bufferSize, int grain )
		:RasterOperator<double>(),
		_iterNum(0), flag(true), _sample_nums(0), _cellSize(cellsize), _nbrPoints(nbrPoints), _idw_power(idw_power), _idw_buffer(bufferSize), _blockGrain(grain), _noData(NODATA_DEFINE)
	{}

	~IDWOperator();

    int getBlockGrain(){return _blockGrain;}

	int readSampleNums( const char* filename, char** pSpatialRefWkt );
	bool readSamples( const char* filename, int fieldIdx, char** pSpatialRefWkt, double **Sample_Array );
	void creatSampleBlocks( double **pSamples );
    Sample_block* getSampleBlocks(){return _pSampleBlocks;}

	void idwLayer(RasterLayer<double> &layerD, char** pSpatialRefWkt);
	virtual bool isTermination();

	int searchNbrSamples( const int subMinRow, int cellRow, int cellCol, double *nbrSamples);

	virtual bool Operator(const CellCoord &coord, bool operFlag);

    int getBlockRowIndex(double x);
    int getBlockColIndex(double y);

    int getBlockRowIndex(int iRow);
    int getBlockColIndex(int iCol);

    int getBlockCols(){return _blockCols;}
    int getBlockRows(){return _blockRows;}

    int getNbrPoints(){return _nbrPoints;}
private:
	float _cellSize;
	int _xSize, _ySize;	//��ǰ������DEM�������
	int _nRows, _nCols; //����ͼ����������������
	double _noData;
	int _myRank;
	int _iterNum;//��������
	bool flag;
protected:
	int _nbrPoints;	//��ֵ����������
	int _idw_power;
	int _idw_buffer;
	extent_info _glb_extent;	//�������Ĳ�ֵդ��ͼ��ȫ����Χ�������������ݷ�Χ����ȷ��
	extent_info _sub_extent;	//������դ��Χ
	int _sample_nums;	//ȫ������������
	int _blockGrain;
	Sample_block* _pSampleBlocks;	//�Դ������,������֯������;ÿ�����̶�����ȫ���飬����������Ҫ�ͷţ�
    int _blockRows;
    int _blockCols;
	RasterLayer<double> *_pIDWLayer;


};

#endif