#ifndef IDWOPERATOR_H
#define IDWOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define EPS 0.0000001
#define  NODATA_DEFINE -9999

struct SamplePoint{
	double x;
	double y;
	double value;
    SamplePoint() {
        // cout<<"I'm inited!"<<endl;
    }
};

struct extent_info {
	double minX;
	double maxX;
	double minY;
	double maxY;
};

struct SampleBlock{
	//CoordBR _MBR;	//�Ƿ���Ҫ����
	//extent_info blockExtent;
	vector<SamplePoint> samplePoints;
};

class IDWOperator : public RasterOperator<double> 
{
public:
	IDWOperator()
		:RasterOperator<double>(),
		 _iterNum(0), flag(true), _sample_nums(0)
		{}

	//�����ĳ�ʼ��
	IDWOperator( float cellsize, int nbrPoints, int idw_power, double bufferSize, double blockSize )
		:RasterOperator<double>(),
		_iterNum(0), flag(true), _sample_nums(0), _cellSize(cellsize), _nbrPoints(nbrPoints), _idw_power(idw_power), _idw_buffer(bufferSize), _blockSize(blockSize), _noData(NODATA_DEFINE)
	{}

	~IDWOperator();

    double getBlockSize(){return _blockSize;}

	int readSampleNums( const char* filename, char** pSpatialRefWkt );
	//bool readSamples( const char* filename, int fieldIdx, char** pSpatialRefWkt, double **Sample_Array );
	//void creatSampleBlocks( double **pSamples );
	bool readSamples( const char* filename, int fieldIdx, char** pSpatialRefWkt, vector<SamplePoint> &samples );
	void creatSampleBlocks(vector<SamplePoint> &samples);
    const vector<SampleBlock>* getSampleBlocks(){return &_pSampleBlocks;}

	void idwLayer(RasterLayer<double> &layerD, char** pSpatialRefWkt, DomDcmpType dcmpType=ROWWISE_DCMP);
	void idwLayer(RasterLayer<double> &layerD, char** pSpatialRefWkt,CoordBR& subWorkBR);
	void maskLayer(RasterLayer<int> &layerD);
	void maskLayer(RasterLayer<int> &layerD, CoordBR& subWorkBR);
	void idwLayerSerial(RasterLayer<double> &layerD, char** pSpatialRefWkt);
	virtual bool isTermination();

	int searchNbrSamples( const int subMinRow, int cellRow, int cellCol, double *nbrSamples);

	virtual bool Operator(const CellCoord &coord, bool operFlag);

    int getBlockRowIndexByCoord(double y){  return (_glb_extent.maxY - y) / _blockSize; }
    int getBlockColIndexByCoord(double x){return (x - _glb_extent.minX) / _blockSize;}
    int getBlockRowIndexByCellIndex(int iRow,double cellSize){return min(int((iRow+0.5)*cellSize/_blockSize),_blockRows-1);}
    int getBlockColIndexByCellIndex(int iCol,double cellSize){ return min(int((iCol+0.5)*cellSize/_blockSize),_blockCols-1);}
    double getYByCellIndex(int iRow,double cellSize){return double(_glb_extent.maxY - (iRow + 0.5) * cellSize);}
    double getXByCellIndex(int iCol,double cellSize){return double((iCol + 0.5) * cellSize + _glb_extent.minX);}

    int getBlockCols(){return _blockCols;}
    int getBlockRows(){return _blockRows;}
    float getCellSize(){return _cellSize;}
    double getGlobalMinX(){return _glb_extent.minX;}
    double getGlobalMaxX(){return _glb_extent.maxX;}
    double getGlobalMinY(){return _glb_extent.minY;}
    double getGlobalMaxY(){return _glb_extent.maxY;}
    RasterLayer<int>* getMaskLayer(){return _pMaskLayer;}

    inline double getMinDistanceToBlockBound(double x,double y);
    int getNbrPoints(){return _nbrPoints;}
    void initIdwLayerGlobalInfo(RasterLayer<double>& layerD, char** pSpatialRefWkt);
private:
    int getBlockRowIndexByCellIndex(int iRow){return (iRow+0.5)*_cellSize/_blockSize;}
    int getBlockColIndexByCellIndex(int iCol){return (iCol+0.5)*_cellSize/_blockSize;}
    
    double getYByCellIndex(int iRow){return double(_glb_extent.maxY - (iRow + 0.5) * _cellSize);}
    double getXByCellIndex(int iCol){return double((iCol + 0.5) * _cellSize + _glb_extent.minX);}

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
	double _idw_buffer;
	extent_info _glb_extent;	//�������Ĳ�ֵդ��ͼ��ȫ����Χ�������������ݷ�Χ����ȷ��
	extent_info _sub_extent;	//������դ��Χ
	int _sample_nums;	//ȫ������������
	double _blockSize;
	vector<SampleBlock> _pSampleBlocks;	//�Դ������,������֯������;ÿ�����̶�����ȫ���飬����������Ҫ�ͷţ�
    int _blockRows;
    int _blockCols;
	RasterLayer<double> *_pIDWLayer;
	RasterLayer<int> *_pMaskLayer;


};

#endif