#ifndef RUNOPERATOR_H
#define RUNOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include "transformation.h"
#include "utility.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.0000001

class FCMOperator : public RasterOperator<double> 
{
public:
	FCMOperator()
		:RasterOperator<double>(),
		 _iterNum(0), flag(true),
		block(0),nodataNums(0), dFlag(true),subval(0.0),totval(0.0),oldtval(0.0),partitionCoef(0.0),entropy(0),totpartitionCoef(0),totentropy(0)
		{}
	//�����ĳ�ʼ��

	~FCMOperator();

	void initialization(int iNum, int cNum, int maxIter,double toler, double m);
	void inputLayer(vector<RasterLayer<double> *> layerD);
	void fcmLayer(RasterLayer<double> &layerD);
	void degLayer(vector<RasterLayer<double> *> layerD);
	
	virtual bool isTermination();

	void createRandomIdx( int nums, int range, int* randomIdx );	//��range��Χ�ڲ���nums���������randomIdx����
	void fnDistance(int curRow, int curCol, double* pInputVal); //�������
	void InitDegree(int curRow, int curCol);//����������
	void initRandomClusterCenters(double *clusterCenters);
	void assignMaxMembershipDegrees();
	virtual bool Operator(const CellCoord &coord, bool operFlag);


private:
	int _cellSize;
	int _xSize, _ySize;	//��ǰ������DEM�������
	int _nRows, _nCols; //����ͼ����������������
	double _noData;
	int _rank;
	int _iterNum;//��������
	bool flag;
	double starttime, endtime;
	double tmpSumTime1, tmpSumTime2;
protected:
	int clusterNum; //������Ŀ
	int maxIteration; //����������
	double tolerance;//������ֵ
	double weight;//��Ȩָ��
	int imageNum;//�����Ӱ����Ŀ

	int block;	//�����̷ֿ�������
	int nodataNums;
	bool dFlag;		//��������Ƿ��һ�ε���
	double* centerVal;	//������������ֵ��һά����
	int* centerIndex;	//��������λ������ clusterNum * 1��С���飬xy����ӳ���һά
	double*** degree;	//���������飬��ά
	double*** dist;	//��Ÿ�������������ĵľ���
	double subval;	//��ű��ε�����Լ������ֵ
	double totval;	//�����̻��ܸ����̣��õ���ֵ��
	double oldtval;	//�ϴε�������ֵ��
	double *sumNumerator;	//�������
	double *sumDenominator;		//��ĸ���
	double *totNumerator;	//���ӹ鲢
	double *totDenominator;		//��ĸ�鲢
	double partitionCoef;	//�ָ�ϵ��
	double entropy;	//��
	double totpartitionCoef;	//��Լ��ķָ�ϵ��
	double totentropy;	//��Լ�����

	vector<RasterLayer<double> *> _vInputLayer;
	RasterLayer<double> *_pFCMLayer;
	vector<RasterLayer<double> *> _vDegLayer;

	Neighborhood<double> *_pDEMNbrhood;

};

#endif