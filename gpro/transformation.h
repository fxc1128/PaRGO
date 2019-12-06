#ifndef Transformation_H
#define Transformation_H

/***************************************************************************
* transformation.h
*
* Project: GPRO, v 2.0
* Purpose: Header file for class GPRO::Transformation
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2017. Ai Beibei
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/
#include <vector>
#include "basicTypes.h"
#include "cellSpace.h"
#include "application.h"
#include "computLayer.h"
#include <iostream>

using namespace std;

namespace GPRO
{
	template <class elemType>
	class Transformation
	{
	public:
		Transformation();
		Transformation( ComputLayer<elemType>* pLayer );	//�̳����Զ���ʵ��ʱ���ô�����
		Transformation( elemType load1, elemType load2, ComputLayer<elemType>* pLayer );
		//������ݸ���ʵ����һ��Ĭ�ϰ汾���������ݷֲ������ȶ�����ֲ����ȵ����

		virtual ~Transformation(){}	//��������Ϊ���࣬���д��������;��ֹ����ָ��ָ����������ͷŲ���

		virtual bool isTermination(bool isIter=false) {return isIter;}	//�̳���ֱ�Ӹİ�termination��Աֵ
		bool Configure(ComputLayer<elemType>* pLayer, bool isCommunication);
		bool paramInit();
		virtual bool Operator(const CellCoord &coord);
		bool run();

	private:
		elemType minLoad;
		elemType maxLoad;
		int _myRank;

		double _noData;	//����ͼ��Ŀ�ֵ
		int _computGrain;
		CoordBR _dataMBR;	//����ͼ������к�
		CoordBR _dataWorkBR;	//�����̴���Ŀ��Ӧ������ͼ�㷶Χ
		//int maxDataRow;
		//int maxDataCol;
		//int glbBeginRow;
		//int glbBeginCol;

	protected:
		//vector<RasterLayer<elemType>* > CommVec;	//��������Ϊ������⣬�����Ҫ
		ComputLayer<elemType>* _pComptLayer;
		//vector<RasterLayer<elemType>* > _pDataLayersV;	//���������Ϊ��ʵ��һЩĬ��ǿ�Ȱ汾
		CoordBR* _pWorkBR;
		int Termination;
	};
};

template <class elemType>
inline GPRO::Transformation<elemType>::
Transformation()
	:minLoad(0),
	maxLoad(0),
	_pComptLayer(NULL),
	_pWorkBR(NULL),
	Termination(1)
{
}


template <class elemType>
inline GPRO::Transformation<elemType>::
Transformation( ComputLayer<elemType>* pLayer )	//�̳����Զ���ʵ��ʱ���ô�����
	:minLoad(0),
	maxLoad(0),
	_pComptLayer(pLayer),
	_pWorkBR(NULL),
	Termination(1)
{
	Configure(pLayer,false);
	_myRank = pLayer->_pMetaData->myrank;
}


template <class elemType>
inline GPRO::Transformation<elemType>::
Transformation( elemType load1, elemType load2, ComputLayer<elemType>* pLayer )	//������ݸ���ָ��Ĭ�ϰ汾
	:minLoad(load1),
	maxLoad(load2),
	_pComptLayer(pLayer),
	_pWorkBR(NULL),
	Termination(1)
{
	Configure(pLayer,false);
	_myRank = pLayer->_pMetaData->myrank;
}


template <class elemType>
bool GPRO::Transformation<elemType>::
Configure(ComputLayer<elemType>* pLayer, bool isCommunication)
{
	//Ŀǰ�ܼ��ԣ��Ժ���ܻ����ͨ�ż��������ݳ�Ա
	if(_pWorkBR == NULL)
	{
		_pWorkBR = &pLayer->_pMetaData->_localworkBR;
	}
	if(isCommunication)
	{
		//do sth.
	}

	return true;
}

template <class elemType>
bool GPRO::Transformation<elemType>::
paramInit()
{
	//��Ҫ��private���ݳ�Աһ���Գ�ʼ����operator�����л��õ�
	if( _pComptLayer->_pDataLayers.empty() ){
		cerr<<"Datalayers used for calculating compute layer should not be null."<<endl;
		return false;
	}
	_noData = _pComptLayer->_pDataLayers[0]->_pMetaData->noData;
	_computGrain = (int)(_pComptLayer->_pMetaData->cellSize / _pComptLayer->_pDataLayers[0]->_pMetaData->cellSize);
	_dataMBR = _pComptLayer->_pDataLayers[0]->_pMetaData->_MBR;
	_dataWorkBR = _pComptLayer->_pDataLayers[0]->_pMetaData->_localworkBR;

	return true;
}

//������ݸ���ʵ����һ��Ĭ�ϰ汾���������ݷֲ������ȶ�����ֲ����ȵ����;�Ҽ���ֲ��ǿ�ָ���ģ����Ǹ��Ӽ���
template <class elemType>
bool GPRO::Transformation<elemType>::
Operator(const CellCoord &coord)
{
	//���ݸ������طֲ������м���
	int cRow = coord.iRow();
	int cCol = coord.iCol();
	CellSpace<elemType> &computL = *(_pComptLayer->cellSpace());

	if( cRow == _pWorkBR->minIRow() && cCol == _pWorkBR->minICol() ){
		cout<<"Transformation operator() function called."<<endl;
		//����һ�����ݳ�Ա��ʼ�������������з�Χ��һ���Գ�ʼ��
		if( minLoad == maxLoad ){
			cerr<<"The load is balanced. No need to use compute layer."<<endl;
			return false;
		}
		if( minLoad<0 || maxLoad<0 ){
			cerr<<"The load specified cannot be negative."<<endl;
			return false;
		}

		if( !paramInit() ){
			return false;
		}
	}
	computL[cRow][cCol] = 0.0;	//��ʼ��
	for( typename vector<RasterLayer<elemType>* >::iterator iter = _pComptLayer->_pDataLayers.begin(); iter!=_pComptLayer->_pDataLayers.end(); ++iter ){
		//��ÿ��ͼ��������㣬�ۻ���������ͼ��ֵ
		CellSpace<elemType> &dataL = *((*iter)->cellSpace());	//ģ�������ָ���������Ƿ���ȷ
		for( int dRow = cRow*_computGrain+_dataWorkBR.minIRow(); dRow<(cRow+1)*_computGrain+_dataWorkBR.minIRow(); ++dRow ){
			for( int dCol = cCol*_computGrain+_dataWorkBR.minICol(); dCol<(cCol+1)*_computGrain+_dataWorkBR.minICol(); ++dCol ){
				if( dRow > _dataMBR.maxIRow()||dRow>=dataL.nRows() || dCol > _dataMBR.maxICol() ||dCol>=dataL.nCols()){
					continue;
				}else{
					if( fabs(dataL[dRow][dCol] - _noData)>Eps && fabs(dataL[dRow][dCol] + 9999)>Eps){	//9999��������ǵĲ������ݶ���д�ģ���ʵû��Ҫ�������������⣬�����ڳ�������
						computL[cRow][cCol] += maxLoad;//wyj 2019-11-12 noData��Ӧ��+=minLoad��
					}else{
						computL[cRow][cCol] += minLoad;
					}
				}
			}
		}
	}

	return true;
}

template <class elemType>
bool GPRO::Transformation<elemType>::
run()
{
	if( Application::_programType != MPI_Type && Application::_programType != MPI_OpenMP_Type ){
		cerr<<"not supported yet."<<endl;
	}
	bool flag = true;
	//MPI_Barrier(MPI_COMM_WORLD);
	if( _myRank == 0 ){	//Ŀǰ��֧�ִ������
		int iterNum = 0;	//��������
		int termSum = 1;
		do
		{
			Termination = 1;
			for(int iRow = _pWorkBR->minIRow(); iRow <= _pWorkBR->maxIRow(); iRow++) 
			{
				for(int iCol = _pWorkBR->minICol();	iCol <= _pWorkBR->maxICol(); iCol++) 
				{
					CellCoord coord(iRow, iCol);
					if(!Operator(coord)) 
					{
						cout<<"Operator is not successes!"<<endl;
						flag = false;
						break;
					}
				}
			}
			//MPI_Allreduce(&Termination, &termSum, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
			termSum = Termination;	//���ṩ��Ŀǰ�Ĵ��а汾
		} while (!termSum);
	}
	//MPI_Barrier(MPI_COMM_WORLD);

	return flag;
}


#endif
