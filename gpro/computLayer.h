#ifndef ComputLayer_H
#define ComputLayer_H

/***************************************************************************
* computLayer.h
*
* Project: GPRO, v 2.0
* Purpose: Header file for class GPRO::ComputLayer
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2017. Ai Beibei
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include "rasterLayer.h"
#include <iostream>

using namespace std;
#define Eps 0.0000001

namespace GPRO {
    template<class elemType>
    class ComputLayer : public RasterLayer<elemType> {
    public:
        ComputLayer();
        ComputLayer( const string RasterLayerName = "Untitled" );
        ComputLayer( vector<RasterLayer<elemType> *> dataLayers,
                     const int tmpGrain,
                     const string RasterLayerName = "Untitled" );
        ~ComputLayer();

        int getGrain() { return _comptGrain; }

        void cleanDataLayers();
        vector<RasterLayer<elemType> *> *dataLayers();
        const vector<RasterLayer<elemType> *> *dataLayers() const;    //���ﷵ�ص�ַ��������const,why
        //��һ��const,���η���ֵ���ڶ���const,���Ա����Ϊconst�����˺��������޸����ݳ�Ա

        //bool newMetaData( const MetaData& rhs, int compuSize );
        bool newMetaData( int comptGrain );    //�������࣬�������ݳ�Ա�����
        //bool transformation();
        bool getCompuLoad( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR );
        bool getCompuLoad2( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR );

        bool readComptFile( const char *outputfile );    //��ʵ�����
        bool writeComptFile( const char *outputfile );

    public:
        vector<RasterLayer<elemType> *> _pDataLayers;
    protected:
        int _comptGrain;
    };
};

template<class elemType>
inline GPRO::ComputLayer<elemType>::
ComputLayer()
    :RasterLayer<elemType>() {
}

template<class elemType>
inline GPRO::ComputLayer<elemType>::
ComputLayer( const string RasterLayerName )
    :RasterLayer<elemType>( RasterLayerName ) {
    //_comptGrain = 10;
    //_pDataLayers.push_back(this);
}

template<class elemType>
inline GPRO::ComputLayer<elemType>::
ComputLayer( vector<RasterLayer<elemType> *> dataLayers, const int tmpGrain, const string RasterLayerName )
    : RasterLayer<elemType>( RasterLayerName ),
      _pDataLayers( dataLayers ), _comptGrain( tmpGrain ) {
    //newMetaData(_comptGrain);	//û��������Ϣ���ʶ������޷�new��workBR
}

template<class elemType>
inline GPRO::ComputLayer<elemType>::
~ComputLayer() {
    //������Զ����û��������������ֻ��Ҫ���ͷ�_pDataLayers��Ա����
    cleanDataLayers();
}

template<class elemType>
void GPRO::ComputLayer<elemType>::
cleanDataLayers() {
    //dataLayers�д�ŵ���ͼ���ָ�룬����ֵ�ͷ�vector�����������ͷ���Щָ����ָ��ͼ��
    //��Щͼ����������ʹ�ã�ֱ���������������������������ȥ�ͷ�
    vector<RasterLayer<elemType> *> vTemp;
    vTemp.swap( _pDataLayers );
}

template<class elemType>
inline vector<GPRO::RasterLayer<elemType> *> *GPRO::ComputLayer<elemType>::
dataLayers() {
    return &_pDataLayers;
}

//template <class elemType>
//inline vector<GPRO::RasterLayer<elemType>* >* GPRO::ComputLayer<elemType>::
//	dataLayers() const
//{
//	return &_pDataLayers;
//}

//template <class elemType>
//bool GPRO::ComputLayer<elemType>::
//newMetaData( const MetaData& rhs, int compuSize ){
//		//_pMetaData = new MetaData();	//VS�����ķ�ʽͨ������gc++������ֱ��ʹ�����Ի��ಿ�ֵ����ݳ�Ա
//		RasterLayer<elemType>::_pMetaData = new MetaData();
//		MetaData* &pMetaData = RasterLayer<elemType>::_pMetaData;	//ָ������ã�ֻ��Ϊ�˼���д,���ñ������ͷ�
//
//		pMetaData->cellSize = rhs.cellSize * compuSize;
//		pMetaData->row = rhs._localworkBR.nRows() / compuSize;	//�����Ԫ���ݶ���Ҫ�������Ȼ���
//		pMetaData->row += (rhs._localworkBR.nRows() % compuSize) ? 1 : 0;
//		pMetaData->column = rhs._localworkBR.nCols() / compuSize;
//		pMetaData->column += (rhs._localworkBR.nCols() % compuSize) ? 1 : 0;
//		pMetaData->format = rhs.format;
//		pMetaData->projection = rhs.projection;
//		pMetaData->noData = rhs.noData;
//		pMetaData->myrank = rhs.myrank;
//		//pMetaData->processor_number = rhs.processor_number;
//		pMetaData->processor_number = 0;	//Ŀǰֻ�Ǵ��й���
//		pMetaData->_domDcmpType = rhs._domDcmpType;	//������Ļ��ַ�ʽδ������������ͬ;Ŀǰ�Ǵ��е�
//		SpaceDims sdim(pMetaData->row, pMetaData->column);
//		pMetaData->_glbDims = sdim;
//		if( pMetaData->_domDcmpType == NON_DCMP ){
//			CoordBR _glbWorkBR;
//			//Neighborhood<elemType>* &pNbrhood = RasterLayer<elemType>::_pNbrhood;
//			RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims );	//���ݼ����������Χȥ�����ռ�
//			//����������Ҳֻ�������ݷ�Χ-����Χ���ķ�Χ
//			pMetaData->_localworkBR = _glbWorkBR;
//			//cout<<"comptLayer L113 "<<*(RasterLayer<elemType>::_pNbrhood)<<endl;
//			//int glbBegin = _glbWorkBR.nwCorner().iRow();
//			//int glbEnd = _glbWorkBR.seCorner().iRow();
//			//CellCoord nwCorner(glbBegin + pNbrhood->minIRow(),
//			//	0);
//			//CellCoord seCorner(glbEnd + pNbrhood->maxIRow(),
//			//	pMetaData->_glbDims.nCols() - 1);
//			CellCoord nwCorner(0, 0);
//			CellCoord seCorner(pMetaData->_glbDims.nRows()-1, pMetaData->_glbDims.nCols()-1);
//			CoordBR subMBR(nwCorner, seCorner);
//			pMetaData->_MBR = subMBR;
//			pMetaData->_localdims = pMetaData->_glbDims;
//			//cout<<"comptLayer L127"<<" dcmpType "<<pMetaData->_domDcmpType<<" glbDims "<<pMetaData->_glbDims<<" localDims "<<pMetaData->_localdims<<" MBR "<<pMetaData->_MBR<<endl;
//			//cout<<"comptLayer L128"<<" _localworkBR "<<pMetaData->_localworkBR<<endl;
//		}else{
//			cerr<<"not support computLayer decomposition now."<<endl;
//			return false;
//		}
//
//		pMetaData->dataType = RasterLayer<elemType>::getGDALType();
//
//		for(int i = 0; i < 6; i++)
//		{
//			pMetaData->pTransform[i] = rhs.pTransform[i];
//		}
//		pMetaData->pTransform[0] += rhs._localworkBR.minICol()*rhs.cellSize;//���������Ͻ������ǹ����ռ䷶Χ��ʼ��
//		pMetaData->pTransform[3] -= rhs._localworkBR.minIRow()*rhs.cellSize;
//		pMetaData->pTransform[1] *= compuSize;//�������ϱ�����һ�����ض�Ӧ�ľ��룬�����
//		pMetaData->pTransform[5] *= compuSize;
//
//		//newCellSpace(pMetaData->_localdims,pMetaData->noData); //allocate
//		RasterLayer<elemType>::newCellSpace(pMetaData->_localdims,0); //allocate,������դ��ֵ��ʼ��Ϊ0
//
//		return true;
//}

template<class elemType>
bool GPRO::ComputLayer<elemType>::
newMetaData( int comptGrain ) {
    //Ŀǰ���н���0�������ô˺���
    //_pMetaData = new MetaData();	//VS�����ķ�ʽͨ������gc++������ֱ��ʹ�����Ի��ಿ�ֵ����ݳ�Ա
    if ( _pDataLayers.empty()) {
        return false;
    }
    //_comptGrain = comptGrain;
    const MetaData &rhs = *( _pDataLayers[0]->_pMetaData );

    RasterLayer<elemType>::_pMetaData = new MetaData();    //newԪ����
    MetaData *&pMetaData = RasterLayer<elemType>::_pMetaData;    //ָ������ã�ֻ��Ϊ�˼���д,���ñ������ͷ�

    pMetaData->cellSize = rhs.cellSize * comptGrain;
    pMetaData->row = rhs._localworkBR.nRows() / comptGrain;    //�����Ԫ���ݶ���Ҫ�������Ȼ���
    pMetaData->row += ( rhs._localworkBR.nRows() % comptGrain ) ? 1 : 0;
    pMetaData->column = rhs._localworkBR.nCols() / comptGrain;
    pMetaData->column += ( rhs._localworkBR.nCols() % comptGrain ) ? 1 : 0;
    pMetaData->format = rhs.format;
    pMetaData->projection = rhs.projection;
    pMetaData->noData = rhs.noData;
    pMetaData->myrank = rhs.myrank;
    //pMetaData->processor_number = rhs.processor_number;
    pMetaData->processor_number = 0;    //Ŀǰֻ�Ǵ��й���
    pMetaData->_domDcmpType = rhs._domDcmpType;    //������Ļ��ַ�ʽδ������������ͬ;Ŀǰ�Ǵ��е�,����non_dcmp
    SpaceDims sdim( pMetaData->row, pMetaData->column );
    pMetaData->_glbDims = sdim;
    if ( pMetaData->_domDcmpType == NON_DCMP ) {
        CoordBR _glbWorkBR;
        //Neighborhood<elemType>* &pNbrhood = RasterLayer<elemType>::_pNbrhood;
        RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims );    //���ݼ����������Χȥ�����ռ�
        //����������Ҳֻ�������ݷ�Χ-����Χ���ķ�Χ
        pMetaData->_localworkBR = _glbWorkBR;
        //cout<<"comptLayer L113 "<<*(RasterLayer<elemType>::_pNbrhood)<<endl;
        //int glbBegin = _glbWorkBR.nwCorner().iRow();
        //int glbEnd = _glbWorkBR.seCorner().iRow();
        //CellCoord nwCorner(glbBegin + pNbrhood->minIRow(),
        //	0);
        //CellCoord seCorner(glbEnd + pNbrhood->maxIRow(),
        //	pMetaData->_glbDims.nCols() - 1);
        CellCoord nwCorner( 0, 0 );
        CellCoord seCorner( pMetaData->_glbDims.nRows() - 1, pMetaData->_glbDims.nCols() - 1 );
        CoordBR subMBR( nwCorner, seCorner );
        pMetaData->_MBR = subMBR;
        pMetaData->_localdims = pMetaData->_glbDims;
        //cout<<"comptLayer L127"<<" dcmpType "<<pMetaData->_domDcmpType<<" glbDims "<<pMetaData->_glbDims<<" localDims "<<pMetaData->_localdims<<" MBR "<<pMetaData->_MBR<<endl;
        //cout<<"comptLayer L128"<<" _localworkBR "<<pMetaData->_localworkBR<<endl;
    } else {
        cerr << "not support computLayer parallel construct now." << endl;
        return false;
    }

    pMetaData->dataType = RasterLayer<elemType>::getGDALType();

    for ( int i = 0; i < 6; i++ ) {
        pMetaData->pTransform[i] = rhs.pTransform[i];
    }
    pMetaData->pTransform[0] += rhs._localworkBR.minICol() * rhs.cellSize;//���������Ͻ������ǹ����ռ䷶Χ��ʼ��
    pMetaData->pTransform[3] -= rhs._localworkBR.minIRow() * rhs.cellSize;
    pMetaData->pTransform[1] *= comptGrain;//�������ϱ�����һ�����ض�Ӧ�ľ��룬�����
    pMetaData->pTransform[5] *= comptGrain;

    //newCellSpace(pMetaData->_localdims,pMetaData->noData); //allocate
    RasterLayer<elemType>::newCellSpace( pMetaData->_localdims, 0 ); //allocate,������դ��ֵ��ʼ��Ϊ0

    return true;
}


////�û��ӿں������ɿ��ǵ�дһ��transformation�����̳У������û�����֪��Ϣ��������͸��
////��Ǩ����transformation��
//template <class elemType>
//bool GPRO::ComputLayer<elemType>::
//transformation(){
//	//����д����ǿ�Ⱥ�������������ͼ��դ��ֵ
//	//cout<<RasterLayer<elemType>::_pMetaData->_localworkBR<<endl;
//	CoordBR &workBR = RasterLayer<elemType>::_pMetaData->_localworkBR;
//	//cout<<"L159 done"<<workBR.minIRow()<<" "<<workBR.maxIRow()<<" "<<workBR.minICol()<<" "<<workBR.maxICol()<<endl;
//	CellSpace<elemType> &computL = *(RasterLayer<elemType>::_pCellSpace);
//	CellSpace<elemType> &dataL = *(_pDataLayers[0]->cellSpace());
//	double dataNoData = _pDataLayers[0]->_pMetaData->noData;
//	int computSize = RasterLayer<elemType>::_pMetaData->cellSize / _pDataLayers[0]->_pMetaData->cellSize;
//	//������ͼ�㲻���ڿ�ֵ��ֻ��0ֵ
//	//cout<<dataNoData<<endl;
//	//dataNoData = -9999;
//	int maxDRow = _pDataLayers[0]->_pMetaData->row;
//	int maxDCol = _pDataLayers[0]->_pMetaData->column;
//	int glbBeginRow = _pDataLayers[0]->_pMetaData->_localworkBR.minIRow();
//	int glbBeginCol = _pDataLayers[0]->_pMetaData->_localworkBR.minICol();
//	//cout<<"L171 done"<<maxDRow<<" "<<maxDCol<<" "<<glbBeginRow<<" "<<glbBeginCol<<endl;
//	for(int cRow = workBR.minIRow(); cRow <= workBR.maxIRow(); cRow++) 
//	{
//		for(int cCol = workBR.minICol(); cCol <= workBR.maxICol(); cCol++) 
//		{
//			//cout<<computL[cRow][cCol]<<endl;
//			//�Ա�������cell����Ч������Ϊ���ؼ���
//			for( int dRow = cRow*computSize+glbBeginRow; dRow < (cRow+1)*computSize+glbBeginRow; ++dRow ){
//				for( int dCol = cCol*computSize+glbBeginCol; dCol < (cCol+1)*computSize+glbBeginCol; ++dCol ){
//					if( dRow > maxDRow-1 || dCol > maxDCol-1 ){
//						continue;
//					}else{
//						if( fabs(dataL[dRow][dCol] - dataNoData)>Eps && fabs(dataL[dRow][dCol] + 9999)>Eps){	//9999��������ǵĲ������ݶ���д�ģ���ʵû��Ҫ�������������⣬�����ڳ�������
//							computL[cRow][cCol] += 10;
//						}else{
//							computL[cRow][cCol] += 1;
//						}
//					}
//				}
//			}
//		}
//	}
//	return true;
//}

//template <class elemType>
//bool GPRO::ComputLayer<elemType>::
//getCompuLoad( int* pDcmpIdx, DomDcmpType dcmpType, const int computSize, const int nSubSpcs ){
//		vector<CoordBR> vDcmpBR;
//		//���������������ȶ�����������ͼ�����ݣ����������ռ䷶Χ������ͼ�����������Ա�ģ�������ͼ����ʱû�л�����һ�£����ٸ������ȣ�����������ͼ�㣻���ܷ�Χ���ܹ����ռ䷶Χ���϶��룬����һ�£�
//		//������⺯������forѭ��compuLayer��ÿ��դ���������ֵ������compuLayer�ĸ�դ��ֵ
//		//����decomp���󣬵��û��ֺ���valRowDcmp()������compuLayerֵ������з�Χ���֣������vector<CoordBR>& ����
//		//����compuLayerͼ��Ļ��ֽ����ӳ�䵽���ݵĹ����ռ䷶Χ�����ظ���������vDcmpIdx
//		//�ȴ�����⣬�����̸�����metedata��ͨ�Ÿ����ӽ���
//		if( _pDataLayers.empty() ){
//			return false;
//		}
//
//		newMetaData( *(_pDataLayers[0]->_pMetaData), computSize );	//��ʼ���˻���rasterLayer�������ݳ�Ա
//		//cout<<"L207 done"<<RasterLayer<elemType>::_pMetaData->cellSize<<" "<<RasterLayer<elemType>::_pMetaData->row<<" "<<RasterLayer<elemType>::_pMetaData->column<<endl;
//		//comptLayer���ǵķ�Χֻ��metaLayer��glbWorkBR��Χ
//		transformation();	//���computLayer._pCellSpace���ݳ�Ա
//		vector<CoordBR> vComptDcmpBR;
//		DeComposition<elemType> deComp(RasterLayer<elemType>::_pMetaData->_glbDims, *(RasterLayer<elemType>::_pNbrhood));
//		if( dcmpType == ROWWISE_DCMP ){
//			deComp.valRowDcmp( vComptDcmpBR, *this, nSubSpcs);	//��ֵ���֣�����Ҫͼ��Ϊ����;���ֽ���������ô��ظ�vComptDcmpBR
//			_pDataLayers[0]->_pMetaData->_domDcmpType = ROWWISE_DCMP;
//		}else{
//			cerr<<"not support until now."<<endl;
//		}
//		//for( vector<CoordBR>::iterator iter = vComptDcmpBR.begin(); iter!=vComptDcmpBR.end(); ++iter ){
//		//	cout<<*iter<<endl;
//		//}
//		//�����ֽ��ӳ������ݿռ���ӷ�Χ
//		CoordBR _glbWorkBR;
//		Neighborhood<elemType> *pDataNbrhood = _pDataLayers[0]->nbrhood();
//		pDataNbrhood->calcWorkBR( _glbWorkBR, _pDataLayers[0]->_pMetaData->_glbDims );	//����ͼ���ȫ�ֹ����ռ�
//		int subBegin = _glbWorkBR.minIRow(), subEnd = _glbWorkBR.minIRow()-1;
//		int i = 0;
//		for( ;i<nSubSpcs-1; ++i ){
//			subBegin = vComptDcmpBR[i].minIRow()*computSize + _glbWorkBR.minIRow();
//			subEnd = (vComptDcmpBR[i+1].minIRow())*computSize + _glbWorkBR.minIRow()-1;
//			//cout<<i<<" "<<subBegin<<" , "<<subEnd<<endl;
//			CellCoord nwCorner(subBegin, _glbWorkBR.minICol());
//			CellCoord seCorner(subEnd, _glbWorkBR.maxICol());
//			CoordBR subMBR(nwCorner, seCorner);
//			vDcmpBR.push_back(subMBR);
//			pDcmpIdx[4*i] = subBegin;
//			pDcmpIdx[4*i+1] = _glbWorkBR.minICol();
//			pDcmpIdx[4*i+2] = subEnd;
//			pDcmpIdx[4*i+3] = _glbWorkBR.maxICol();
//		}
//		CellCoord nwCorner(subEnd+1, _glbWorkBR.minICol());
//		CellCoord seCorner(_glbWorkBR.maxIRow(), _glbWorkBR.maxICol());
//		CoordBR subMBR(nwCorner, seCorner);
//		vDcmpBR.push_back(subMBR);
//		pDcmpIdx[4*i] = subEnd+1;
//		pDcmpIdx[4*i+1] = _glbWorkBR.minICol();
//		pDcmpIdx[4*i+2] = _glbWorkBR.maxIRow();
//		pDcmpIdx[4*i+3] = _glbWorkBR.maxICol();
//		//for( vector<CoordBR>::iterator iter = vDcmpBR.begin(); iter != vDcmpBR.end(); ++iter ){
//		//	cout<<*iter<<endl;	//it's ok here
//		//}
//		//��ΪĿǰMPI��֧��ͨ���Զ������ͣ���vDcmpBR��ʱʵ���ϲ�û��
//
//		////������η��ʸ�Ԫ���ݼ�դ��ֵ
//		//CoordBR* pWorkBR = &(_pMetaData->_localworkBR);
//		//cout<<pWorkBR->minIRow()<<" "<<pWorkBR->maxIRow()<<" "<<pWorkBR->minICol()<<" "<<pWorkBR->maxICol()<<endl;
//		//for(int iRow = pWorkBR->minIRow(); iRow <= pWorkBR->minIRow()+1; iRow++) 
//		//{
//		//	for(int iCol = pWorkBR->minICol();	iCol <= pWorkBR->minICol()+5; iCol++) 
//		//	{
//		//		CellSpace<double> &demL = *(_pCellSpace);
//		//		cout<<demL[iRow][iCol]<<endl;
//		//	}
//		//}
//
//		cout<<"getCompuLoad"<<endl;
//
//		return true;
//}

template<class elemType>
bool GPRO::ComputLayer<elemType>::
getCompuLoad( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR ) {
    int myRank, process_nums;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &process_nums );

    int *pDcmpIdx = new int[nSubSpcs * 4];    //�洢�黮�ֽ�������㲥
    if ( myRank == 0 ) {    //Ŀǰ���صĻ������������̶������
        vector<CoordBR> vDcmpBR;    //���л��ֿ�,��nSubSpcs��;��ΪMPI��֧��ͨ���Զ������ͣ���ʵ������
        //���������������ȶ�����������ͼ�����ݣ����������ռ䷶Χ������ͼ�����������Ա�ģ�������ͼ����ʱû�л�����һ�£����ٸ������ȣ�����������ͼ�㣻���ܷ�Χ���ܹ����ռ䷶Χ���϶��룬����һ�£�
        //������⺯������forѭ��compuLayer��ÿ��դ���������ֵ������compuLayer�ĸ�դ��ֵ
        //����decomp���󣬵��û��ֺ���valRowDcmp()������compuLayerֵ������з�Χ���֣������vector<CoordBR>& ����
        //����compuLayerͼ��Ļ��ֽ����ӳ�䵽���ݵĹ����ռ䷶Χ�����ظ���������vDcmpIdx
        //�ȴ�����⣬�����̸�����metedata��ͨ�Ÿ����ӽ���

        vector<CoordBR> vComptDcmpBR;
        cout << "hold " << RasterLayer<elemType>::_pMetaData->_glbDims << endl;
        DeComposition<elemType>
            deComp( RasterLayer<elemType>::_pMetaData->_glbDims, *( RasterLayer<elemType>::_pNbrhood ));
        if ( dcmpType == ROWWISE_DCMP ) {
            deComp.valRowDcmp( vComptDcmpBR, *this, nSubSpcs );    //��ֵ���֣�����Ҫͼ��Ϊ����;���ֽ���������ô��ظ�vComptDcmpBR
            //_pDataLayers[0]->_pMetaData->_domDcmpType = ROWWISE_DCMP;	//��ʱ;����û��//0316ɾ������֪���Ƿ�Ӱ��
        } else {
            cerr << "computLayer L388: not support until now." << endl;
        }
        //�����ֽ��ӳ������ݿռ���ӷ�Χ
        CoordBR _glbWorkBR;
        Neighborhood<elemType> *pDataNbrhood = _pDataLayers[0]->nbrhood();
        pDataNbrhood->calcWorkBR( _glbWorkBR, _pDataLayers[0]->_pMetaData->_glbDims );    //����ͼ���ȫ�ֹ����ռ�
        int subBegin = _glbWorkBR.minIRow(), subEnd = _glbWorkBR.minIRow() - 1;
        int i = 0;
        for ( ; i < nSubSpcs - 1; ++i ) {
            subBegin = vComptDcmpBR[i].minIRow() * _comptGrain + _glbWorkBR.minIRow();
            subEnd = ( vComptDcmpBR[i + 1].minIRow()) * _comptGrain + _glbWorkBR.minIRow() - 1;
            //cout<<i<<" "<<subBegin<<" , "<<subEnd<<endl;
            CellCoord nwCorner( subBegin, _glbWorkBR.minICol());
            CellCoord seCorner( subEnd, _glbWorkBR.maxICol());
            CoordBR subMBR( nwCorner, seCorner );
            vDcmpBR.push_back( subMBR );
            pDcmpIdx[4 * i] = subBegin;
            pDcmpIdx[4 * i + 1] = _glbWorkBR.minICol();
            pDcmpIdx[4 * i + 2] = subEnd;
            pDcmpIdx[4 * i + 3] = _glbWorkBR.maxICol();
        }
        //�洢ת��������
        CellCoord nwCorner( subEnd + 1, _glbWorkBR.minICol());
        CellCoord seCorner( _glbWorkBR.maxIRow(), _glbWorkBR.maxICol());
        CoordBR subMBR( nwCorner, seCorner );
        vDcmpBR.push_back( subMBR );
        pDcmpIdx[4 * i] = subEnd + 1;
        pDcmpIdx[4 * i + 1] = _glbWorkBR.minICol();
        pDcmpIdx[4 * i + 2] = _glbWorkBR.maxIRow();
        pDcmpIdx[4 * i + 3] = _glbWorkBR.maxICol();
    }

    MPI_Barrier( MPI_COMM_WORLD );    //Ŀǰ����У������̻�ȥ��ɼ�����ͼ��ļ��㣬������������Ⱥ򣬼�����ɺ��ȡ���Ե�subWorkBR����
    //Ŀǰ�Ǵ��У�����Ҫ�㲥���������̶�������һ��or���л��֣������޸�
    MPI_Bcast( pDcmpIdx, process_nums * 4, MPI_INT, 0, MPI_COMM_WORLD );
    CellCoord nwCorner2( pDcmpIdx[4 * myRank], pDcmpIdx[4 * myRank + 1] );
    CellCoord seCorner2( pDcmpIdx[4 * myRank + 2], pDcmpIdx[4 * myRank + 3] );
    CoordBR tmpWorkBR( nwCorner2, seCorner2 );
    subWorkBR = tmpWorkBR;    //�Ѳ��ԣ�������ֵû����
    //cout<<"tmpWorkBR "<<tmpWorkBR<<" subWorkBR "<<subWorkBR<<" shoule be the same"<<endl;
    MPI_Barrier( MPI_COMM_WORLD );
    delete[]pDcmpIdx;

    return true;
}

template<class elemType>
bool GPRO::ComputLayer<elemType>::
getCompuLoad2( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR ) {
    int myRank, process_nums;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &process_nums );

    int *pDcmpIdx = new int[nSubSpcs * 4];    //�洢�黮�ֽ�������㲥
    if ( myRank == 0 ) {    //Ŀǰ���صĻ������������̶������
        vector<CoordBR> vDcmpBR;    //���л��ֿ�,��nSubSpcs��;��ΪMPI��֧��ͨ���Զ������ͣ���ʵ������
        //���������������ȶ�����������ͼ�����ݣ����������ռ䷶Χ������ͼ�����������Ա�ģ�������ͼ����ʱû�л�����һ�£����ٸ������ȣ�����������ͼ�㣻���ܷ�Χ���ܹ����ռ䷶Χ���϶��룬����һ�£�
        //������⺯������forѭ��compuLayer��ÿ��դ���������ֵ������compuLayer�ĸ�դ��ֵ
        //����decomp���󣬵��û��ֺ���valRowDcmp()������compuLayerֵ������з�Χ���֣������vector<CoordBR>& ����
        //����compuLayerͼ��Ļ��ֽ����ӳ�䵽���ݵĹ����ռ䷶Χ�����ظ���������vDcmpIdx
        //�ȴ�����⣬�����̸�����metedata��ͨ�Ÿ����ӽ���

        vector<CoordBR> vComptDcmpBR;
        DeComposition<elemType>
            deComp( RasterLayer<elemType>::_pMetaData->_glbDims, *( RasterLayer<elemType>::_pNbrhood ));
        if ( dcmpType == ROWWISE_DCMP ) {
            deComp.valRowDcmp( vComptDcmpBR, *this, nSubSpcs );    //��ֵ���֣�����Ҫͼ��Ϊ����;���ֽ���������ô��ظ�vComptDcmpBR
            //_pDataLayers[0]->_pMetaData->_domDcmpType = ROWWISE_DCMP;	//��ʱ;����û��//0316ɾ������֪���Ƿ�Ӱ��
        } else {
            cerr << "computLayer L388: not support until now." << endl;
        }
        //�����ֽ��ӳ������ݿռ���ӷ�Χ
        CoordBR _glbWorkBR;
        Neighborhood<elemType> *pDataNbrhood = _pDataLayers[0]->nbrhood();
        pDataNbrhood->calcWorkBR( _glbWorkBR, _pDataLayers[0]->_pMetaData->_glbDims );    //����ͼ���ȫ�ֹ����ռ�
        int subBegin = _glbWorkBR.minIRow(), subEnd = _glbWorkBR.minIRow() - 1;
        int i = 0;
        for ( ; i < nSubSpcs - 1; ++i ) {
            subBegin = vComptDcmpBR[i].minIRow() * _comptGrain + _glbWorkBR.minIRow();
            subEnd = ( vComptDcmpBR[i + 1].minIRow()) * _comptGrain + _glbWorkBR.minIRow() - 1;
            //cout<<i<<" "<<subBegin<<" , "<<subEnd<<endl;
            CellCoord nwCorner( subBegin, _glbWorkBR.minICol());
            CellCoord seCorner( subEnd, _glbWorkBR.maxICol());
            CoordBR subMBR( nwCorner, seCorner );
            vDcmpBR.push_back( subMBR );
            pDcmpIdx[4 * i] = subBegin;
            pDcmpIdx[4 * i + 1] = _glbWorkBR.minICol();
            pDcmpIdx[4 * i + 2] = subEnd;
            pDcmpIdx[4 * i + 3] = _glbWorkBR.maxICol();
        }
        //�洢ת��������
        CellCoord nwCorner( subEnd + 1, _glbWorkBR.minICol());
        CellCoord seCorner( _glbWorkBR.maxIRow(), _glbWorkBR.maxICol());
        CoordBR subMBR( nwCorner, seCorner );
        vDcmpBR.push_back( subMBR );
        pDcmpIdx[4 * i] = subEnd + 1;
        pDcmpIdx[4 * i + 1] = _glbWorkBR.minICol();
        pDcmpIdx[4 * i + 2] = _glbWorkBR.maxIRow();
        pDcmpIdx[4 * i + 3] = _glbWorkBR.maxICol();
    }

    MPI_Barrier( MPI_COMM_WORLD );    //Ŀǰ����У������̻�ȥ��ɼ�����ͼ��ļ��㣬������������Ⱥ򣬼�����ɺ��ȡ���Ե�subWorkBR����
    //Ŀǰ�Ǵ��У�����Ҫ�㲥���������̶�������һ��or���л��֣������޸�
    MPI_Bcast( pDcmpIdx, process_nums * 4, MPI_INT, 0, MPI_COMM_WORLD );
    CellCoord nwCorner2( pDcmpIdx[4 * myRank], pDcmpIdx[4 * myRank + 1] );
    CellCoord seCorner2( pDcmpIdx[4 * myRank + 2], pDcmpIdx[4 * myRank + 3] );
    CoordBR tmpWorkBR( nwCorner2, seCorner2 );
    subWorkBR = tmpWorkBR;    //�Ѳ��ԣ�������ֵû����
    //cout<<"tmpWorkBR "<<tmpWorkBR<<" subWorkBR "<<subWorkBR<<" shoule be the same"<<endl;
    MPI_Barrier( MPI_COMM_WORLD );
    delete[]pDcmpIdx;

    return true;
}

template<class elemType>
bool GPRO::ComputLayer<elemType>::
writeComptFile( const char *outputfile ) {
    //Ŀǰ��֧�ִ���д��
    GDALAllRegister();

    if ( !RasterLayer<elemType>::createFile( outputfile )) {
        cout << "create file is not correct!" << endl;
        MPI_Finalize();
    }
    GDALDataset *poDataset = NULL;
    poDataset = (GDALDataset *) GDALOpen( outputfile, GA_Update );
    if ( poDataset == NULL /*����Ƿ��������ļ�*/) {
        //do something
        cout << "data file is not open correct" << endl;
        exit( 1 );
    }
    GDALRasterBand *poBanddest = poDataset->GetRasterBand( 1 );
    if ( poBanddest == NULL ) {
        //do something
        cout << "poBanddest is NULL" << endl;
        exit( 1 );
    }
    poBanddest->SetNoDataValue( RasterLayer<elemType>::_pMetaData->noData );

    if ( RasterLayer<elemType>::_pMetaData->myrank == 0 ) {
        poBanddest->RasterIO( GF_Write,
                              0,
                              0,
                              RasterLayer<elemType>::_pMetaData->_glbDims.nCols(),
                              RasterLayer<elemType>::_pMetaData->_glbDims.nRows(),
                              RasterLayer<elemType>::_pCellSpace->_matrix,
                              RasterLayer<elemType>::_pMetaData->_glbDims.nCols(),
                              RasterLayer<elemType>::_pMetaData->_glbDims.nRows(),
                              RasterLayer<elemType>::_pMetaData->dataType,
                              0,
                              0 );
    }

    if ( poDataset != NULL ) {
        GDALClose((GDALDatasetH) poDataset );
        poDataset = NULL;
    }

    return true;
}

#endif