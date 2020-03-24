#ifndef ComputeLayer_H
#define ComputeLayer_H

/***************************************************************************
* computLayer.h
*
* Project: GPRO, v 2.0
* Purpose: Header file for class GPRO::ComputeLayer
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
#include "utility.h"
#include <iostream>

using namespace std;
#define Eps 0.0000001

inline bool ifDecomposeBySpace(const string& arg) {
    if (StringMatch(arg, "space")) {
        return true;
    }
    if (StringMatch(arg, "compute")) {
        return false;
    }
    return true;
}

namespace GPRO {
    template<class elemType>
    class ComputeLayer : public RasterLayer<elemType> {
    public:
        ComputeLayer();
		ComputeLayer( const string RasterLayerName = "Untitled" );
		ComputeLayer( RasterLayer<elemType> * dataLayers,
			         const int tmpGrain,
			         const string RasterLayerName = "Untitled" );
        ComputeLayer( vector<RasterLayer<elemType> *> dataLayers,
                     const int tmpGrain,
                     const string RasterLayerName = "Untitled" );
        
        ~ComputeLayer();
        int getGrain() { return _comptGrain; }
        void setComputGrain(int comptGrain) { _comptGrain = comptGrain;}
        ComputeLayerType getComputeLayerType() {return _computeLayerType;}
        void setComputeLayerType(ComputeLayerType type) {_computeLayerType=type;}
        

        bool init(const RasterLayer<elemType>* pDataLayer, const char* neighborFile,int comptGrain=1);
        bool init(int comptGrain=1);

        void cleanDataLayers();
        vector<RasterLayer<elemType> *> *dataLayers();
        const vector<RasterLayer<elemType> *> *dataLayers() const;
        void addRasterLayer(RasterLayer<elemType> &dataLayer);
        void addRasterLayers(vector<RasterLayer<elemType> *> dataLayers);

        bool getCompuLoad( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR );

        bool readComputeFile( const char *loadFile, const char* nbrFile);
        bool writeComputeFile( const char *outputfile );
    public:
        vector<RasterLayer<elemType> *> _vDataLayers;
    private:
        double _comptGrain;
        ComputeLayerType _computeLayerType;
    };
};

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer()
    :RasterLayer<elemType>() {}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer( const string RasterLayerName )
    :RasterLayer<elemType>( RasterLayerName ) {}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer( RasterLayer<elemType> * dataLayer, const int tmpGrain, const string RasterLayerName )
	: RasterLayer<elemType>( RasterLayerName ), _comptGrain( tmpGrain ) {
		_vDataLayers.push_back(dataLayer);
    //GPRO::RasterLayer<elemType>::_pMetaData=dataLayer->metaData(); // wyj 12-19 ̫�ң���ʱ��һ�ʿ��ܲ����ܡ����ֲ��ܣ������˼ҵ����õ��¶���������
		//newMetaData(_comptGrain);	//û��������Ϣ���ʶ������޷�new��workBR ///wyj ����Ҫ���ֶ�����newMetaData...
}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer( vector<RasterLayer<elemType> *> dataLayers, const int tmpGrain, const string RasterLayerName )
    : RasterLayer<elemType>( RasterLayerName ),
      _vDataLayers( dataLayers ), _comptGrain( tmpGrain ) {
    //newMetaData(_comptGrain);	//û��������Ϣ���ʶ������޷�new��workBR ///wyj ����Ҫ���ֶ�����newMetaData...
}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
~ComputeLayer() {
    //������Զ����û��������������ֻ��Ҫ���ͷ�_pDataLayers��Ա����
    cleanDataLayers();
}

template<class elemType>
void GPRO::ComputeLayer<elemType>::
cleanDataLayers() {
    //dataLayers�д�ŵ���ͼ���ָ�룬����ֵ�ͷ�vector�����������ͷ���Щָ����ָ��ͼ��
    //��Щͼ����������ʹ�ã�ֱ���������������������������ȥ�ͷ�
    vector<RasterLayer<elemType> *> vTemp;
    vTemp.swap( _vDataLayers );
}

template<class elemType>
inline vector<GPRO::RasterLayer<elemType> *> *GPRO::ComputeLayer<elemType>::
dataLayers() {
    return &_vDataLayers;
}

template<class elemType>
inline void GPRO::ComputeLayer<elemType>::
addRasterLayer(RasterLayer<elemType> &dataLayer) {
    _vDataLayers.push_back(&dataLayer);
}

template<class elemType>
inline void GPRO::ComputeLayer<elemType>::
addRasterLayers(vector<RasterLayer<elemType> *> dataLayers) {
    _vDataLayers=dataLayers;
}

template<class elemType>
bool GPRO::ComputeLayer<elemType>::
init(int comptGrain) {
    // It is a SERIAL function. Only invoked by process 0.
    // Implicitly using members from base class is valid in Visual Studio but not allowed in gc++. i.e. _pMetaData = new MetaData() arises an error.
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if ( pDataLayer == nullptr ) {
        return false;
    }
    
    RasterLayer<elemType>::readNeighborhood(neighborFile);

    const MetaData &rhs = *( pDataLayer->_pMetaData );

    RasterLayer<elemType>::_pMetaData = new MetaData();
    MetaData *&pMetaData = RasterLayer<elemType>::_pMetaData; //Pointer as a reference. No need to delete/free.

    pMetaData->cellSize = rhs.cellSize * comptGrain;
    pMetaData->row = rhs._localworkBR.nRows() / comptGrain;
    pMetaData->row += ( rhs._localworkBR.nRows() % comptGrain ) ? 1 : 0;
    pMetaData->column = rhs._localworkBR.nCols() / comptGrain;
    pMetaData->column += ( rhs._localworkBR.nCols() % comptGrain ) ? 1 : 0;
    pMetaData->format = rhs.format;
    pMetaData->projection = rhs.projection;
    pMetaData->noData = rhs.noData;
    pMetaData->myrank = rhs.myrank;
    pMetaData->processor_number = myRank; //Ŀǰֻ�Ǵ��й���
    pMetaData->_domDcmpType = rhs._domDcmpType; //������Ļ��ַ�ʽδ������������ͬ;Ŀǰ�Ǵ��е�,����non_dcmp
    SpaceDims sdim( pMetaData->row, pMetaData->column );
    pMetaData->_glbDims = sdim
    if ( pMetaData->_domDcmpType == NON_DCMP ) {
        CoordBR _glbWorkBR;
        RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims ); //���ݼ����������Χȥ�����ռ�
        pMetaData->_localworkBR = _glbWorkBR;
        CellCoord nwCorner( 0, 0 );
        CellCoord seCorner( pMetaData->_glbDims.nRows() - 1, pMetaData->_glbDims.nCols() - 1 );
        CoordBR subMBR( nwCorner, seCorner );
        pMetaData->_MBR = subMBR;
        pMetaData->_localdims = pMetaData->_glbDims;
    } else {
        cerr << "not support computeLayer parallel construct now." << endl;
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

    RasterLayer<elemType>::newCellSpace( pMetaData->_localdims, 0 ); //allocate,������դ��ֵ��ʼ��Ϊ0

    return true;
}

template<class elemType>
bool GPRO::ComputeLayer<elemType>::
init(const RasterLayer<elemType>* pDataLayer, const char* neighborFile,int comptGrain) {
    // It is a SERIAL function. Only invoked by process 0.
    // Implicitly using members from base class is valid in Visual Studio but not allowed in gc++. i.e. _pMetaData = new MetaData() arises an error.
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if ( pDataLayer == nullptr ) {
        return false;
    }
    
    RasterLayer<elemType>::readNeighborhood(neighborFile);

    const MetaData &rhs = *( pDataLayer->_pMetaData );

    RasterLayer<elemType>::_pMetaData = new MetaData();
    MetaData *&pMetaData = RasterLayer<elemType>::_pMetaData; //Pointer as a reference. No need to delete/free.

    pMetaData->cellSize = rhs.cellSize * comptGrain;
    pMetaData->row = rhs._localworkBR.nRows() / comptGrain;
    pMetaData->row += ( rhs._localworkBR.nRows() % comptGrain ) ? 1 : 0;
    pMetaData->column = rhs._localworkBR.nCols() / comptGrain;
    pMetaData->column += ( rhs._localworkBR.nCols() % comptGrain ) ? 1 : 0;
    pMetaData->format = rhs.format;
    pMetaData->projection = rhs.projection;
    pMetaData->noData = rhs.noData;
    pMetaData->myrank = rhs.myrank;
    pMetaData->processor_number = myRank; //Ŀǰֻ�Ǵ��й���
    pMetaData->_domDcmpType = rhs._domDcmpType; //������Ļ��ַ�ʽδ������������ͬ;Ŀǰ�Ǵ��е�,����non_dcmp
    SpaceDims sdim( pMetaData->row, pMetaData->column );
    pMetaData->_glbDims = sdim;
    if ( pMetaData->_domDcmpType == NON_DCMP ) {
        CoordBR _glbWorkBR;
        RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims ); //���ݼ����������Χȥ�����ռ�
        pMetaData->_localworkBR = _glbWorkBR;
        CellCoord nwCorner( 0, 0 );
        CellCoord seCorner( pMetaData->_glbDims.nRows() - 1, pMetaData->_glbDims.nCols() - 1 );
        CoordBR subMBR( nwCorner, seCorner );
        pMetaData->_MBR = subMBR;
        pMetaData->_localdims = pMetaData->_glbDims;
    } else {
        cerr << "not support computeLayer parallel construct now." << endl;
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

    RasterLayer<elemType>::newCellSpace( pMetaData->_localdims, 0 ); //allocate,������դ��ֵ��ʼ��Ϊ0

    return true;
}

template<class elemType>
bool GPRO::ComputeLayer<elemType>::
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
        DeComposition<elemType> deComp( RasterLayer<elemType>::_pMetaData->_glbDims, *( RasterLayer<elemType>::_pNbrhood ));
        if ( dcmpType == ROWWISE_DCMP ) {
            deComp.valRowDcmp( vComptDcmpBR, *this, nSubSpcs );    //��ֵ���֣�����Ҫͼ��Ϊ����;���ֽ���������ô��ظ�vComptDcmpBR
            //_vDataLayers[0]->_pMetaData->_domDcmpType = ROWWISE_DCMP;	//��ʱ;����û��//0316ɾ������֪���Ƿ�Ӱ��
        } else {
            cerr << "computLayer L388: not support until now." << endl;
        }
        //�����ֽ��ӳ������ݿռ���ӷ�Χ
        CoordBR glbWorkBR;
        Neighborhood<elemType> *pDataNbrhood = _vDataLayers[0]->nbrhood();
        pDataNbrhood->calcWorkBR( glbWorkBR, _vDataLayers[0]->_pMetaData->_glbDims );    //����ͼ���ȫ�ֹ����ռ�
        int subBegin = glbWorkBR.minIRow(), subEnd = glbWorkBR.minIRow() - 1;
        int outputRows=_vDataLayers[0]->_pMetaData->_glbDims.nRows();
        int loadFileRows=this->metaData()->_glbDims.nRows();
        _comptGrain=outputRows/double(loadFileRows);//wyj 2019-12-6 ������һ��...���������ͼʱ��ͼԴ���Ȳ�ͳһ�����⣨����ռ���⻮�ֵĸ���ͼ�ǵȱ����ģ������㸺�ػ��ֵĸ���ͼ��1:10�ģ�
        int i = 0;
        for ( ; i < nSubSpcs - 1; ++i ) {
            subBegin = vComptDcmpBR[i].minIRow() * _comptGrain + glbWorkBR.minIRow();
            subEnd = ( vComptDcmpBR[i + 1].minIRow()) * _comptGrain + glbWorkBR.minIRow() - 1;
            //cout<<i<<" "<<subBegin<<" , "<<subEnd<<endl;
            CellCoord nwCorner( subBegin, glbWorkBR.minICol());
            CellCoord seCorner( subEnd, glbWorkBR.maxICol());
            CoordBR subMBR( nwCorner, seCorner );
            vDcmpBR.push_back( subMBR );
            pDcmpIdx[4 * i] = subBegin;
            pDcmpIdx[4 * i + 1] = glbWorkBR.minICol();
            pDcmpIdx[4 * i + 2] = subEnd;
            pDcmpIdx[4 * i + 3] = glbWorkBR.maxICol();
        }
        //�洢ת��������
        CellCoord nwCorner( subEnd + 1, glbWorkBR.minICol());
        CellCoord seCorner( glbWorkBR.maxIRow(), glbWorkBR.maxICol());
        CoordBR subMBR( nwCorner, seCorner );
        vDcmpBR.push_back( subMBR );
        pDcmpIdx[4 * i] = subEnd + 1;
        pDcmpIdx[4 * i + 1] = glbWorkBR.minICol();
        pDcmpIdx[4 * i + 2] = glbWorkBR.maxIRow();
        pDcmpIdx[4 * i + 3] = glbWorkBR.maxICol();
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
bool GPRO::ComputeLayer<elemType>::
readComputeFile(const char *loadFile, const char* nbrFile ) {
    RasterLayer<elemType> loadLayer("loadLayer");
    loadLayer.readNeighborhood(nbrFile);
    loadLayer.readFile(loadFile);
    loadLayer.newCellSpace(loadLayer.metaData()->_localdims);
    addRasterLayer(loadLayer);
    return true;
}

template<class elemType>
bool GPRO::ComputeLayer<elemType>::
writeComputeFile( const char *outputfile ) {
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