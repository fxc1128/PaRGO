/**
 * \file basicTypes
 * \author Qin ChengZhi (Qincz@lreis.ac.cn)
 * \brief Header file for class GPRO::RasterLayer
 * \version 1.0
 * 
 * \copyright Copyright (c) 2018
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 * changelog:
 *  - 1. 2019-10 - Yujing Wang - Code reformat
 */


#ifndef RasterLayer_H
#define RasterLayer_H

#include "basicTypes.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "metaData.h"
#include "deComposition.h"
#include <fstream>
#include <sstream>

#include "mpi.h"
#include <gdal_priv.h>

using namespace std;

namespace GPRO {
    /**
     * \ingroup gpro
     * \class RasterLayer
     * \brief Basic class to store PaRGO information of layers
     */
    template<class elemType>
    class RasterLayer {
    public:

        RasterLayer();
        RasterLayer( const string layerName );
        RasterLayer( const RasterLayer<elemType> &rhs );
        virtual ~RasterLayer();

        RasterLayer<elemType> &operator=( const RasterLayer<elemType> &rhs );

        const string &name() const;
        void name( const string &layerName );
        unsigned int id() const;
        void id(const int layerID);
        const string title() const;

        void cleanCellSpace();
        void cleanNbrhood();
        void cleanMetaData();

        bool hasCellSpace() const;
        bool hasNbrhood() const;
        bool hasMetaData() const;

        bool newCellSpace();
        bool newCellSpace( const SpaceDims &dims );
        bool newCellSpace( const SpaceDims &dims, const elemType &initVal );
        bool newCellSpace( int nRows, int nCols );
        bool newCellSpace( int nRows, int nCols, const elemType &initVal );

        bool newNbrhood();
        bool newNbrhood( const vector<CellCoord> &vNbrCoords, double weight = 1.0 );
        bool newNbrhood( const vector<CellCoord> &vNbrCoords, const vector<double> &vNbrWeights );
        template<class elemType2>
        bool newNbrhood( const Neighborhood<elemType2> &nbr );

        CellSpace<elemType> *cellSpace();
        const CellSpace<elemType> *cellSpace() const;
        Neighborhood<elemType> *nbrhood();
        const Neighborhood<elemType> *nbrhood() const;
        MetaData *metaData();

        GDALDataType getGDALType();
        MPI_Datatype getMPIType();

        //IO function
        bool readNeighborhood( const char *neighborfile );
		//wyj: ΪʲôĬ��NON_DCMP�ˣ���ʱ�޸ĳ�row dcmp?
        bool readFile( const char *inputfile, DomDcmpType dcmpType = NON_DCMP );
        bool readFile( const char *inputfile, const CoordBR &subWorkBR, DomDcmpType dcmpType = NON_DCMP );
		//ÿ�����̶��������ļ�;ȫ������������̹������������������ã����ߣ������̶��ˣ�Ȼ���ͣ�
        bool readGlobalFile( const char *inputfile, DomDcmpType dcmpType = NON_DCMP );
		//����������Ǵ��ж������ж�����˼Ӧ���������̶�Ȼ�󷢲��������̣�û��Ҫ�����������ʵ�Ƕ�ȫ�����ϲ���subWorkBR��������
		//Todo:��ָ�����з�Χ����

		bool copyLayerInfo( const RasterLayer<elemType> &rhs );

		//void col2row( int &subRows, int &lastSubRows );	//����д���������طֲ�Ϊ����д��
		//IO������+���п鲢��,�����䴮�е�IO
		bool writeFile( const char *outputfile );
        bool updateMetadata( const CoordBR &subWorkBR, DomDcmpType dcmpType = NON_DCMP );
        //bool createFileBlock(const char *outputfile);
        //bool writeFileBlock(const char* outputfile);
    protected:
        bool createFile( const char *outputfile );
        bool rowWriteFile( const char *outputfile );
        bool colWriteFile( const char *outputfile );

    public:
        MetaData *_pMetaData; //wyj 2019-11-12: ����� _pCellSpace��� Metadata��ʲô����

    protected:
        string _strLayerName; /// layer name
        CellSpace<elemType> *_pCellSpace; /// data of raster
        Neighborhood<elemType> *_pNbrhood; /// data of neighborhood
        int _nLayerID; /// layer ID
    };
};

template<class elemType>
inline GPRO::RasterLayer<elemType>::
RasterLayer()
    :
    _strLayerName( "Untitled" ),
    _pCellSpace( 0 ),
    _pNbrhood( 0 ),
    _pMetaData( 0 ),
	_nLayerID(-1) {}

template<class elemType>
inline GPRO::RasterLayer<elemType>::
RasterLayer( const string layerName )
    :_strLayerName( layerName ),
     _pCellSpace( 0 ),
     _pNbrhood( 0 ),
	 _pMetaData( 0 ),
	 _nLayerID(-1) {}

template<class elemType>
inline GPRO::RasterLayer<elemType>::
~RasterLayer() {
    cleanCellSpace();
    cleanNbrhood();
    cleanMetaData();
}

template<class elemType>
inline GPRO::RasterLayer<elemType> &GPRO::RasterLayer<elemType>::
operator=( const RasterLayer<elemType> &rhs ) {
    if ( this == &rhs ) {
        return *this;
    }

    //if( NULL != _strLayerName ){
    if( !_strLayerName.empty() ){
        _strLayerName = rhs._strLayerName + "_copy";
    }

	////Todo:LayerID�ĸ���

    if ( rhs._pCellSpace ) {
        if ( _pCellSpace ) {
			cleanCellSpace();
            *( _pCellSpace ) = *( rhs._pCellSpace );
        } else {
            _pCellSpace = new CellSpace<elemType>( *( rhs._pCellSpace ));
        }
    } else {
        cleanCellSpace();
    }

    if ( rhs._pNbrhood ) {
        if ( _pNbrhood ) {
			cleanNbrhood();
            *( _pNbrhood ) = *( rhs._pNbrhood );
        } else {
            _pNbrhood = new Neighborhood<elemType>( *( rhs._pNbrhood ));
        }
    } else {
        cleanNbrhood();
    }

    return *this;
}

template<class elemType>
inline const string &GPRO::RasterLayer<elemType>::
name() const {
    return _strLayerName;
}

template<class elemType>
inline void GPRO::RasterLayer<elemType>::
name( const string &layerName ) {
    _strLayerName = layerName;
}

template<class elemType>
inline unsigned int GPRO::RasterLayer<elemType>::
id() const {
    //int myID = -1;
    //return myID;
	return _nLayerID;
}

template<class elemType>
inline void GPRO::RasterLayer<elemType>::
id( const int layerID ) {
	if( layerID < 0 ){
		cout << "The layerID should be bigger than zero."<<endl;
	}
	_nLayerID = layerID;
}

template<class elemType>
inline const string GPRO::RasterLayer<elemType>::
title() const {
    ostringstream myTitle;
    //myTitle << _strLayerName << id();
    return myTitle.str();
}

template<class elemType>
void GPRO::RasterLayer<elemType>::
cleanCellSpace() {
    if ( _pCellSpace ) {
        delete _pCellSpace;
        _pCellSpace = 0;
    }
}

template<class elemType>
void GPRO::RasterLayer<elemType>::
cleanNbrhood() {
    if ( _pNbrhood ) {
        delete _pNbrhood;
        _pNbrhood = 0;
    }
}

template<class elemType>
void GPRO::RasterLayer<elemType>::
cleanMetaData() {
    if ( _pMetaData ) {
        delete _pMetaData;
        _pMetaData = 0;
    }
}


template<class elemType>
bool GPRO::RasterLayer<elemType>::
hasCellSpace() const {
    bool hasIt = true;
    if ( !_pCellSpace || _pCellSpace->empty()) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: no CellSpace associated with RasterLayer[" \
			 << title() << "]" << endl;
        hasIt = false;
    }

    return hasIt;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
hasNbrhood() const {
    bool hasIt = true;
    if ( !_pNbrhood ) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: no neighborhood associated with RasterLayer[" \
			 << title() << "]" << endl;
        hasIt = false;
    }

    return hasIt;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
hasMetaData() const {
	bool hasIt = true;
	if ( !_pMetaData ) {
		cerr << __FILE__ << " " << __FUNCTION__ \
			<< " Error: no metaData associated with RasterLayer[" \
			<< title() << "]" << endl;
		hasIt = false;
	}

	return hasIt;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace() {
    cleanCellSpace();

    _pCellSpace = new CellSpace<elemType>();
    if ( !_pCellSpace ) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: unable to new a CellSpace" \
			 << endl;
        return false;
    }

    return true;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace( const SpaceDims &dims ) {
    cleanCellSpace();
    _pCellSpace = new CellSpace<elemType>( dims );
    if ( !_pCellSpace ) {
        cerr << __FILE__ << " " << __FUNCTION__ \
            << " Error: unable to new a CellSpace" \
            << endl;
        return false;
    }

    return true;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace( const SpaceDims &dims, const elemType &initVal ) {
    cleanCellSpace();
    _pCellSpace = new CellSpace<elemType>( dims, initVal );
    if ( !_pCellSpace ) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: unable to new a CellSpace" \
			 << endl;
        return false;
    }

    return true;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace( int nRows, int nCols ) {
    return newCellSpace( SpaceDims( nRows, nCols ));
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace( int nRows, int nCols, const elemType &initVal ) {
    return newCellSpace( SpaceDims( nRows, nCols ), initVal );
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood() {
    cleanNbrhood();
    _pNbrhood = new Neighborhood<elemType>();
    if ( !_pNbrhood ) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: unable to new a Neighborhood on process[" \
			 << title() << "]" << endl;
        return false;
    }
    
	return true;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood( const vector<CellCoord> &vNbrCoords, double weight ) {

    cleanNbrhood();
    _pNbrhood = new Neighborhood<elemType>( vNbrCoords, weight );
    if ( !_pNbrhood ) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: unable to new a Neighborhood on process[" \
			 << title() << "]" << endl;
        return false;
    }

    return true;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood( const vector<CellCoord> &vNbrCoords, const vector<double> &vNbrWeights ) {
	cleanNbrhood();
	_pNbrhood = new Neighborhood<elemType>( vNbrCoords, vNbrWeights );
	if ( !_pNbrhood ) {
		cerr << __FILE__ << " " << __FUNCTION__ \
			<< " Error: unable to new a Neighborhood on process[" \
			<< title() << "]" << endl;
		return false;
	}

	return true;
}

template<class elemType>
template<class elemType2>
bool GPRO::RasterLayer<elemType>::
newNbrhood( const Neighborhood<elemType2> &nbr ) {
	cleanNbrhood();
	_pNbrhood = new Neighborhood<elemType>( nbr );
	if ( !_pNbrhood ) {
		cerr << __FILE__ << " " << __FUNCTION__ \
			<< " Error: unable to new a Neighborhood on process[" \
			<< title() << "]" << endl;
		return false;
	}

	return true;
}

template<class elemType>
inline GPRO::CellSpace<elemType> *GPRO::RasterLayer<elemType>::
cellSpace() {
    return _pCellSpace;
}

template<class elemType>
inline const GPRO::CellSpace<elemType> *GPRO::RasterLayer<elemType>::
cellSpace() const {
    return _pCellSpace;
}

template<class elemType>
inline GPRO::Neighborhood<elemType> *GPRO::RasterLayer<elemType>::
nbrhood() {
    return _pNbrhood;
}

template<class elemType>
inline const GPRO::Neighborhood<elemType> *GPRO::RasterLayer<elemType>::
nbrhood() const {
    return _pNbrhood;
}

template<class elemType>
inline GPRO::MetaData *GPRO::RasterLayer<elemType>::
metaData() {
    return _pMetaData;
}

template<class elemType>
MPI_Datatype GPRO::RasterLayer<elemType>::
getMPIType() {
    elemType style;
    if ( typeid( style ) == typeid( float ) ) {
        return MPI_FLOAT;
    } else if ( typeid( style ) == typeid( double ) ) {
        return MPI_DOUBLE;
        //cout<<"double"<<endl;
    } else if ( typeid( style ) == typeid( int ) ) {
        return MPI_INT;
        //cout<<"int"<<endl;
    } else if ( typeid( style ) == typeid( unsigned int ) ) {
        return MPI_UNSIGNED;
        //cout<<"unsigned int"<<endl;
    } else if ( typeid( style ) == typeid( char ) ) {
        return MPI_CHAR;
        //cout<<"char"<<endl;
    } else {
        return MPI_BYTE;
    }
}

template<class elemType>
GDALDataType GPRO::RasterLayer<elemType>::
getGDALType() {
    elemType style;
    GDALDataType dataType;
    //cout<<typeid(style).name()<<endl;
    if ( typeid( style ) == typeid( float ) ) {
        dataType = GDT_Float32;
        //cout<<"float"<<endl;
    } else if ( typeid( style ) == typeid( double ) ) {
        dataType = GDT_Float64;
        //cout<<"double"<<endl;
    } else if ( typeid( style ) == typeid( int ) ) {
        dataType = GDT_Int32;
        //cout<<"int"<<endl;
    } else if ( typeid( style ) == typeid( unsigned int ) ) {
        dataType = GDT_UInt32;
        //cout<<"unsigned int"<<endl;
    } else if ( typeid( style ) == typeid( char ) ) {
        dataType = GDT_Byte;
        //cout<<"char"<<endl;
    } else {
        dataType = GDT_Unknown;
    }
    return dataType;
}

/* Initialize the Neighborhood object by loading the neighbors stored in a ASCII file */
template<class elemType>
bool GPRO::RasterLayer<elemType>::
readNeighborhood( const char *neighborfile ) {
    newNbrhood();
    fstream nbrFile( neighborfile, ios::in );
    nbrFile >> ( *_pNbrhood );
    nbrFile.close();

    return true;
}

/// ԭserialReadFile
template <class elemType>
bool GPRO::RasterLayer<elemType>::
readGlobalFile(const char* inputfile, DomDcmpType dcmpType)
{
	GDALAllRegister();
	//cout<<inputfile<<endl;
	GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(inputfile, GA_ReadOnly );
	//GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(inputfile, GA_ReadOnly );
	if( poDatasetsrc == NULL /*����Ƿ��������ļ�*/)
	{
		cout<<"[ERROR] data file is not open correct"<<endl;
		exit(1);
	}
	_pMetaData = new MetaData();

	if(_pMetaData == NULL)
	{
		cout<<"[ERROR] MetaData is not allocate correct"<<endl;
		exit(1);
	}

	
	poDatasetsrc->GetGeoTransform(_pMetaData->pTransform);
	_pMetaData->projection = poDatasetsrc->GetProjectionRef();
	GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand( 1 );

	_pMetaData->noData = poBandsrc->GetNoDataValue();
	_pMetaData->row = poBandsrc->GetYSize();
	_pMetaData->column = poBandsrc->GetXSize();
	SpaceDims sdim(_pMetaData->row, _pMetaData->column);
	_pMetaData->_glbDims = sdim;
	_pMetaData->cellSize = _pMetaData->pTransform[1];
	_pMetaData->format = "GTiff";

	MPI_Comm_rank(MPI_COMM_WORLD, &_pMetaData->myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &_pMetaData->processor_number);	//ʵ�ʽ�����������ֻ��������ִ�д˺���

	//���ݿռ���ȫ��,�����ռ�dcmpType����
	_pMetaData->_domDcmpType = dcmpType;	//����������Ǹ�Ϊconst����
	//cout<<"rasterLayer L782 "<<_pMetaData->processor_number<<" dcmpType "<<_pMetaData->_domDcmpType<<endl;

	DeComposition<elemType> deComp(_pMetaData->_glbDims, *_pNbrhood);
	if( dcmpType == ROWWISE_DCMP ){
		//�������ݷ�Χ���л���,��ʼ����metaData._MBR��metaData._localdims��metaData._localworkBR
		deComp.rowDcmp(*_pMetaData, _pMetaData->processor_number);
	}else{
		if( dcmpType == COLWISE_DCMP ){
			deComp.colDcmp(*_pMetaData, _pMetaData->processor_number);
		}else{
			if( dcmpType == BLOCK_DCMP ){			
				cout<<"not support right now."<<endl;//������
			}else{
				if( dcmpType == NON_DCMP ){
					//��Ȿ���̵����ݷ�Χ
					_pMetaData->_localdims = _pMetaData->_glbDims;
					CoordBR _glbWorkBR;
					_pNbrhood->calcWorkBR( _glbWorkBR, _pMetaData->_glbDims );
					_pMetaData->_localworkBR = _glbWorkBR;
					//int glbBegin = _glbWorkBR.nwCorner().iRow();
					//int glbEnd = _glbWorkBR.seCorner().iRow();
					//CellCoord nwCorner(glbBegin + _pNbrhood->minIRow(),
					//	0);
					//CellCoord seCorner(glbEnd + _pNbrhood->maxIRow(),
					//	_pMetaData->_glbDims.nCols() - 1);
					CellCoord nwCorner(0, 0);
					CellCoord seCorner(_pMetaData->_glbDims.nRows()-1, _pMetaData->_glbDims.nCols()-1);
					CoordBR subMBR(nwCorner, seCorner);
					_pMetaData->_MBR = subMBR;
				}else{
					cerr<<"please choose right decomposition type."<<endl;
					return false;
				}
			}
		}
	}
	newCellSpace(_pMetaData->_localdims);

	_pMetaData->dataType = getGDALType();

	poBandsrc->RasterIO(GF_Read, 0, _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pMetaData->dataType, 0, 0);

	if (poDatasetsrc != NULL)
	{
		GDALClose((GDALDatasetH)poDatasetsrc);
		poDatasetsrc = NULL;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	return true;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
copyLayerInfo( const RasterLayer<elemType> &rhs ) {
    _pMetaData = new MetaData();
    if ( _pMetaData == NULL ) {
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit( 1 );
    }

    _pMetaData->cellSize = rhs._pMetaData->cellSize;
    _pMetaData->row = rhs._pMetaData->row;
    _pMetaData->column = rhs._pMetaData->column;
    _pMetaData->format = rhs._pMetaData->format;
    _pMetaData->projection = rhs._pMetaData->projection;
    _pMetaData->noData = rhs._pMetaData->noData;
    _pMetaData->myrank = rhs._pMetaData->myrank;
    _pMetaData->processor_number = rhs._pMetaData->processor_number;
    _pMetaData->_domDcmpType = rhs._pMetaData->_domDcmpType;
    _pMetaData->_glbDims = rhs._pMetaData->_glbDims;
    _pMetaData->_MBR = rhs._pMetaData->_MBR;
    _pMetaData->_localdims = rhs._pMetaData->_localdims;
    _pMetaData->_localworkBR = rhs._pMetaData->_localworkBR;
    _pMetaData->dataType = getGDALType();

    for ( int i = 0; i < 6; i++ ) {
        _pMetaData->pTransform[i] = rhs._pMetaData->pTransform[i];
    }

    newNbrhood( *( rhs.nbrhood()));	//allocate and init
    //����ʼ��Ϊ��ֵ��double���͵�nodata��ֵ��elemtype���ͣ�����Խ��
    //newCellSpace( _pMetaData->_localdims, _pMetaData->noData );	//allocate
    //newCellSpace(_pMetaData->_localdims,0); //allocate
    newCellSpace( _pMetaData->_localdims, rhs._pMetaData->noData );	//allocate

    return true;
}

//����һ��copyLayerInfo�ǲ�ͬ����ͼ����

//Todo:���ڲ��������򣬲����Ƿ�֧��
template<class elemType>
bool GPRO::RasterLayer<elemType>::
readFile( const char *inputfile, DomDcmpType dcmpType ) {
    GDALAllRegister();

    GDALDataset *poDatasetsrc = (GDALDataset *) GDALOpen( inputfile, GA_ReadOnly );
    if ( poDatasetsrc == NULL /*����Ƿ��������ļ�*/) {
        cout << "[ERROR] data file is not open correct" << endl;
        exit( 1 );
    }

    _pMetaData = new MetaData();
    if ( _pMetaData == NULL ) {
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit( 1 );
    }

    poDatasetsrc->GetGeoTransform( _pMetaData->pTransform );
    _pMetaData->projection = poDatasetsrc->GetProjectionRef();
    GDALRasterBand *poBandsrc = poDatasetsrc->GetRasterBand( 1 );

    _pMetaData->noData = poBandsrc->GetNoDataValue();
    _pMetaData->row = poBandsrc->GetYSize();
    _pMetaData->column = poBandsrc->GetXSize();
    _pMetaData->cellSize = _pMetaData->pTransform[1];
    _pMetaData->format = "GTiff";
    _pMetaData->_domDcmpType = dcmpType;

    SpaceDims sdim( _pMetaData->row, _pMetaData->column );
    _pMetaData->_glbDims = sdim;

    MPI_Comm_rank( MPI_COMM_WORLD, &_pMetaData->myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &_pMetaData->processor_number );
    _pMetaData->dataType = getGDALType();

    //�������ݷ�Χ���л��֣����㹤���ռ�ȣ�NON_DCMP��������̶���ȫ����ȫ����orֻ������������д��
    DeComposition<elemType> deComp( _pMetaData->_glbDims, *_pNbrhood );
    if ( ROWWISE_DCMP == _pMetaData->_domDcmpType ) {
        deComp.rowDcmp( *_pMetaData, _pMetaData->processor_number );
    } else {
        if ( COLWISE_DCMP == _pMetaData->_domDcmpType ) {
            deComp.colDcmp( *_pMetaData, _pMetaData->processor_number );
        } else {
            if ( BLOCK_DCMP == _pMetaData->_domDcmpType ) {
                cout << __FILE__ << " " << __FUNCTION__ \
                    << "Error: not support this dcmpType_" << dcmpType \
                    << " right now" << endl;//������
                return false;
            } else {
                if( NON_DCMP == _pMetaData->_domDcmpType ){
                    _pMetaData->_localdims = _pMetaData->_glbDims;
                    CoordBR _glbWorkBR;
                    _pNbrhood->calcWorkBR( _glbWorkBR, _pMetaData->_glbDims );
                    _pMetaData->_localworkBR = _glbWorkBR;
                    CellCoord nwCorner( 0, 0 );
                    CellCoord seCorner( _pMetaData->_glbDims.nRows() - 1, _pMetaData->_glbDims.nCols() - 1 );
                    CoordBR subMBR( nwCorner, seCorner );
                    _pMetaData->_MBR = subMBR;
                }else{
                    cout << __FILE__ << " " << __FUNCTION__ \
                        << "Error: not support this dcmpType_" << dcmpType \
                        << " right now" << endl;
                    return false;
                }
            }
        }
    }
    newCellSpace( _pMetaData->_localdims );

    //���ж�д����ִ������һ��
    //poBandsrc->RasterIO(GF_Read, 0, _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pMetaData->dataType, 0, 0);
    poBandsrc->RasterIO( GF_Read, _pMetaData->_MBR.minICol(), _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
                         _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
                         _pMetaData->dataType, 0, 0 );

    if ( poDatasetsrc != NULL ) {
        //poDatasetsrc->~GDALDataset();
        GDALClose((GDALDatasetH) poDatasetsrc );
        poDatasetsrc = NULL;
    }
    MPI_Barrier( MPI_COMM_WORLD );

    return true;
}

//Todo:���û���˵��CoordBR���ɷ��ʣ���Ϊ����ID
//dcmpType�����Ƿ�����
template<class elemType>
bool GPRO::RasterLayer<elemType>::
readFile( const char *inputfile, const CoordBR &subWorkBR, DomDcmpType dcmpType ) {
    //��֪subWorkBR���Լ��Ĺ����ռ䷶Χ�������Ƶ��������ݷ�Χ��������
    //���ļ�֮ǰ��������֮ǰ���ܲ�������Ϣ;why cleanMetaData�����
    if ( _pMetaData ) {
        cleanMetaData();
    }    //Ҳ���Բ�clean��ֱ�Ӹ���
    //cleanCellSpace();	//����û��Ҫclean����Ϊ����new��ʱ����ȵ���clean��������д�ͻ��ظ�clean��û��Ҫ

    GDALAllRegister();

    GDALDataset *poDatasetsrc = (GDALDataset *) GDALOpen( inputfile, GA_ReadOnly );
    if ( poDatasetsrc == NULL /*����Ƿ��������ļ�*/) {
        cout << "[ERROR] data file is not open correct" << endl;
        exit( 1 );
    }

    _pMetaData = new MetaData();    //Ҳ���Բ���Ҫ��new��ֱ�Ӹ���ԭ��Ϣ
    if ( _pMetaData == NULL ) {
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit( 1 );
    }

    poDatasetsrc->GetGeoTransform( _pMetaData->pTransform );
    _pMetaData->projection = poDatasetsrc->GetProjectionRef();
    GDALRasterBand *poBandsrc = poDatasetsrc->GetRasterBand( 1 );

    _pMetaData->noData = poBandsrc->GetNoDataValue();
    _pMetaData->row = poBandsrc->GetYSize();
    _pMetaData->column = poBandsrc->GetXSize();
	_pMetaData->cellSize = _pMetaData->pTransform[1];
	_pMetaData->format = "GTiff";
	_pMetaData->_domDcmpType = dcmpType;

	SpaceDims sdim( _pMetaData->row, _pMetaData->column );
    _pMetaData->_glbDims = sdim;

    MPI_Comm_rank( MPI_COMM_WORLD, &_pMetaData->myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &_pMetaData->processor_number );

    int glbBegin = subWorkBR.nwCorner().iRow();
    int glbEnd = subWorkBR.seCorner().iRow();
    CellCoord nwCorner( glbBegin + _pNbrhood->minIRow(), 0 );
    CellCoord seCorner( glbEnd + _pNbrhood->maxIRow(), _pMetaData->_glbDims.nCols() - 1 );
    CoordBR subMBR( nwCorner, seCorner );
    _pMetaData->_MBR = subMBR;
    SpaceDims dims( subMBR.nRows(), subMBR.nCols());
    _pMetaData->_localdims = dims;
    //subWorkBR��ȫ�����кţ�����ֲ����кű���_localworkBR
    CoordBR workBR;
    if ( !_pNbrhood->calcWorkBR( workBR, dims )) {
        return false;
    }
    _pMetaData->_localworkBR = workBR;

    newCellSpace( _pMetaData->_localdims );

    _pMetaData->dataType = getGDALType();

	poBandsrc->RasterIO( GF_Read, 0, _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
						 _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
                         _pMetaData->dataType, 0, 0 );

    if ( poDatasetsrc != NULL ) {
        GDALClose((GDALDatasetH) poDatasetsrc );
        poDatasetsrc = NULL;
    }
    MPI_Barrier( MPI_COMM_WORLD );	//�Ƿ���Ҫ

    return true;
}

//template<class elemType>
//void GPRO::RasterLayer<elemType>::
//col2row( int &subRows, int &lastSubRows ){
//	int myRank, process_nums;
//	MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
//	MPI_Comm_size( MPI_COMM_WORLD, &process_nums );	//��������ȡ��������ʽ�͵���pMetaData�ķ�ʽ��������ʲô����������ֽ��̽�����ʱ��
//	MPI_Datatype datatype = getMPIType();
//	
//	//��һ�׶Σ��طֲ�
//	//������nRows(),�ҳ������ػ��ֵĸ���������Ҫ�е�������
//	if ( _pMetaData->_glbDims.nRows() % process_nums == 0 ) {
//		subRows = _pMetaData->_glbDims.nRows() / process_nums;
//		lastSubRows = subRows;
//	} else {
//		subRows = _pMetaData->_glbDims.nRows() / process_nums + 1;
//		lastSubRows = _pMetaData->_glbDims.nRows() - ( _pMetaData->_glbDims.nRows() / process_nums + 1 ) * ( process_nums - 1 );
//	}
//	int glbWorkCols = _pMetaData->_glbDims.nCols() + _pNbrhood->minICol() - _pNbrhood->maxICol();
//
//	//�����¿ռ��������;�Ƚ��Լ����룬������ÿ����һ�����ڴ�һ��
//	int cellCount = 0;
//	if ( myRank < process_nums - 1 ) {
//		cellCount = subRows * _pMetaData->_glbDims.nCols();
//	} else {
//		cellCount = lastSubRows * _pMetaData->_glbDims.nCols();
//	}
//	elemType *pRowMatrix = new elemType[cellCount];
//	//cout<<myRank<<" "<<subRows<<" "<<lastSubRows<<" "<<glbWorkCols<<endl;
//	for ( int rRow = 0; rRow < subRows; ++rRow ) {
//		if ( myRank == 0 ) {
//			for ( int rCol = 0; rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol ) {
//				pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
//			}
//		} else {
//			if ( myRank == process_nums - 1 ) {
//				for ( int rCol = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
//					rCol <= _pMetaData->_MBR.maxICol() + _pNbrhood->maxICol(); ++rCol ) {
//						pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] =
//							_pCellSpace->_matrix[( rRow + myRank * subRows ) * _pMetaData->_localdims.nCols() + rCol - _pMetaData->_MBR.minICol()];
//				}
//			} else {
//				for ( int rCol = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
//					rCol <= _pMetaData->_MBR.maxICol(); ++rCol ) {
//						pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] =
//							_pCellSpace->_matrix[( rRow + myRank * subRows ) * _pMetaData->_localdims.nCols() + rCol - _pMetaData->_MBR.minICol()];
//				}
//			}
//		}
//	}
//	//cout<<myRank<<" "<<_pCellSpace->_matrix[(0+subRows*myRank)*_pMetaData->_localdims.nCols()+_pMetaData->_localworkBR.minICol()]<<endl;
//	//cout<<myRank<<" "<<_pCellSpace->_matrix[(1+subRows*myRank)*_pMetaData->_localdims.nCols()+_pMetaData->_localworkBR.minICol()]<<endl;
//
//	MPI_Status status;
//	//ȷ��rankID�ĵ�myRank���СΪ����*myrank��������������ô�����ʱ�ռ䣬���ͣ����գ����£��ͷ�
//	//Q�ѵ���С���ǹ̶��ģ�Ϊʲô��Ҫ��ͨ����ʱ�ռ�Ĵ�С
//	for ( int objRank = 0; objRank < process_nums; ++objRank ) {
//		if ( objRank == myRank ) {
//			continue;
//		}
//		//ͨ��objRank��workBR.ncols()��MBR.minCol()��myRank
//		int sendnCols, sendColPos;
//		if ( myRank == 0 ) {
//			sendnCols = _pMetaData->_localworkBR.nCols() + 1;
//			sendColPos = 0;
//		} else {
//			if ( myRank == process_nums - 1 ) {
//				sendnCols = _pMetaData->_localworkBR.nCols() + 1;
//				sendColPos = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
//			} else {
//				sendnCols = _pMetaData->_localworkBR.nCols();
//				sendColPos = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
//			}
//		}
//		int recvnCols, recvColPos;
//		MPI_Send( &sendnCols, 1, MPI_INT, objRank, 1, MPI_COMM_WORLD );
//		MPI_Recv( &recvnCols, 1, MPI_INT, objRank, 1, MPI_COMM_WORLD, &status );
//		MPI_Send( &sendColPos, 1, MPI_INT, objRank, 2, MPI_COMM_WORLD );
//		MPI_Recv( &recvColPos, 1, MPI_INT, objRank, 2, MPI_COMM_WORLD, &status );
//		int sendMatrixSize, recvMatrixSize;
//		if ( myRank < process_nums - 1 ) {
//			recvMatrixSize = recvnCols * subRows;
//		} else {
//			recvMatrixSize = recvnCols * lastSubRows;
//		}
//		if ( objRank < process_nums - 1 ) {
//			sendMatrixSize = sendnCols * subRows;
//		} else {
//			sendMatrixSize = sendnCols * lastSubRows;
//		}
//		elemType *sendMatrix = new elemType[sendMatrixSize];
//		elemType *recvMatrix = new elemType[recvMatrixSize];
//		//��myRankҪ��rankID���յĴ�СΪ���ҵ�localRows*rankID�Ĺ������������������Ҫ��ͨ������ȡ��
//		//Ҫ���͵ĺͽ��յ����ݶ����ȴ�����ʱһά���飬ͨ�ź���ӳ��
//		cellCount = 0;
//		for ( int rRow = objRank * subRows; rRow < ( objRank + 1 ) * subRows; ++rRow ) {
//			if ( rRow >= _pMetaData->_glbDims.nRows() ) {
//				break;
//			}
//			if ( myRank == 0 ) {
//				for ( int rCol = 0; rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol ) {
//					sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
//				}
//			} else {
//				if ( myRank == process_nums - 1 ) {
//					for ( int rCol = _pMetaData->_localworkBR.minICol(); rCol <= _pMetaData->_localworkBR.maxICol() + _pNbrhood->maxICol(); ++rCol ) {
//						sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
//					}
//				} else {
//					for ( int rCol = _pMetaData->_localworkBR.minICol(); rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol ) {
//						sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
//					}
//				}
//			}
//		}
//
//		int styleSize;
//		MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize );
//		int bufsize = MPI_BSEND_OVERHEAD + styleSize * ( sendMatrixSize + recvMatrixSize );
//		void *buf = malloc( bufsize );
//		MPI_Buffer_attach( buf, bufsize );
//
//		MPI_Bsend( sendMatrix, sendMatrixSize, datatype, objRank, 1, MPI_COMM_WORLD );
//		MPI_Recv( recvMatrix, recvMatrixSize, datatype, objRank, 1, MPI_COMM_WORLD, &status );
//		cellCount = 0;
//		//�����յ�recvMatrix�������
//		for ( int rRow = 0; rRow < subRows; ++rRow ) {
//			if ( myRank == process_nums - 1 && rRow >= lastSubRows ) {
//				break;
//			}
//			for ( int rCol = recvColPos; rCol < recvColPos + recvnCols; ++rCol ) {
//				pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] = recvMatrix[cellCount++];
//			}
//		}
//		//if( myRank == 1 )
//		//cout<<recvColPos<<" "<<recvCols<<endl;
//		MPI_Buffer_detach( &buf, &bufsize );
//		MPI_Buffer_attach( buf, bufsize );
//	}
//}

// �����̴���դ���ļ�
template<class elemType>
bool GPRO::RasterLayer<elemType>::
createFile( const char *outputfile ) {
    GDALAllRegister();

    if ( _pMetaData->myrank == 0 ) {
        GDALDriver *poDriver = NULL;
        poDriver = GetGDALDriverManager()->GetDriverByName( _pMetaData->format.c_str());
        if ( poDriver == NULL ) {
            cout << "poDriver is NULL." << endl;
            return false;
        }
        char **papszMetadata = NULL;
        papszMetadata = CSLSetNameValue( papszMetadata, "BLOCKXSIZE", "256" );
        papszMetadata = CSLSetNameValue( papszMetadata, "BLOCKYSIZE", "1" );

        if ( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) );
        //cout<< "This format driver supports Create() method."<<endl;
        if ( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE ) );
        //cout<< "This format driver supports CreateCopy() method."<<endl;

        GDALDataset *poDataset = poDriver->Create( outputfile, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                                                   1, _pMetaData->dataType, papszMetadata );
        if ( poDataset == NULL ) {
            cout << "poDatasetdest is NULL" << endl;
            return false;
        }
        poDataset->SetGeoTransform( _pMetaData->pTransform );
        poDataset->SetProjection( _pMetaData->projection.c_str() );

        if ( poDataset != NULL ) {
            GDALClose((GDALDatasetH) poDataset );
            poDataset = NULL;
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );//wyj 2019-11-12:��д�����������ļ���ʱ�򣬺����������̵��ã�Barrier�������

    return true;
}
 

//ֻд���˹����ռ�
//�����Բ��Գ�������-1/0/1/2������Χ�ģ�
template<class elemType>
bool GPRO::RasterLayer<elemType>::
writeFile( const char *outputfile ){
    //������Ĭ��ʵ�ְ���д��������������������ͣ����ö�Ӧ��д����
    if( COLWISE_DCMP == _pMetaData->_domDcmpType ){
        return colWriteFile( outputfile );
    }else if( BLOCK_DCMP == _pMetaData->_domDcmpType ){
        cout << __FILE__ << " " << __FUNCTION__ \
            << "Error: not support this dcmpType_" << _pMetaData->_domDcmpType \
            << " right now" << endl;	//�����
        return false;
    }else if( NON_DCMP == _pMetaData->_domDcmpType ){
        //����ɣ�ֻ����������д?or��֧��
    }else{
        //cout << __FILE__ << " " << __FUNCTION__ \
            << "Warning: not support this dcmpType_" << _pMetaData->_domDcmpType \
            << " right now, using row-wise decomposition as default." << endl;
    }
    return rowWriteFile(outputfile);
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
rowWriteFile( const char *outputfile ) {
   GDALAllRegister();

   if ( !createFile( outputfile ) ) {
		cout << __FILE__ << " " << __FUNCTION__ \
			<< " Error: create file failed!" << endl;
       MPI_Finalize();
   }
   ////����ǲ��ֽ���������д�����أ������޵ȴ��� // wyj:��
   GDALDataset *poDataset = NULL;
   poDataset = (GDALDataset *) GDALOpen( outputfile, GA_Update );
   if ( poDataset == NULL /*����Ƿ��������ļ�*/) {
       cout << "data file is not open correct" << endl;
       exit( 1 );
   }

   CellSpace<elemType> &computL = *cellSpace();
   //for (int i=0;i<RasterLayer<elemType>::_pCellSpace->nRows();i++)
   //{
	  // for (int j=0;j<RasterLayer<elemType>::_pCellSpace->nCols();j++)
	  // {
		 //  if(computL[i][j]>10){
			//   cout<<i<<" "<<j<<" "<<computL[i][j]<<" - rank"<<_pMetaData->myrank;
			//   cout<<endl;
		 //  }
	  // }
   //}
   GDALRasterBand *poBanddest = poDataset->GetRasterBand( 1 );
   if ( poBanddest == NULL ) {
       cout << "poBanddest is NULL" << endl;
       exit( 1 );
   }
   if ( _pMetaData->myrank == 0 ) {
       poBanddest->SetNoDataValue( _pMetaData->noData );	//Ϊɶ�ǽ��������̣�����ÿ�����̶���
   }

   if ( _pMetaData->processor_number == 1 ) {
       poBanddest->RasterIO( GF_Write, 0, 0, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                             _pCellSpace->_matrix, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                             _pMetaData->dataType, 0, 0 );
   } else {
       if ( _pMetaData->myrank == 0 ) {
           poBanddest->RasterIO( GF_Write, 0, 0, _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow() + 1,
                                 _pCellSpace->_matrix, _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow() + 1,
                                 _pMetaData->dataType, 0, 0 );
       } else if ( _pMetaData->myrank == ( _pMetaData->processor_number - 1 ) ) {
           poBanddest->RasterIO( GF_Write, 0, _pMetaData->_MBR.minIRow() - _pNbrhood->minIRow(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow() + 1,
                                 _pCellSpace->_matrix - _pNbrhood->minIRow() * _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow() + 1,
                                 _pMetaData->dataType, 0, 0 );
       } else {
           poBanddest->RasterIO( GF_Write, 0, _pMetaData->_MBR.minIRow() - _pNbrhood->minIRow(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow() + _pNbrhood->minIRow() + 1,
                                 _pCellSpace->_matrix - _pNbrhood->minIRow() * _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow() + _pNbrhood->minIRow() + 1,
                                 _pMetaData->dataType, 0, 0 );
       }
   }

   MPI_Barrier( MPI_COMM_WORLD );	//�Ƿ���Ҫ
   if ( poDataset != NULL ) {
       GDALClose((GDALDatasetH) poDataset );
       poDataset = NULL;
   }
	return true;
}

//ֻд���˹����ռ�
//parallel IO��refer to Qin13_TransactionsInGIS
template<class elemType>
bool GPRO::RasterLayer<elemType>::
colWriteFile( const char *outputfile ) {
    int myRank, process_nums;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &process_nums );
    MPI_Datatype datatype = getMPIType();

    GDALAllRegister();
    if ( !createFile( outputfile )) {
        cout << "create file is not correct!" << endl;
        MPI_Finalize();
    }

    GDALDataset *poDataset = NULL;
    poDataset = (GDALDataset *) GDALOpen( outputfile, GA_Update );
    if ( poDataset == NULL /*����Ƿ��������ļ�*/) {
        cout << "data file is not open correct" << endl;
        exit( 1 );
    }
    GDALRasterBand *poBanddest = poDataset->GetRasterBand( 1 );
    if ( poBanddest == NULL ) {
        cout << "poBanddest is NULL" << endl;
        exit( 1 );
    }
    if ( _pMetaData->myrank == 0 ) {
        poBanddest->SetNoDataValue( _pMetaData->noData );
    }

    if ( _pMetaData->processor_number == 1 ) {
        poBanddest->RasterIO( GF_Write, 0, 0, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                              _pCellSpace->_matrix, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                              _pMetaData->dataType, 0, 0 );
    } else {
        //��һ�׶Σ��طֲ�
        //������nRows(),�ҳ������ػ��ֵĸ���������Ҫ�е�������
        int subRows = 0, lastSubRows = 0;
        if ( _pMetaData->_glbDims.nRows() % process_nums == 0 ) {
            subRows = _pMetaData->_glbDims.nRows() / process_nums;
            lastSubRows = subRows;
        } else {
            subRows = _pMetaData->_glbDims.nRows() / process_nums + 1;
            lastSubRows = _pMetaData->_glbDims.nRows() - ( _pMetaData->_glbDims.nRows() / process_nums + 1 ) * ( process_nums - 1 );
        }
        int glbWorkCols = _pMetaData->_glbDims.nCols() + _pNbrhood->minICol() - _pNbrhood->maxICol();
        //�����¿ռ��������;�Ƚ��Լ����룬������ÿ����һ�����ڴ�һ��
        int cellCount = 0;
        if ( myRank < process_nums - 1 ) {
            cellCount = subRows * _pMetaData->_glbDims.nCols();
        } else {
            cellCount = lastSubRows * _pMetaData->_glbDims.nCols();
        }

        elemType *pRowMatrix = new elemType[cellCount];
        //cout<<myRank<<" "<<subRows<<" "<<lastSubRows<<" "<<glbWorkCols<<endl;
        for ( int rRow = 0; rRow < subRows; ++rRow ) {
            if ( myRank == 0 ) {
                for ( int rCol = 0; rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol ) {
                    pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                }
            } else {
                if ( myRank == process_nums - 1 ) {
                    for ( int rCol = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                          rCol <= _pMetaData->_MBR.maxICol() + _pNbrhood->maxICol(); ++rCol ) {
                              pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] =
                                  _pCellSpace->_matrix[( rRow + myRank * subRows ) * _pMetaData->_localdims.nCols() + rCol - _pMetaData->_MBR.minICol()];
                    }
                } else {
                    for ( int rCol = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                          rCol <= _pMetaData->_MBR.maxICol(); ++rCol ) {
                              pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] =
                                  _pCellSpace->_matrix[( rRow + myRank * subRows ) * _pMetaData->_localdims.nCols() + rCol - _pMetaData->_MBR.minICol()];
                    }
                }
            }
        }
        //cout<<myRank<<" "<<_pCellSpace->_matrix[(0+subRows*myRank)*_pMetaData->_localdims.nCols()+_pMetaData->_localworkBR.minICol()]<<endl;
        //cout<<myRank<<" "<<_pCellSpace->_matrix[(1+subRows*myRank)*_pMetaData->_localdims.nCols()+_pMetaData->_localworkBR.minICol()]<<endl;
        if( NULL == pRowMatrix ){
            cout<<"c "<<myRank<<endl;
        }
        MPI_Status status;
        //ȷ��rankID�ĵ�myRank���СΪ����*myrank��������������ô�����ʱ�ռ䣬���ͣ����գ����£��ͷ�
        //Q�ѵ���С���ǹ̶��ģ�Ϊʲô��Ҫ��ͨ����ʱ�ռ�Ĵ�С
        for ( int objRank = 0; objRank < process_nums; ++objRank ) {
            if ( objRank == myRank ) {
                continue;
            }
            //ͨ��objRank��workBR.ncols()��MBR.minCol()��myRank
            int sendnCols, sendColPos;
            if ( myRank == 0 ) {
                sendnCols = _pMetaData->_localworkBR.nCols() + 1;
                sendColPos = 0;
            } else {
                if ( myRank == process_nums - 1 ) {
                    sendnCols = _pMetaData->_localworkBR.nCols() + 1;
                    sendColPos = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                } else {
                    sendnCols = _pMetaData->_localworkBR.nCols();
                    sendColPos = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                }
            }
            int recvnCols, recvColPos;
            MPI_Send( &sendnCols, 1, MPI_INT, objRank, 1, MPI_COMM_WORLD );
            MPI_Recv( &recvnCols, 1, MPI_INT, objRank, 1, MPI_COMM_WORLD, &status );
            MPI_Send( &sendColPos, 1, MPI_INT, objRank, 2, MPI_COMM_WORLD );
            MPI_Recv( &recvColPos, 1, MPI_INT, objRank, 2, MPI_COMM_WORLD, &status );
            int sendMatrixSize, recvMatrixSize;
            if ( myRank < process_nums - 1 ) {
                recvMatrixSize = recvnCols * subRows;
            } else {
                recvMatrixSize = recvnCols * lastSubRows;
            }
            if ( objRank < process_nums - 1 ) {
                sendMatrixSize = sendnCols * subRows;
            } else {
                sendMatrixSize = sendnCols * lastSubRows;
            }
            elemType *sendMatrix = new elemType[sendMatrixSize];
            elemType *recvMatrix = new elemType[recvMatrixSize];
            //��myRankҪ��rankID���յĴ�СΪ���ҵ�localRows*rankID�Ĺ������������������Ҫ��ͨ������ȡ��
            //Ҫ���͵ĺͽ��յ����ݶ����ȴ�����ʱһά���飬ͨ�ź���ӳ��
            cellCount = 0;
            for ( int rRow = objRank * subRows; rRow < ( objRank + 1 ) * subRows; ++rRow ) {
                if ( rRow >= _pMetaData->_glbDims.nRows() ) {
                    break;
                }
                if ( myRank == 0 ) {
                    for ( int rCol = 0; rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol ) {
                        sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                    }
                } else {
                    if ( myRank == process_nums - 1 ) {
                        for ( int rCol = _pMetaData->_localworkBR.minICol(); rCol <= _pMetaData->_localworkBR.maxICol() + _pNbrhood->maxICol(); ++rCol ) {
                            sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                        }
                    } else {
                        for ( int rCol = _pMetaData->_localworkBR.minICol(); rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol ) {
                            sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                        }
                    }
                }
            }

            int styleSize;
            MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize );
            int bufsize = MPI_BSEND_OVERHEAD + styleSize * ( sendMatrixSize + recvMatrixSize );
            void *buf = malloc( bufsize );
            MPI_Buffer_attach( buf, bufsize );

            MPI_Bsend( sendMatrix, sendMatrixSize, datatype, objRank, 1, MPI_COMM_WORLD );
            MPI_Recv( recvMatrix, recvMatrixSize, datatype, objRank, 1, MPI_COMM_WORLD, &status );
            cellCount = 0;
            //�����յ�recvMatrix�������
            for ( int rRow = 0; rRow < subRows; ++rRow ) {
                if ( myRank == process_nums - 1 && rRow >= lastSubRows ) {
                    break;
                }
                for ( int rCol = recvColPos; rCol < recvColPos + recvnCols; ++rCol ) {
                    pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] = recvMatrix[cellCount++];
                }
            }
            //if( myRank == 1 )
            //cout<<recvColPos<<" "<<recvCols<<endl;
            MPI_Buffer_detach( &buf, &bufsize );
            MPI_Buffer_attach( buf, bufsize );	//����Ҫ��
        }

        if ( myRank < process_nums - 1 ) {
            poBanddest->RasterIO( GF_Write, 0, myRank * subRows, _pMetaData->_glbDims.nCols(), subRows,
                                  pRowMatrix, _pMetaData->_glbDims.nCols(), subRows,
                                  _pMetaData->dataType, 0, 0 );
        } else {
            cout<<"a "<<myRank<<" "<<pRowMatrix<<endl;
            poBanddest->RasterIO( GF_Write, 0, myRank * subRows, _pMetaData->_glbDims.nCols(), lastSubRows,
                                  pRowMatrix, _pMetaData->_glbDims.nCols(), lastSubRows,
                                  _pMetaData->dataType, 0, 0 );
        }
    }
    cout<<"b "<<myRank<<endl;

    MPI_Barrier( MPI_COMM_WORLD );	//ͬ�ϣ�ֻ�в��ֽ��̵��ã������ʲô���
    if ( poDataset != NULL ) {
        GDALClose((GDALDatasetH) poDataset );
        poDataset = NULL;
    }

    return true;
}

template<class elemType>
bool GPRO::RasterLayer<elemType>::
updateMetadata( const CoordBR &subWorkBR, DomDcmpType dcmpType ) {
    //����subWorkBR����Ԫ���ݣ�����֮�󣬻���ô������»�ȡԪ���ݣ�
    if ( _pMetaData == NULL ) {
        cerr << "origional metadata does not exist." << endl;
    }
    int glbBegin = subWorkBR.nwCorner().iRow();
    int glbEnd = subWorkBR.seCorner().iRow();
    CellCoord nwCorner( glbBegin + _pNbrhood->minIRow(), 0 );
    CellCoord seCorner( glbEnd + _pNbrhood->maxIRow(), _pMetaData->_glbDims.nCols() - 1 );
    CoordBR subMBR( nwCorner, seCorner );
    _pMetaData->_MBR = subMBR;
    SpaceDims dims( subMBR.nRows(), subMBR.nCols());
    _pMetaData->_localdims = dims;
    CoordBR workBR;
    if ( !_pNbrhood->calcWorkBR( workBR, dims )) {
        return false;
    }
    _pMetaData->_localworkBR = workBR;

    cleanCellSpace();    //�ͷ�ԭ���Ŀռ�
    newCellSpace( _pMetaData->_localdims );

    return true;
}

#endif
