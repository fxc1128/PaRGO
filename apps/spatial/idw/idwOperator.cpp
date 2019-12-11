#include <ogrsf_frmts.h>
#include "utility.h"
//#include <gdal_priv.h>

#include "idwOperator.h"

int IDWOperator::getBlockColIndex(double x) {
    return (x - _glb_extent.minX) / _cellSize / _blockGrain;
}

int IDWOperator::getBlockRowIndex(double y) {
    return (_glb_extent.maxY - y) / _cellSize / _blockGrain;
}

IDWOperator::~IDWOperator() {
    //delete _pSampleBlocks
}

int IDWOperator::readSampleNums(const char* filename, char** pSpatialRefWkt) {
    //��ȡʸ�������Ԫ���ݣ���ȡ��Χ
#if GDAL_VERSION_MAJOR >= 2
	GDALAllRegister();
	GDALDataset* poDatasetsrc = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);
#else
    OGRRegisterAll();
    OGRDataSource* poDatasetsrc = OGRSFDriverRegistrar::Open(filename, FALSE);
#endif
    if (poDatasetsrc == NULL) {
        printf("[ERROR] Open failed.\n");
        exit(1);
    }

    string file = filename;
    string f2 = file.substr(0, file.length() - 4);
    int pos = f2.find_last_of(SEP); //ע�⣬linux��'/',windows��'\\'
    string f3 = f2.substr(pos + 1);

    OGRLayer* poLayer = poDatasetsrc->GetLayerByName(f3.c_str()); //f3���ļ�����������׺
    OGRSpatialReference* sref = poLayer->GetSpatialRef();
    sref->exportToWkt(pSpatialRefWkt);

    OGRFeature* poFeature;

    poLayer->ResetReading();
    while ((poFeature = poLayer->GetNextFeature()) != NULL) {
        _sample_nums++;
        OGRFeature::DestroyFeature(poFeature);
    }
    //_sample_nums = poLayer->GetFeatureCount();	//Ϊʲô��ֱ�����������
#if GDAL_VERSION_MAJOR >= 2
	GDALClose(poDatasetsrc);
#else
    OGRDataSource::DestroyDataSource(poDatasetsrc);
#endif
    return _sample_nums;
}

bool IDWOperator::readSamples(const char* filename, int fieldIdx, char** pSpatialRefWkt, double** Sample_Array) {

    //��λ����Ϣ��������Ϣ���������Sample_Array��
#if GDAL_VERSION_MAJOR >= 2
	GDALAllRegister();
	GDALDataset* poDatasetsrc = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);
#else
    OGRRegisterAll();
    OGRDataSource* poDatasetsrc = OGRSFDriverRegistrar::Open(filename, FALSE);
#endif
    if (poDatasetsrc == NULL) {
        printf("[ERROR] Open failed.\n");
        exit(1);
    }

    string file = filename;
    string f2 = file.substr(0, file.length() - 4);
    int pos = f2.find_last_of(SEP);
    string f3 = f2.substr(pos + 1);
    OGRLayer* poLayer = poDatasetsrc->GetLayerByName(f3.c_str());
    poLayer->ResetReading();

    int idx = 0;
    double x = 0.0;
    double y = 0.0;
    OGRFeature* poFeature;
    while ((poFeature = poLayer->GetNextFeature()) != NULL) {
        Sample_Array[idx][2] = poFeature->GetFieldAsDouble(fieldIdx); //��ȡ����ֵ
        //cout<<"value:"<<showpoint<<Sample_Array[idx][2]<<endl;	//����û���⣬��С��
        OGRGeometry* poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        if (poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint) {
            //����Χ��һ��Ϊʲô����poLayer->GetExtent��
            OGRPoint* poPoint = (OGRPoint *)poGeometry;
            x = poPoint->getX();
            y = poPoint->getY();
            //�洢λ����Ϣ
            Sample_Array[idx][0] = x;
            Sample_Array[idx][1] = y;
            if (idx == 0) {
                _glb_extent.minX = x;
                _glb_extent.maxX = x;
                _glb_extent.minY = y;
                _glb_extent.maxY = y;
            }
            else {
                if (x > _glb_extent.maxX)
                    _glb_extent.maxX = x;
                if (x < _glb_extent.minX)
                    _glb_extent.minX = x;
                if (y > _glb_extent.maxY)
                    _glb_extent.maxY = y;
                if (y < _glb_extent.minY)
                    _glb_extent.minY = y;
            }
        }
        else {
            printf("[ERROR] No point geometry\n");
            return 1;
        }
        OGRFeature::DestroyFeature(poFeature);
        idx++;
    }

    _glb_extent.minX = _glb_extent.minX - _cellSize / 2;
    int totalcol = (_glb_extent.maxX - _glb_extent.minX) / _cellSize;
    if ((_glb_extent.maxX - _glb_extent.minX) != totalcol * _cellSize) {
        totalcol++;
    }
    _glb_extent.maxX = _glb_extent.minX + _cellSize * totalcol;

    _glb_extent.minY = _glb_extent.minY - _cellSize / 2;
    int totalrow = (_glb_extent.maxY - _glb_extent.minY) / _cellSize;
    if ((_glb_extent.maxY - _glb_extent.minY) != totalrow * _cellSize) {
        totalrow++;
    }
    _glb_extent.maxY = _glb_extent.minY + _cellSize * totalrow;
    //cout<<idw.extent_All.minX<<"	"<<idw.extent_All.maxX<<endl;
    _nRows = totalrow;
    _nCols = totalcol;

#if GDAL_VERSION_MAJOR >= 2
	GDALClose(poDatasetsrc);
#else
    OGRDataSource::DestroyDataSource(poDatasetsrc);
#endif
    return true;
}

void IDWOperator::creatSampleBlocks(double** pSamples) {
    //˼·����idwLayer��������ֿ飬1D�洢,��ȡ�ÿ�ķ�Χ,��ʱ����Ҫ��Χ
    //ע���������귶Χ��դ����ƫ�ƣ����������кţ�
    //ÿ����֪�����xy�����Ƴ����ڿ�;(x-minX)/cellSize/blockGrain�Ϳ�������к�,ͬ�����кţ�
    _blockRows = _nRows / _blockGrain;
    _blockCols = _nCols / _blockGrain;
    _blockRows += (_nRows % _blockGrain) ? 1 : 0;
    _blockCols += (_nCols % _blockGrain) ? 1 : 0;
    _pSampleBlocks = new Sample_block[_blockRows * _blockCols];
    for (int i = 0; i < _sample_nums; ++i) {
        double x = pSamples[i][0];
        double y = pSamples[i][1];
        double z = pSamples[i][2];
        int iCol = getBlockColIndex(x);
        int iRow = getBlockRowIndex(y);
        Sample_Point tmpPoint = {x, y, z};
        _pSampleBlocks[iRow * _blockCols + iCol].sample_Points.push_back(tmpPoint); //pushback�ǿ���ֵ��
    }
}
 
void IDWOperator::idwLayer(RasterLayer<double>& layerD, char** pSpatialRefWkt) {
    //����_pIDWLayer/layerD�Ļ���Ԫ����;������extent��_cellSize��Ϣ������դ��ͼ��
    //MetaData **pMetaData = &(layerD._pMetaData);	//���Կ�����ָ���ָ���д
    layerD._pMetaData = new MetaData();
    if (layerD._pMetaData == NULL) {
        //do something
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit(1);
    }

    layerD._pMetaData->noData = _noData;
    //cout<<NODATA_DEFINE<<" = "<<layerD._pMetaData->noData<<endl;
    layerD._pMetaData->row = _nRows;
    layerD._pMetaData->column = _nCols;
    SpaceDims sdim(layerD._pMetaData->row, layerD._pMetaData->column);
    layerD._pMetaData->_glbDims = sdim;
    // _pMetaData->pTransform
    layerD._pMetaData->cellSize = _cellSize;
    layerD._pMetaData->format = "GTiff";
    layerD._pMetaData->_domDcmpType = ROWWISE_DCMP; //�Ƿ���Ҫ������ָ��
    MPI_Comm_rank(MPI_COMM_WORLD, &layerD._pMetaData->myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &layerD._pMetaData->processor_number);
    _myRank = layerD._pMetaData->myrank;

    DeComposition<double> deComp(layerD._pMetaData->_glbDims, *(layerD.nbrhood()));
    deComp.rowDcmp(*(layerD._pMetaData), layerD._pMetaData->processor_number); //�������ݷ�Χ���л���,���÷�ʽ���ظ�*(layerD._pMetaData)
    layerD.newCellSpace(layerD._pMetaData->_localdims); //ÿ�����̶���ȫ���������ݣ���ֻ���Լ���workBR
    layerD._pMetaData->_glbDims.nRows(_nRows);
    layerD._pMetaData->_glbDims.nCols(_nCols);

    layerD._pMetaData->dataType = layerD.getGDALType();
    //pSpatialRefWkt Ŀǰָ��main�����е�char* pSpatialRefWkt�ĵ�ַ
    layerD._pMetaData->projection = *pSpatialRefWkt; //char* to string,ֱ�Ӹ�ֵ���ɣ�string to char*,����c_str()
    //cout<<"layerD._pMetaData->projection "<<layerD._pMetaData->projection<<endl;
    layerD._pMetaData->pTransform[0] = _glb_extent.minX;
    layerD._pMetaData->pTransform[1] = _cellSize;
    layerD._pMetaData->pTransform[2] = 0;
    layerD._pMetaData->pTransform[3] = _glb_extent.maxY;
    layerD._pMetaData->pTransform[4] = 0;
    layerD._pMetaData->pTransform[5] = -_cellSize;

    //�����ӿռ����ݷ�Χ
    _sub_extent.minX = _glb_extent.minX;
    _sub_extent.maxX = _glb_extent.maxX;
    _sub_extent.maxY = _glb_extent.maxY - layerD._pMetaData->_MBR.minIRow() * _cellSize;
    _sub_extent.minY = _glb_extent.maxY - layerD._pMetaData->_MBR.maxIRow() * _cellSize - _cellSize;
    //cout<<"myrank "<<layerD._pMetaData->myrank<<" "<<_sub_extent.maxY<<" "<<_sub_extent.minY<<endl;
    _xSize = layerD._pMetaData->_localdims.nCols();
    _ySize = layerD._pMetaData->_localdims.nRows();

    _pIDWLayer = &layerD;
    Configure(_pIDWLayer, false);

}
bool IDWOperator::isTermination() {
    return flag;
}

int IDWOperator::searchNbrSamples(const int subMinRow, int cellRow, int cellCol, double* nbrSamples) {
    //double *nbrSamples = new double [_nbrPoints*2];	//���δ�ž��������ֵ��
    int blockRow = (cellRow + subMinRow) / _blockGrain; //ȷ����ǰդ�����ڿ�
    int blockCol = cellCol / _blockGrain;
    int blockRows = _nRows / _blockGrain;
    blockRows += (_nRows % _blockGrain) ? 1 : 0;
    int blockCols = _nCols / _blockGrain;
    blockCols += (_nCols % _blockGrain) ? 1 : 0;

    double cellX = (cellCol + 0.5) * _cellSize + _sub_extent.minX; //��ǰ����ֵդ������
    double cellY = _sub_extent.maxY - (cellRow + 0.5) * _cellSize;
    double maxDist = 0.0; //Ŀǰ��������������ֵ;Ҳ����ܻ���������
    int maxDistIdx = -1; //Ŀǰ����������������������λ��
    int tailIdx = -1; //Ŀǰ������������β�����У��������������������-1
    int searchRad = 0; //�������������뾶;1����3*3����
    //int searchRadLeast = 0;	//�ڴ˲��ϣ���������_nbrPoints������
    bool isSearch = true;
    do {
        //�ռ����������ĺ�ѡblock idx
        int* block2search;
        if (searchRad == 0) {
            block2search = new int[1];
        }
        else {
            block2search = new int[2 * searchRad * 4];
        }
        int blockCount = 0;
        for (int tRow = blockRow - searchRad; tRow <= blockRow + searchRad; ++tRow) {
            if (tRow < 0 || tRow >= blockRows) {
                continue;
            }
            if (tRow == blockRow - searchRad || tRow == blockRow + searchRad) {
                //��ĩ���д�ȫ��
                for (int tCol = blockCol - searchRad; tCol <= blockCol + searchRad; ++tCol) {
                    if (tCol < 0 || tCol >= blockCols) {
                        continue;
                    }
                    block2search[blockCount++] = tRow * blockCols + tCol;
                }
            }
            else {
                //�������������ұ߽�����,���߽�������Ч��棬��Ч������
                if (blockCol - searchRad >= 0 && blockCol - searchRad < blockCols) {
                    block2search[blockCount++] = tRow * blockCols + blockCol - searchRad;
                }
                if (blockCol + searchRad >= 0 && blockCol + searchRad < blockCols) {
                    block2search[blockCount++] = tRow * blockCols + blockCol + searchRad;
                }
            }
        }
        //cout<<"myrank "<<_myRank<<" "<<tRow<<" "<<tCol<<" "<<_pSampleBlocks[tRow*blockCols+tCol].sample_Points.size()<<" "<<endl;
        //����int block2search[2*searchRad*4]�д洢��block idx;
        //��_pSampleBlocks[i].samplePoints��������
        for (int i = 0; i < blockCount; ++i) {
            int blockIdx = block2search[i];
            //cout<<blockIdx<<" "<<_pSampleBlocks[blockIdx].sample_Points.size()<<endl;
            for (vector<Sample_Point>::iterator iter = _pSampleBlocks[blockIdx].sample_Points.begin(); iter != _pSampleBlocks[blockIdx].sample_Points.end(); ++iter) {
                double tmpDist = sqrt((iter->x - cellX) * (iter->x - cellX) + (iter->y - cellY) * (iter->y - cellY));
                if (tailIdx < _nbrPoints - 1) {
                    //�������ĵ㻹����ָ���ĸ�������ֱ�ӷ���β��
                    tailIdx++;
                    nbrSamples[tailIdx * 2] = tmpDist;
                    nbrSamples[tailIdx * 2 + 1] = iter->value;
                    if (tmpDist > maxDist) {
                        maxDist = tmpDist;
                        maxDistIdx = tailIdx;
                    }
                }
                else {
                    tailIdx++;
                    //�Ѿ����㹻�ڽ����㣬�����ѵ����µ������������滻Ŀǰ��Զ�Ǹ���,��������Զ���뼰ID
                    if (tmpDist < maxDist) {
                        nbrSamples[maxDistIdx * 2] = tmpDist;
                        nbrSamples[maxDistIdx * 2 + 1] = iter->value;
                        maxDist = nbrSamples[0];
                        maxDistIdx = 0;
                        //����maxDist,���Ǹ����������ݽṹ������ʡȥ����
                        for (int i = 1; i < _nbrPoints; ++i) {
                            if (nbrSamples[i * 2] > maxDist) {
                                maxDist = nbrSamples[i * 2];
                                maxDistIdx = i;
                            }
                        }
                    }
                }
            }
        }
        if (tailIdx >= _nbrPoints - 1) {
            //if( _myRank==0 ){
            //	cout<<cellRow<<" "<<cellCol<<" "<<(searchRad+0.5)*_cellSize*_blockGrain<<" "<<maxDist<<endl;
            //}
            if ((searchRad + 0.5) * _cellSize * _blockGrain >= maxDist) {
                //maxDist��Խ��ԽС��searchRad��Խ��Խ��
                isSearch = false;
            }
            else {
                ++searchRad;
            }
            //isSearch = false;
        }
        else {
            ++searchRad;
        }
        delete block2search;
    }
    while (isSearch);
    return tailIdx;
}

bool IDWOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& idwL = *(_pIDWLayer->cellSpace());
    const int minRow = _pIDWLayer->_pMetaData->_MBR.minIRow();
    int iRow = coord.iRow();
    int iCol = coord.iCol();
    double startTime, endTime;
    startTime = MPI_Wtime();
    //ÿ���㶼�Ǵ���ֵ�㣬ֻ��������Χ��ͬ����
    double* pNbrSamples = new double [_nbrPoints * 2]; //���δ�ž��������ֵ��
    searchNbrSamples(minRow, iRow, iCol, pNbrSamples); //������ǰդ�������ֵ����nbrSamples�з���
    //�����ֵ���
    double weightSum = 0.0;
    double* pWeight = new double[_nbrPoints];
    for (int i = 0; i < _nbrPoints; ++i) {
        pWeight[i] = 1 / pow(pNbrSamples[i * 2], _idw_power);
        weightSum += pWeight[i];
    }
    idwL[iRow][iCol] = 0;
    for (int i = 0; i < _nbrPoints; ++i) {
        idwL[iRow][iCol] += pNbrSamples[i * 2 + 1] * pWeight[i] / weightSum;
    }
    //if( iRow==20 && iCol==20 && _myRank==0 ){...}	//���Բ�ͬ���ȣ��߽�դ�������������ڽ������Ƿ�һ��

    endTime = MPI_Wtime();
    double time = (endTime - startTime) * 1000;
    if(_myRank==0 && _pComptLayer) {
        (*_pComptLayer->cellSpace())[iRow][iCol] += time;        
    }
    delete pNbrSamples;
    delete pWeight;

    return true;
}
