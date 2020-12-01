#include "fcmOperator.h"

FCMOperator::~FCMOperator() {
    delete centerVal;
    delete centerIndex;
    delete sumNumerator;
    delete sumDenominator;
    delete totNumerator;
    delete totDenominator;
    for (int i = 0; i < clusterNum; i++) {
        for (int j = 0; j < _xSize; j++) {
            delete[] dist[i][j];
            delete[] degree[i][j];
        }
        delete[] dist[i];
        delete[] degree[i];
    }
    delete[] dist;
    delete[] degree;
}

void FCMOperator::initialization(int iNum, int cNum, int maxIter, double toler, double m) {
    imageNum = iNum;
    clusterNum = cNum;
    maxIteration = maxIter;
    tolerance = toler;
    weight = m;
}

void FCMOperator::inputLayer(vector<RasterLayer<double> *> layerD) {
    for (size_t i = 0; i < layerD.size(); ++i) {
        _vInputLayer.push_back(layerD[i]);
        Configure(layerD[i], false);
    }
    _noData = layerD[0]->_pMetaData->noData;
    _cellSize = layerD[0]->_pMetaData->cellSize;
    _xSize = layerD[0]->_pMetaData->_localdims.nRows();
    _ySize = layerD[0]->_pMetaData->_localdims.nCols();
    _nRows = layerD[0]->_pMetaData->row;
    _nCols = layerD[0]->_pMetaData->column;
    _rank = layerD[0]->_pMetaData->myrank;
    _pComptLayer = nullptr;
}

void FCMOperator::fcmLayer(RasterLayer<double>& layerD) {
    _pFCMLayer = &layerD;
    Configure(_pFCMLayer, false);
}

void FCMOperator::degLayer(vector<RasterLayer<double> *> layerD) {
    for (size_t i = 0; i < layerD.size(); ++i) {
        _vDegLayer.push_back(layerD[i]);
        Configure(layerD[i], false);
    }
}

bool FCMOperator::isTermination() {
    return flag;
}


bool FCMOperator::isNoData(double value, int rasterLayerIndex){
    return fabs(value - _vInputLayer[rasterLayerIndex]->metaData()->noData) <Eps;
}

bool FCMOperator::allNoDataAt(int row, int col) {
    for (int i=0;i<_vInputLayer.size();i++) {
        if(!isNoData((*_vInputLayer[i]->cellSpace())[row][col],i)) {
            return false;
        }
    }
    return true;
}

void FCMOperator::createRandomIdx(int nums, int range, int* randomIdx) {
    //在range范围内产生nums个随机数，randomIdx返回
    srand((unsigned int)time(NULL)); //初始化随机种子
    for (int i = 0; i < nums; i++) {
        int tmp = rand() % range; //产生随机数0~n-1
        int j;
        for (j = 0; j < i; j++) {
            if (randomIdx[j] == tmp)
                break;
        }
        if (j >= i) {
            //新产生的随机数与前面产生的不重复
            randomIdx[i] = tmp;
        }
        else {
            i--; //如果有重复，则重新产生，直至不重复
        }
    }
    ////测试时暂用固定聚类中心,针对分5类,已确定单进程不会为空，多进程需判断再调整
    //centerIndex[0] = (range/5)+(_ySize/3)*2;
    //centerIndex[1] = (range/3)+(_ySize/3);
    //centerIndex[2] = (range/3)+(_ySize/3)*2;
    //centerIndex[3] = (range/3)*2+(_ySize/3);
    //centerIndex[4] = (range/3)*2+(_ySize/3)*2;
}

void FCMOperator::fnDistance(int curRow, int curCol, double* pInputVal) {
    for (int i = 0; i < clusterNum; i++) {
        dist[i][curRow][curCol] = 0.0;
        int count = 0;
        for (int j = 0; j < imageNum; j++) {
            if (isNoData(pInputVal[j], j)) {
                cout<<"found a nodata! ("<< curRow<<","<<curCol<<")"<<endl;
                continue;
            }
            dist[i][curRow][curCol] += (pInputVal[j] - centerVal[i * imageNum + j]) * (pInputVal[j] - centerVal[i *imageNum + j]);
            count++;//2020-10-21更新，多层时计算平均值
        }
        dist[i][curRow][curCol] = sqrt(dist[i][curRow][curCol]/count);
    }
}

void FCMOperator::InitDegree(int curRow, int curCol) {
    //求取该cell到各类中心的隶属度
    for (int p = 0; p < clusterNum; p++) {
        double sumDistance = 0.0; //每一个样本到各个聚类中心的距离之和
        for (int q = 0; q < clusterNum; q++) {
            if (dist[q][curRow][curCol] == 0) {
                dist[q][curRow][curCol] = Eps;
            }
            sumDistance += pow(dist[p][curRow][curCol] / dist[q][curRow][curCol], 2 / (weight - 1));
        }
        degree[p][curRow][curCol] = sumDistance == 0.0 ? 1.0 : 1.0 / sumDistance;
        //用于阈值判断,先叠加到val中
        subval += pow(degree[p][curRow][curCol], weight) * pow(dist[p][curRow][curCol], 2);
    }
}

void FCMOperator::initRandomClusterCenters(double* clusterCenters) {
    cout << "initialization " << _noData << " is ok" << endl;
    //先找出本进程所有非空栅格，再在这些栅格里随机取clusterNum个点，作为初始聚类中心
    vector<int> tmpIdx;
    for (int i = 1; i < _xSize - 1; ++i) {
        for (int j = 1; j < _ySize - 1; ++j) {
            const double value = (*_vInputLayer[0]->cellSpace())[i][j];
            if (!isNoData(value,0)) {
                tmpIdx.push_back(i * _ySize + j);
            }
        }
    }
    int* randomIdx = new int[clusterNum];
    cout << "cells with value to init cluster center are " << tmpIdx.size() << endl;
    createRandomIdx(clusterNum, tmpIdx.size(), randomIdx); //产生的聚类中心直接放入centerIndex
    cout << "init center idx and value is done. " << endl;
    for (int i = 0; i < clusterNum; ++i) {
        centerIndex[i] = tmpIdx[randomIdx[i]];
        //先还原为行列号，再获取值等
        int tmpRow = centerIndex[i] / _ySize;
        int tmpCol = centerIndex[i] % _ySize; //从产生的随机数映射回对应的行列号,肯定不为空
        //cout<<tmpRow<<" "<<tmpCol<<" ";
        //获取聚类中心各图层属性值
        for (int j = 0; j < imageNum; ++j) {
            centerVal[i * imageNum + j] = (*_vInputLayer[j]->cellSpace())[tmpRow][tmpCol];
            //cout<<centerVal[i*imageNum+j]<<" ";
        }
        //cout<<endl;
    }
}

///归类,确定每个cell的最大隶属类，并将该类编号赋给fcmL[i][j]
void FCMOperator::assignMaxMembershipDegrees() {
    double startTime=0,endTime=0;
    CellSpace<double>& fcmL = *(_pFCMLayer->cellSpace());
    for (int i = 1; i < _xSize - 1; i++) {
        for (int j = 1; j < _ySize - 1; j++) {
            startTime = MPI_Wtime();
            if (!allNoDataAt(i,j)) {
                int cNum = -1; //所属类号
                double degreeMax = 0.0; //最大隶属度值
                for (int p = 0; p < clusterNum; p++) {
                    if (degree[p][i][j] > degreeMax) {
                        degreeMax = degree[p][i][j];
                        cNum = p;
                    }
                    //赋值对各类的隶属度给degreeLayer
                    (*_vDegLayer[p]->cellSpace())[i][j] = degree[p][i][j];
                    //目前测试中还是输出了熵信息，如果不输出，不知道影响有多大
                    if ((degree[p][i][j] - _vDegLayer[p]->metaData()->noData) > Eps) {
                        //这判断应该是冗余的
                        int nRows = _vInputLayer[0]->_pMetaData->row;
                        int nCols = _vInputLayer[0]->_pMetaData->column;
                        partitionCoef += degree[p][i][j] * degree[p][i][j] / (nRows * nCols - nodataNums);
                        if (degree[p][i][j] > 0) {
                            entropy += -1.0 * degree[p][i][j] * log(degree[p][i][j]) / (nRows * nCols -
                                nodataNums);
                        }
                    }
                }
                fcmL[i][j] = cNum;
            }else {
                fcmL[i][j] = _noData;
                for (int p = 0; p < clusterNum; p++) {
                    (*_vDegLayer[p]->cellSpace())[i][j] = _noData;
                }
            }
        }
    }
}

bool FCMOperator::Operator(const CellCoord& coord, bool operFlag) {
    //cout<<"rank"<<_rank<<" ("<<coord.iRow()<<","<<coord.iCol()<<")"<<endl;
    double time;
    int iRow = coord.iRow();
    int iCol = coord.iCol();
    if (_pComptLayer) {
        if (_iterNum == 0) {
            (*_pComptLayer->cellSpace())[iRow][iCol] = 0.0; //引入了额外代价，影响估计的准确程度;但不这样对-9999会计算有误，对空值和非空值都没影响
        }
    }

    double* pInputVal = new double[imageNum];
    for (int i = 0; i < imageNum; ++i) {
        pInputVal[i] = (*_vInputLayer[i]->cellSpace())[iRow][iCol]; //*二级，[]一级优先
    }
    if ((iRow == 1) && (iCol == 1)) {
        if ((_iterNum == 0)) {
            //第一次迭代,访问第一个栅格时，初始化各变量,并由主进程产生聚类中心并广播
            centerVal = new double[clusterNum * imageNum]; //聚类中心 number * imageNum
            centerIndex = new int[clusterNum]; //聚类中心索引 clusterNum * 1
            sumNumerator = new double[clusterNum * imageNum];
            sumDenominator = new double[clusterNum];
            totNumerator = new double[clusterNum * imageNum];
            totDenominator = new double[clusterNum];
            //存放各个点与各类中心的距离
            dist = new double **[clusterNum]; //dist三维，类数*行数*列数
            for (int i = 0; i < clusterNum; ++i) {
                dist[i] = new double *[_xSize];
                for (int j = 0; j < _xSize; ++j) {
                    dist[i][j] = new double[_ySize];
                }
            }
            degree = new double **[clusterNum]; //隶属度数组,初始化为空值
            for (int i = 0; i < clusterNum; i++) {
                degree[i] = new double *[_xSize];
                for (int j = 0; j < _xSize; j++) {
                    degree[i][j] = new double[_ySize];
                    for (int p = 0; p < _ySize; p++)
                        degree[i][j][p] = _noData;
                }
            }

            //第一次迭代主进程产生随机聚类中心,初始聚类中心都在主进程的数据范围内
            if (_rank == 0) {
                initRandomClusterCenters(centerVal);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(centerVal, clusterNum * imageNum, MPI_DOUBLE, 0, MPI_COMM_WORLD); //进程0负责广播聚类中心
        }
        else {
            //本次迭代开始时，广播更新的聚类中心
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(centerVal, clusterNum * imageNum, MPI_DOUBLE, 0, MPI_COMM_WORLD); //进程0负责广播聚类中心
        }
        starttime=MPI_Wtime();
    }
    
    time = MPI_Wtime();

    if (fabs(pInputVal[0] + 9999) <= Eps || fabs(pInputVal[0] - _noData) <= Eps) {
        //空值不做处理，最后一并赋空值
    }
    else {
        fnDistance(iRow, iCol, pInputVal); //非nodata的cell到各类中心的距离值计算
        InitDegree(iRow, iCol); //隶属度计算
    }
    // if(allNoDataAt(iRow,iCol)) {
    //     for (int i = 0; i < clusterNum; i++) {
    //         dist[i][iRow][iCol] = -1;
    //     }
    // }else{
    //     fnDistance(iRow, iCol, pInputVal); //cell到各类中心的距离值计算
    //     InitDegree(iRow, iCol); //隶属度计算
    // }

    if (_pComptLayer) (*_pComptLayer->cellSpace())[iRow][iCol] += (MPI_Wtime() - time) * 1000;


    if ((iRow == _xSize - 2) && (iCol == _ySize - 2)) {
        //MPI_Barrier(MPI_COMM_WORLD);
        cout<<"rank"<<_rank<<" compute time before reduce is "<<MPI_Wtime()-starttime<<"s"<<endl;
        time = MPI_Wtime();
        MPI_Allreduce(&subval, &totval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        cout<<"rank"<<_rank<<" reduce1 time is "<<MPI_Wtime()-time<<"s"<<endl;
        _iterNum++;
        subval = 0.0;
        if (_rank == 0) {
            cout << "_iterNum is " << _iterNum << endl;
        }
        //阈值判断,到最后一个栅格计算后，val归约到进程0，判断两次绝对差是否在阈值范围内
        if ((fabs(totval - oldtval) <= tolerance) || (_iterNum >= maxIteration)) {
            
            double tmpStart = MPI_Wtime();
            assignMaxMembershipDegrees();
            time=MPI_Wtime();
            MPI_Allreduce(&partitionCoef, &totpartitionCoef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&entropy, &totentropy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (_rank == 0) {
                double tmpEnd = MPI_Wtime();
                cout << "totpartitionCoef " << totpartitionCoef << " totentropy " << totentropy << " nodataNums " << nodataNums << endl;
                cout << _rank << " time is " << tmpEnd-tmpStart << endl;
            }
            cout<<"rank"<<_rank<<" reduce2 time is "<<MPI_Wtime()-time<<"s"<<endl;

        }
        else {
            time=MPI_Wtime();
            oldtval = totval;
            for (int p = 0; p < clusterNum; p++) {
                //归零分子分母
                sumDenominator[p] = 0.0;
                for (int q = 0; q < imageNum; q++) {
                    sumNumerator[p * imageNum + q] = 0.0;
                }
                //计算分子分母,求新聚类中心
                int valCount = 0;
                for (int i = 1; i < _xSize - 1; i++) {
                    //一定要注意这里只计算有效空间
                    for (int j = 1; j < _ySize - 1; j++) {
                        //if (fabs((*_vInputLayer[0]->cellSpace())[i][j] - _noData) > Eps && fabs((*_vInputLayer[0]->cellSpace())[i][j] + 9999) > Eps) {
                        //    sumDenominator[p] += pow(degree[p][i][j], weight);
                        //    for (int q = 0; q < imageNum; q++) {
                        //        sumNumerator[p * imageNum + q] += (pow(degree[p][i][j], weight) * (*_vInputLayer[q]->cellSpace())[i][j]);
                        //    }
                        //    ++valCount;
                        //}

                        //如果也计算含nodata的多图层呢
                        if(allNoDataAt(i,j)) {
                            continue;
                        }
                        sumDenominator[p] += pow(degree[p][i][j], weight);
                        for (int q = 0; q < imageNum; q++)
                        {
                            const double value = (*_vInputLayer[q]->cellSpace())[i][j];
                            if (isNoData(value,q)) {
                                continue;
                            }
                            sumNumerator[p * imageNum + q] += pow(degree[p][i][j], weight) * value;
                            ++valCount;
                        }
                    }
                }
            }
            cout<<"rank"<<_rank<<" disp time is "<<MPI_Wtime()-time<<"s"<<endl;
            
            time=MPI_Wtime();
            MPI_Allreduce(sumDenominator, totDenominator, clusterNum, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(sumNumerator, totNumerator, clusterNum * imageNum, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            cout<<"rank"<<_rank<<" reduce3 time is "<<MPI_Wtime()-time<<"s"<<endl;
            //更新聚类中心center
            for (int p = 0; p < clusterNum; p++) {
                for (int q = 0; q < imageNum; q++) {
                    centerVal[p * imageNum + q] = totNumerator[p * imageNum + q] / totDenominator[p];
                }
            }
            Termination = 0;
        }
    }

    delete pInputVal;


    //if (_pComptLayer && 
    //    !(iRow == _xSize - 2 && iCol == _ySize - 2) &&
    //    !(iRow == 1 && iCol == 1)) {
    //    (*_pComptLayer->cellSpace())[iRow][iCol] += (endTime - startTime) * 1000;
    //}
    //*_pComptLayer->cellSpace()[iRow][iCol] += (endTime - startTime) * 1000;
    // if (comptL[iRow][iCol] > 30)
    // {
    // 	cout<<"("<<iRow<<","<<iCol<<":"<<comptL[iRow][iCol]<<") after compute.";
    // }

    return true;
}
