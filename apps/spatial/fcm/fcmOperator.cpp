#include "fcmOperator.h"

FCMOperator::~FCMOperator()
{
	delete centerVal;
	delete centerIndex;
	delete sumNumerator;
	delete sumDenominator;
	delete totNumerator;
	delete totDenominator;
	for (int i = 0; i < clusterNum; i++)
	{
		for (int j = 0; j < _xSize; j++)
		{
			delete[] dist[i][j];
			delete[] degree[i][j];
		}
		delete[] dist[i];
		delete[] degree[i];
	}
	delete[] dist;
	delete[] degree;
}

void FCMOperator::initialization(int iNum, int cNum, int maxIter, double toler, double m)
{
	imageNum = iNum;
	clusterNum = cNum;
	maxIteration = maxIter;
	tolerance = toler;
	wm = m;
}

void FCMOperator::inputLayer(vector<RasterLayer<double> *> layerD)
{
	for (size_t i = 0; i < layerD.size(); ++i)
	{
		_vInputLayer.push_back(layerD[i]);
		Configure(layerD[i], false);
	}
	//_pDEMNbrhood = layerD[0]->nbrhood();

	_noData = layerD[0]->_pMetaData->noData;
	_cellSize = layerD[0]->_pMetaData->cellSize;
	_xSize = layerD[0]->_pMetaData->_localdims.nRows();
	_ySize = layerD[0]->_pMetaData->_localdims.nCols();
	_nRows = layerD[0]->_pMetaData->row;
	_nCols = layerD[0]->_pMetaData->column;
	_rank = layerD[0]->_pMetaData->myrank;
}

void FCMOperator::fcmLayer(RasterLayer<double> &layerD)
{
	_pFCMLayer = &layerD;
	Configure(_pFCMLayer, false);
}

void FCMOperator::degLayer(vector<RasterLayer<double> *> layerD)
{
	for (size_t i = 0; i < layerD.size(); ++i)
	{
		_vDegLayer.push_back(layerD[i]);
		Configure(layerD[i], false);
	}
}

void FCMOperator::comptLayer(RasterLayer<double> &layerD)
{
	_pComptLayer = &layerD;
	Configure(_pComptLayer, false);
}

bool FCMOperator::isTermination()
{
	return flag;
}

void FCMOperator::createRandomIdx(int nums, int range, int *randomIdx)
{
	//��range��Χ�ڲ���nums���������randomIdx����
	srand((unsigned int)time(NULL)); //��ʼ���������
	for (int i = 0; i < nums; i++)
	{
		int tmp = rand() % range; //���������0~n-1
		int j;
		for (j = 0; j < i; j++)
		{
			if (randomIdx[j] == tmp)
				break;
		}
		if (j >= i)
		{ //�²������������ǰ������Ĳ��ظ�
			randomIdx[i] = tmp;
		}
		else
		{
			i--; //������ظ��������²�����ֱ�����ظ�
		}
	}
	////����ʱ���ù̶���������,��Է�5��,��ȷ�������̲���Ϊ�գ���������ж��ٵ���
	//centerIndex[0] = (range/5)+(_ySize/3)*2;
	//centerIndex[1] = (range/3)+(_ySize/3);
	//centerIndex[2] = (range/3)+(_ySize/3)*2;
	//centerIndex[3] = (range/3)*2+(_ySize/3);
	//centerIndex[4] = (range/3)*2+(_ySize/3)*2;
}

void FCMOperator::fnDistance(int curRow, int curCol, double *pInputVal)
{
	for (int i = 0; i < clusterNum; i++)
	{
		dist[i][curRow][curCol] = 0.0;
		for (int j = 0; j < imageNum; j++)
		{
			dist[i][curRow][curCol] += (pInputVal[j] - centerVal[i * imageNum + j]) * (pInputVal[j] - centerVal[i * imageNum + j]);
		}
		dist[i][curRow][curCol] = sqrt(dist[i][curRow][curCol]);
	}
}

void FCMOperator::InitDegree(int curRow, int curCol)
{
	//��ȡ��cell���������ĵ�������
	for (int p = 0; p < clusterNum; p++)
	{
		double sumDistance = 0.0; //ÿһ�������������������ĵľ���֮��
		for (int q = 0; q < clusterNum; q++)
		{
			if (dist[q][curRow][curCol] == 0)
			{
				dist[q][curRow][curCol] = Eps;
			}
			sumDistance += pow((dist[p][curRow][curCol] / dist[q][curRow][curCol]), (2 / (wm - 1)));
		}
		degree[p][curRow][curCol] = (sumDistance == 0.0) ? 1.0 : 1.0 / sumDistance;
		//������ֵ�ж�,�ȵ��ӵ�val��
		subval += pow(degree[p][curRow][curCol], wm) * pow(dist[p][curRow][curCol], 2);
	}
}

bool FCMOperator::Operator(const CellCoord &coord, bool operFlag)
{
	//cout<<"rank"<<_rank<<" ("<<coord.iRow()<<","<<coord.iCol()<<")"<<endl;
	starttime = MPI_Wtime();
	int nRows = _vInputLayer[0]->_pMetaData->row;
	int nCols = _vInputLayer[0]->_pMetaData->column;
	int iRow = coord.iRow();
	int iCol = coord.iCol();

	vector<CellSpace<double> *> vInputL, vDegL;
	for (int i = 0; i < imageNum; ++i)
	{
		vInputL.push_back(_vInputLayer[i]->cellSpace());
	}
	CellSpace<double> &comptL = *(_pComptLayer->cellSpace());
	if (_iterNum == 0)
	{
		comptL[iRow][iCol] = 0.0; //�����˶�����ۣ�Ӱ����Ƶ�׼ȷ�̶�;����������-9999��������󣬶Կ�ֵ�ͷǿ�ֵ��ûӰ��
	}
	if (!((iRow == 1) && (iCol == 1)) && (fabs((*(vInputL[0]))[iRow][iCol] + 9999) <= Eps || fabs((*(vInputL[0]))[iRow][iCol] - _noData) <= Eps) && !((iRow == _xSize - 2) && (iCol == _ySize - 2)))
	{
		endtime = MPI_Wtime();
		comptL[iRow][iCol] += (endtime-starttime)*1000;
		tmpSumTime1 += (endtime - starttime) * 1000;
		return true; //��ֵդ��û��Ҫ������������ֱ������
	}
	for (int i = 0; i < clusterNum; ++i)
	{
		vDegL.push_back(_vDegLayer[i]->cellSpace());
	}
	CellSpace<double> &fcmL = *(_pFCMLayer->cellSpace());

	double *pInputVal = new double[imageNum];
	for (int i = 0; i < imageNum; ++i)
	{
		pInputVal[i] = (*(vInputL[i]))[iRow][iCol]; //*������[]һ������
	}
	if ((iRow == 1) && (iCol == 1))
	{
		if ((_iterNum == 0))
		{
			//��һ�ε���,���ʵ�һ��դ��ʱ����ʼ��������,���������̲����������Ĳ��㲥
			centerVal = new double[clusterNum * imageNum]; //�������� number * imageNum
			centerIndex = new int[clusterNum];			   //������������ clusterNum * 1
			sumNumerator = new double[clusterNum * imageNum];
			sumDenominator = new double[clusterNum];
			totNumerator = new double[clusterNum * imageNum];
			totDenominator = new double[clusterNum];
			//��Ÿ�������������ĵľ���
			dist = new double **[clusterNum]; //dist��ά������*����*����
			for (int i = 0; i < clusterNum; ++i)
			{
				dist[i] = new double *[_xSize];
				for (int j = 0; j < _xSize; ++j)
				{
					dist[i][j] = new double[_ySize];
				}
			}
			degree = new double **[clusterNum]; //����������,��ʼ��Ϊ��ֵ
			for (int i = 0; i < clusterNum; i++)
			{
				degree[i] = new double *[_xSize];
				for (int j = 0; j < _xSize; j++)
				{
					degree[i][j] = new double[_ySize];
					for (int p = 0; p < _ySize; p++)
						degree[i][j][p] = _noData;
				}
			}

			//��һ�ε��������̲��������������,��ʼ�������Ķ��������̵����ݷ�Χ��
			if (_rank == 0)
			{
				cout << "initialization " << _noData << " is ok" << endl;
				//���ҳ����������зǿ�դ��������Щդ�������ȡclusterNum���㣬��Ϊ��ʼ��������
				vector<int> tmpIdx;
				for (int i = 1; i < _xSize - 1; ++i)
				{
					for (int j = 1; j < _ySize - 1; ++j)
					{
						if (fabs((*(vInputL[0]))[i][j] - _noData) > Eps && fabs((*(vInputL[0]))[i][j] + 9999) > Eps)
						{
							tmpIdx.push_back(i * _ySize + j);
						}
					}
				}
				int *randomIdx = new int[clusterNum];
				cout << "cells with value to init cluster center are " << tmpIdx.size() << endl;
				createRandomIdx(clusterNum, tmpIdx.size(), randomIdx); //�����ľ�������ֱ�ӷ���centerIndex
				cout << "init center idx and value is done. " << endl;
				for (int i = 0; i < clusterNum; ++i)
				{
					centerIndex[i] = tmpIdx[randomIdx[i]];
					//�Ȼ�ԭΪ���кţ��ٻ�ȡֵ��
					int tmpRow = centerIndex[i] / _ySize;
					int tmpCol = centerIndex[i] % _ySize; //�Ӳ����������ӳ��ض�Ӧ�����к�,�϶���Ϊ��
					//cout<<tmpRow<<" "<<tmpCol<<" ";
					//��ȡ�������ĸ�ͼ������ֵ
					for (int j = 0; j < imageNum; ++j)
					{
						centerVal[i * imageNum + j] = (*(vInputL[j]))[tmpRow][tmpCol];
						//cout<<centerVal[i*imageNum+j]<<" ";
					}
					//cout<<endl;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(centerVal, clusterNum * imageNum, MPI_DOUBLE, 0, MPI_COMM_WORLD); //����0����㲥��������
		}
		else
		{
			//���ε�����ʼʱ���㲥���µľ�������
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(centerVal, clusterNum * imageNum, MPI_DOUBLE, 0, MPI_COMM_WORLD); //����0����㲥��������
		}
	}

	//starttime = MPI_Wtime();	//����ʱ�䲶׽
	if (fabs(pInputVal[0] + 9999) <= Eps || fabs(pInputVal[0] - _noData) <= Eps)
	{
		//��ֵ�����������һ������ֵ
	}
	else
	{
		//��nodata��cell���������ĵľ���ֵ����
		fnDistance(iRow, iCol, pInputVal);
		//�����ȼ���
		InitDegree(iRow, iCol);
	}
	endtime = MPI_Wtime(); //����ʱ�䲶׽
	comptL[iRow][iCol] += (endtime-starttime)*1000;
	tmpSumTime1 += (endtime - starttime) * 1000;
	if ((iRow == _xSize - 2) && (iCol == _ySize - 2))
	{
		//MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&subval, &totval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		_iterNum++;
		subval = 0.0;
		if (_rank == 0)
		{
			cout << "_iterNum is " << _iterNum << endl;
			//cout<<"totval "<<totval<<endl;
			//cout<<"oldtval "<<oldtval<<endl;
			//cout<<"fabs(totval-oldtval) "<<fabs(totval-oldtval)<<endl;
		}
		//��ֵ�ж�,�����һ��դ������val��Լ������0���ж����ξ��Բ��Ƿ�����ֵ��Χ��
		if ((fabs(totval - oldtval) <= tolerance) || (_iterNum >= maxIteration))
		{
			//����,ȷ��ÿ��cell����������࣬���������Ÿ���fcm[i][j]
			for (int i = 1; i < _xSize - 1; i++)
			{
				for (int j = 1; j < _ySize - 1; j++)
				{
					starttime = MPI_Wtime();
					if ((*(vInputL)[0])[i][j] != _noData)
					{
						int cNum = -1;			//�������
						double degreeMax = 0.0; //���������ֵ
						for (int p = 0; p < clusterNum; p++)
						{
							if (degree[p][i][j] > degreeMax)
							{
								degreeMax = degree[p][i][j];
								cNum = p;
							}
							//��ֵ�Ը���������ȸ�degreeLayer
							(*vDegL[p])[i][j] = degree[p][i][j];
							//Ŀǰ�����л������������Ϣ��������������֪��Ӱ���ж��
							if ((degree[p][i][j] - _noData) > Eps)
							{ //���ж�Ӧ���������
								partitionCoef += degree[p][i][j] * degree[p][i][j] / (nRows * nCols - nodataNums);
								if (degree[p][i][j] > 0)
								{
									entropy += -1.0 * degree[p][i][j] * log(degree[p][i][j]) / (nRows * nCols - nodataNums);
								}
							}
						}
						fcmL[i][j] = cNum;
					}
					else
					{
						fcmL[i][j] = _noData;
						for (int p = 0; p < clusterNum; p++)
						{
							(*vDegL[p])[i][j] = _noData;
						}
					}
					endtime = MPI_Wtime();
					comptL[iRow][iCol] += (endtime-starttime)*1000;
					tmpSumTime1 += (endtime - starttime) * 1000;
				}
			}
			MPI_Allreduce(&partitionCoef, &totpartitionCoef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&entropy, &totentropy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			if (_rank == 0)
				cout << "totpartitionCoef " << totpartitionCoef << " totentropy " << totentropy << " nodataNums " << nodataNums << endl;
			cout << _rank << " time is " << tmpSumTime1 << endl;
		}
		else
		{
			oldtval = totval;
			for (int p = 0; p < clusterNum; p++)
			{
				double tmpSum = 0;
				//������ӷ�ĸ
				sumDenominator[p] = 0.0;
				for (int q = 0; q < imageNum; q++)
				{
					sumNumerator[p * imageNum + q] = 0.0;
				}
				//������ӷ�ĸ,���¾�������
				int valCount = 0;
				for (int i = 1; i < _xSize - 1; i++)
				{ //һ��Ҫע������ֻ������Ч�ռ�
					for (int j = 1; j < _ySize - 1; j++)
					{
						starttime = MPI_Wtime();
						if (fabs((*vInputL[0])[i][j] - _noData) > Eps && fabs((*vInputL[0])[i][j] + 9999) > Eps)
						{
							sumDenominator[p] += pow(degree[p][i][j], wm);
							tmpSum += degree[p][i][j];
							for (int q = 0; q < imageNum; q++)
							{
								sumNumerator[p * imageNum + q] += (pow(degree[p][i][j], wm) * (*vInputL[q])[i][j]);
							}
							++valCount;
						}
						endtime = MPI_Wtime();
						//ScomptL[iRow][iCol] += (endtime-starttime)*1000;
						tmpSumTime1 += (endtime - starttime) * 1000;
					}
				}
			}
			MPI_Allreduce(sumDenominator, totDenominator, clusterNum, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(sumNumerator, totNumerator, clusterNum * imageNum, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			//���¾�������center
			for (int p = 0; p < clusterNum; p++)
			{
				for (int q = 0; q < imageNum; q++)
				{
					centerVal[p * imageNum + q] = totNumerator[p * imageNum + q] / totDenominator[p];
				}
			}
			Termination = 0;
		}
	}

	delete pInputVal;

	endtime = MPI_Wtime();
	comptL[iRow][iCol] += (endtime - starttime) * 1000;
	if (comptL[iRow][iCol] > 30)
	{
		cout<<"("<<iRow<<","<<iCol<<":"<<comptL[iRow][iCol]<<") after compute.";
	}

	return true;
}