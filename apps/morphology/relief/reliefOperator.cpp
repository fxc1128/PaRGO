#include"reliefOperator.h"


//�������ͼ��
void reliefOperator::
demLayer(RasterLayer<double> &layerD) 
{
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
  cellSize = _pDEMLayer->_pMetaData->cellSize;
  noData = _pDEMLayer->_pMetaData->noData;
    //pWorkBR = &_pDEMLayer->_pMetaData->_localworkBR;
  Configure(_pDEMLayer, false);
}

//������ͼ��
void reliefOperator::reliefLayer(RasterLayer<double> &layerD) 
{
  _pReliefLayer = &layerD;
  Configure(_pReliefLayer,false);
}

//�жϵ���
bool reliefOperator::isTermination()
{
	num--;
	if(num > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//��������ȵ��㷨ʵ��
bool reliefOperator::Operator(const CellCoord &coord)
{
	CellSpace<double> &dem = *(_pDEMLayer->cellSpace());//����ͼ���դ������

	CellSpace<double> &relief = *(_pReliefLayer->cellSpace());//���ͼ���դ������
	  
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);//��������
	
	int iRow = coord.iRow();
	int iCol = coord.iCol();
	
	double d[9];//���3*3���򴰿�
	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	int dCellSize = _pDEMLayer->_pMetaData->cellSize;//DEM��������
	int nodata = _pDEMLayer->_pMetaData->noData;
	int iRow1, iCol1;
	
	//�洢�������ڷ�Χ�ڵ�DEMֵ
	int k = 0;
	int tag=0;
	for(iRow1 = iRow - iNeighborCells; iRow1 <= iRow + iNeighborCells; iRow1++)
	{
		for(iCol1 = iCol - iNeighborCells; iCol1 <= iCol + iNeighborCells; iCol1++)
		{
			d[k] = dem[iRow1][iCol1];
			if(d[k]==nodata)
			{
				tag=1;
			}
			k++;
		}
	}
			if(tag==1)
		{
			relief[iRow][iCol] = nodata;
			return true;
		}
		else
		{
			//����ȣ��߲�㷨
			double max,min;
			max=d[0];	min=d[0];
			for(int m=1;m<9;m++)
			{
				if(d[m]>max)
					max=d[m];
				else if(d[m]<min)
					min=d[m];
			}
			relief[iRow][iCol] =max-min;
			return true;
		}
}

