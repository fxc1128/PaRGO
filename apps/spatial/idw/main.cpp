/***************************************************************************
* main.cpp
*
* Project: PaRGO_0210
* Purpose: developing computation load balancing decomposition for GPRO. 
*			
* Author:  Zhan Lijun;Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2017. Ai Beibei
* 
****************************************************************************/


#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <sstream>
#include <omp.h>
#include "mpi.h"
#include "neighborhood.h"
#include "cellSpace.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "computLayer.h"
#include "application.h"
#include "idwOperator.h"
#include "communication.h"
#include "deComposition.h"
#include "transformation.h"

using namespace std;
using namespace GPRO;


int main(int argc, char *argv[]) 
{
	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type,
				   Serial_Type};*/
	/*  enum DomDcmpType{NON_DCMP = 0,
		ROWWISE_DCMP,
		COLWISE_DCMP,
		BLOCK_DCMP};*/
	/*  enum DomDcmpObj{SPACE_DIM = 0,
				   DATA_LOAD,
				   COMPT_LOAD};*/
	Application::START(MPI_Type, argc, argv); //init

	char *samplefilename, *outputfilename;
	char *dataNeighbor, *compuNeighbor;
	float cellSize;
	int fldIdx, idw_nbrPoints, idw_power, idw_buffer;	//��ʱ������Ϊint�ɸ�Ϊ������
	int dcmpObjType;
	if( argc != 11 )
	{
		cout<<"Please Check the input parameters!"<<endl;
		return 0;
	}
	samplefilename = argv[1];
	outputfilename = argv[2];
	dataNeighbor = argv[3];	//1*1����
	compuNeighbor = argv[4]; 
	cellSize = atof(argv[5]);	//����ֵդ��ֱ���
	fldIdx = atoi(argv[6]);	//ʸ����������ֵ������
	idw_power = atoi(argv[7]);	//�������Ȩ�ݣ�ͨ��ȡ2
	idw_nbrPoints = atoi(argv[8]);	//�����ڽ�����
	idw_buffer = atof(argv[9]);	//��������뾶
	dcmpObjType = atoi(argv[10]);

	int myRank, process_nums;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
	double starttime;
	double endtime;

	int blockGrain = 100;	//granularity,�����Կ��ŵĴ��������ȣ���դ��ֱ���Ϊ������λ���û���������ָ��
	IDWOperator idwOper(cellSize, idw_nbrPoints, idw_power, idw_buffer, blockGrain);
	char* spatialrefWkt;	//ͶӰ��Ϣ
	int sample_nums;
	sample_nums = idwOper.readSampleNums( samplefilename, &spatialrefWkt );	//��ȡ������Ŀ��idwOper.sample_nums
	double **pAllSamples=(double **)malloc(sample_nums*sizeof(double *));//Ϊ��������������洢�ռ䣬�����ӿ����ͷ�
	for (int k=0; k<sample_nums; k++)
		pAllSamples[k]=(double *)malloc(3*sizeof(double));
	if (pAllSamples==NULL)
	{
		cout<<"Faliure memory request!"<<endl;
		return 0;
	}
	idwOper.readSamples( samplefilename, fldIdx, &spatialrefWkt, pAllSamples );	//��ȡ���㣬��������idwOper.glb_extent
	//�ɻ�ȡidwLayer�����귶Χ,����idwOper.sample_extent
	idwOper.creatSampleBlocks(pAllSamples);	//����pAllSamples���ֿ����idwOper._pSampleBlocks��Ա
	//cout<<"creatSampleBlocks() done."<<endl;

	//�Դ�������ʽ��֯���㣬���ݳ�Ա��������ÿ��դ������һϵ������
	RasterLayer<double> idwLayer("idwLayer");
	idwLayer.readNeighborhood(dataNeighbor);

	if(dcmpObjType==0){
		//equal row dcmp based on region
		idwOper.idwLayer(idwLayer, &spatialrefWkt);	//�Ƚ�idwOperator�����ݳ�Աָ��idwLayerͼ�㣬�ٽ�˴���idwLayer�Ļ���Ԫ����
		//�������������ʱ���󣬸��ݱ�ͼ���Ԫ����ֱ�ӻ���,�Ƿ���д�����
		starttime = MPI_Wtime();
		idwOper.Run();	//���У����д��idwLayer��cellspace��
	}
	else if(dcmpObjType==2){
		//balanced row dcmp based on compute burden
		idwOper.idwLayer(idwLayer, &spatialrefWkt);	//�Ƚ�idwOperator�����ݳ�Աָ��idwLayerͼ�㣬�ٽ�˴���idwLayer�Ļ���Ԫ����
		
		idwLayer._pMetaData->_domDcmpType = NON_DCMP;	//wyj 2019-11-12:��ʱд������
		CoordBR subWorkBR;
		if( myRank==0 )	//��һ����Ӧ�ý����û����������õ�computLayer��ȥ
		{
			starttime = MPI_Wtime();
			ComputLayer<double> comptLayer(&idwLayer,200,"computLayer");
			comptLayer.readNeighborhood(compuNeighbor);
			const int compuSize = 10;	//������ͼ��ֱ���������ͼ���10��,�����û�ָ���������ݶ�Ϊ10
			comptLayer.newMetaData( compuSize );
			Transformation<double> transOper( 1, 15, &comptLayer );	//ָ����ֵ��ֵդ��ļ���ǿ��
			transOper.run();	//�������м��㺯��,����comptLayer._pCellSpace
			comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);	
			comptLayer.writeComptFile("D:\\arcgis-data\\pargo\\idw\\output\\comp.tif");
			endtime = MPI_Wtime();
			cout<<myRank<<" dcmp time is "<<endtime-starttime<<endl;
		}else{
			ComputLayer<double> comptLayer("untitled");
			comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );
			MPI_Barrier(MPI_COMM_WORLD);
		}

		//cout<<"idwLayer metadata initialized."<<endl;
		//�������������ʱ���󣬸��ݱ�ͼ���Ԫ����ֱ�ӻ���,�Ƿ���д�����
		starttime = MPI_Wtime();
		idwOper.Run();	//���У����д��idwLayer��cellspace��
	}
	//cout<<"idw compute done."<<endl;
	idwLayer.writeFile(outputfilename);

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	if (myRank==0)
		cout<<"compute time is "<<endtime-starttime<<endl;

	cout<<"write done."<<endl;

	Application::END();
	//system("pause");
	return 0;
}
