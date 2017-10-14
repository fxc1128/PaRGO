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
	//int threadNUM;
	if( argc != 10 )
	{
		cout<<"Please Check the input paarameters!"<<endl;
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
	//threadNUM = atoi(argv[9]);

	//omp_set_num_threads(threadNUM);

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
	//equal row dcmp based on region
	idwOper.idwLayer(idwLayer, &spatialrefWkt);	//�Ƚ�idwOperator�����ݳ�Աָ��idwLayerͼ�㣬�ٽ�˴���idwLayer�Ļ���Ԫ����
	//cout<<"idwLayer metadata initialized."<<endl;
	//�������������ʱ���󣬸��ݱ�ͼ���Ԫ����ֱ�ӻ���,�Ƿ���д�����
	starttime = MPI_Wtime();
	idwOper.Run();	//���У����д��idwLayer��cellspace��
	cout<<"idw compute done."<<endl;
	idwLayer.rowWriteFile(outputfilename, true);

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	if (myRank==0)
		cout<<"compute time is "<<endtime-starttime<<endl;

	//MPI_Barrier(MPI_COMM_WORLD);

	////balanced row dcmp based on compute burden
	//vInputLayers[0]->readNeighborhood(dataNeighbor);
	//int* pDcmpIdx = new int[process_nums*4];	//MPI������㲥�Զ������ͣ�ֻ����д�ɽ����к�
	//if( myRank==0 )	//��һ����Ӧ�ý����û����������õ�computLayer��ȥ
	//{
	//	starttime = MPI_Wtime();
	//	vInputLayers[0]->readGlobalFile(vInputnames[0]);	//����ͳ�Ƽ����������ͼ�������������Ԫ����,�ظ�����readFileʱҪ��֤clear������new
	//	vector<RasterLayer<double>* > inputLayers;
	//	inputLayers.push_back( vInputLayers[0] );		//���ж��ͼ���������򹹽�������push_back���
	//	ComputLayer<double> comptLayer( inputLayers, "computLayer" );
	//	comptLayer.readNeighborhood(compuNeighbor);
	//	const int compuSize = 10;	//������ͼ��ֱ���������ͼ���10��,�����û�ָ���������ݶ�Ϊ10
	//	//��ȡ���ؾ���Ļ��֣�������ظ�vDcmpIdx,���ַ�ʽ�ɵڶ�������ָ��
	//	comptLayer.getCompuLoad( pDcmpIdx, ROWWISE_DCMP,compuSize, process_nums );	
	//	comptLayer.writeComptFile(outputfilename);
	//	endtime = MPI_Wtime();
	//	cout<<myRank<<" dcmp time is "<<endtime-starttime<<endl;
	//}
	////MPI_Barrier(MPI_COMM_WORLD);
	////��vDcmpIdx�������̸��������̹㲥�乤���ռ䷶Χ
	//MPI_Bcast(pDcmpIdx,process_nums*4,MPI_INT,0,MPI_COMM_WORLD);
	//CellCoord nwCorner(pDcmpIdx[4*myRank], pDcmpIdx[4*myRank+1]);
	//CellCoord seCorner(pDcmpIdx[4*myRank+2], pDcmpIdx[4*myRank+3]);
	//CoordBR subWorkBR(nwCorner, seCorner);
	//delete []pDcmpIdx;
	////MPI_Barrier(MPI_COMM_WORLD);
	//cout<<myRank<<" "<<subWorkBR<<" "<<subWorkBR.maxIRow() - subWorkBR.minIRow()<<endl;	//rowComptDcmp based on compt-burden has been ok.

	//vInputLayers[0]->readFile(vInputnames[0], subWorkBR, ROWWISE_DCMP);
	//for(int i=1; i<vInputnames.size(); i++){
	//	vInputLayers[i]->readNeighborhood(dataNeighbor);
	//	vInputLayers[i]->readFile(vInputnames[i], subWorkBR, ROWWISE_DCMP);
	//}
	//fcmLayer.copyLayerInfo(*vInputLayers[0]); //�������ͼ��
	//for(int i=0; i<clusterNum; i++){
	//	vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
	//}

	cout<<"write done."<<endl;

	Application::END();
	//system("pause");
	return 0;
}
