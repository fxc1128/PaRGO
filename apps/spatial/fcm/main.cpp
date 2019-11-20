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
#include "fcmOperator.h"
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

	//...
	char* inputfilenames;
	char* dataNeighbor;
	char* compuNeighbor;
	char* outputfilename;
	int clusterNum; //������Ŀ
	double maxIteration; //����������
	double tolerance;//������ֵ
	int wm;//��Ȩָ��
	int domDcmpObj;// 0:SPACE_DIM 2:COMPT_LOAD
	if (argc>0 && argc < 11)
	{
		inputfilenames = argv[1];
		dataNeighbor = argv[2];
		compuNeighbor = argv[3]; 
		outputfilename = argv[4];
		clusterNum =atoi(argv[5]);
		maxIteration =atoi(argv[6]);
		tolerance = atof(argv[7]);
		wm = atof(argv[8]);
		domDcmpObj = atoi(argv[9]);
	}

	int myRank, process_nums;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
	double starttime;
	double endtime;

	//�ַ������������ļ���
	vector<char *> vInputnames;	//�����ļ�����������
	vector<RasterLayer<double> *> vInputLayers;
	//int imageNum;	//�����Ӱ����Ŀ
	char* token = strtok(inputfilenames,",");
	char *filename;
	while(NULL != token)
	{
		filename=token;
		vInputnames.push_back(filename);	//�������ź����ȡ��
		RasterLayer<double> *pLayer = new RasterLayer<double>("none");
		vInputLayers.push_back(pLayer);
		token=strtok(NULL,",");
	}
	if( vInputnames.empty() || clusterNum==0 || maxIteration==0 ){
		return 1;
	}
	RasterLayer<double> fcmLayer("fcmLayer");//�����������ͼ��fcmLayer
	//Ԥ�������ͼ��
	char** pDegLayerName = new char*[clusterNum];
	//string degLayerName;
	vector<RasterLayer<double> *> vDegreeLayer;
	for(int i=0;i<clusterNum;i++)
	{
		pDegLayerName[i] = new char[50];
		sprintf(pDegLayerName[i],"degreeLayer%d.tif",i);	//���������
		RasterLayer<double> *pLayer = new RasterLayer<double>(pDegLayerName[i]);
		vDegreeLayer.push_back(pLayer);
	}

	if(domDcmpObj==0){
		//equal row dcmp based on region
		for(int i=0; i<vInputnames.size(); i++){
			vInputLayers[i]->readNeighborhood(dataNeighbor);
			vInputLayers[i]->readFile(vInputnames[i], ROWWISE_DCMP);
		}
		fcmLayer.copyLayerInfo(*vInputLayers[0]); //�������ͼ��
		for(int i=0; i<clusterNum; i++){
			vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
		}

		RasterLayer<double> comptLayer("fcmLayer");//��ʱ�����ã���׽��ʵ����ǿ�ȣ��Ժ�ķ�װ͸��
		comptLayer.copyLayerInfo(*vInputLayers[0]);

		starttime = MPI_Wtime();
		FCMOperator fcmOper;
		fcmOper.initialization(vInputLayers.size(), clusterNum, maxIteration, tolerance, wm);
		fcmOper.inputLayer(vInputLayers);
		fcmOper.fcmLayer(fcmLayer);
		fcmOper.degLayer(vDegreeLayer);
		fcmOper.comptLayer(comptLayer);//������
		fcmOper.Run();

		char* comptfilename = "D:\\arcgis-data\\pargo\\fcm\\nenjiang_out\\realComputeTime.tif";
		comptLayer.writeFile(comptfilename);	//�����ã�д����׽���ļ���ʱ��
	}else if(domDcmpObj==2){
		////balanced row dcmp based on compute burden

		/// ԭfcm�㷨�汾

		// RasterLayer<double> demLayer("demLayer");
		// RasterLayer<double> pitLayer("pitLayer"); //�������ͼ��
		// demLayer.readNeighborhood(dataNeighbor);
		
		// //vector<CoordBR> vDcmpIdx;	//��Ż���λ����������������ú�㲥��������;MPI������㲥�Զ������ͣ�����;
		// //vDcmpIdx��ʲô�������ͺ��ʣ�����computLayer��ת����ʽ�����û�͸����Ӧ������֤�û�����
		// int* pDcmpIdx = new int[process_nums*4];	//MPI������㲥�Զ������ͣ�ֻ����д�ɽ����к�
		// if( myRank==0 )	//��һ����Ӧ�ý����û����������õ�computLayer��ȥ
		// {
		// 	starttime = MPI_Wtime();
		// 	demLayer.readGlobalFile(inputfilename0);	//����ͳ�Ƽ����������ͼ�������������Ԫ����,�ظ�����readFileʱҪ��֤clear������new
		// 	vector<RasterLayer<double>* > inputLayers;
		// 	inputLayers.push_back( &demLayer );		//���ж��ͼ���������򹹽�������push_back���
			
		// 	ComputLayer<double> comptLayer( inputLayers, "computLayer" );
		// 	comptLayer.readNeighborhood(compuNeighbor);
		// 	const int compuSize = 10;	//������ͼ��ֱ���������ͼ���10��,�����û�ָ���������ݶ�Ϊ10
		// 	//��ȡ���ؾ���Ļ��֣�������ظ�vDcmpIdx,���ַ�ʽ�ɵڶ�������ָ��
		// 	comptLayer.getCompuLoad( pDcmpIdx, ROWWISE_DCMP,compuSize, process_nums );	
		// 	//for( int i=0; i<process_nums*4; i+=4 ){
		// 	//	//it's ok here
		// 	//	cout<<pDcmpIdx[i]<<" "<<pDcmpIdx[i+1]<<" "<<pDcmpIdx[i+2]<<" "<<pDcmpIdx[i+3]<<endl;
		// 	//}
		// 	//comptLayer.writeComptFile(outputfilename);	//��ѡ,Ŀǰ��֧�ִ���д

		// 	endtime = MPI_Wtime();
		// 	cout<<myRank<<" dcmp time is "<<endtime-starttime<<endl;
		// }
		// MPI_Barrier(MPI_COMM_WORLD);
		// //��vDcmpIdx�������̸��������̹㲥�乤���ռ䷶Χ
		// MPI_Bcast(pDcmpIdx,process_nums*4,MPI_INT,0,MPI_COMM_WORLD);
		// CellCoord nwCorner(pDcmpIdx[4*myRank], pDcmpIdx[4*myRank+1]);
		// CellCoord seCorner(pDcmpIdx[4*myRank+2], pDcmpIdx[4*myRank+3]);
		// CoordBR subWorkBR(nwCorner, seCorner);
		// delete []pDcmpIdx;
		// MPI_Barrier(MPI_COMM_WORLD);
		// cout<<myRank<<" "<<subWorkBR<<endl;	//rowComptDcmp based on compt-burden has been ok.

		// demLayer.readFile(inputfilename0, subWorkBR, ROWWISE_DCMP);
		// pitLayer.copyLayerInfo(demLayer);

		/// ���������İ�
		starttime = MPI_Wtime();
		vInputLayers[0]->readNeighborhood(dataNeighbor);
		CoordBR subWorkBR;
		if( myRank==0 ){
			vInputLayers[0]->readGlobalFile(vInputnames[0]);
			vector<RasterLayer<double>* > inputLayers;
			inputLayers.push_back( vInputLayers[0] );//���ж��ͼ���������򹹽�������push_back���
			const int comptGrain = 10;	//������ͼ��ֱ���������ͼ���10��,�����û�ָ��
			ComputLayer<double> comptLayer( inputLayers, comptGrain, "computLayer" );
			comptLayer.readNeighborhood(compuNeighbor);
			comptLayer.newMetaData( comptGrain );
			Transformation<double> transOper( 1, 15, &comptLayer );	//ָ����ֵ��ֵդ��ļ���ǿ��
			transOper.run();	//�������м��㺯��,����comptLayer._pCellSpace
			comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );	//�ں���ͬ��
			comptLayer.writeComptFile("D:\\arcgis-data\\pargo\\fcm\\nenjiang_out\\comp.tif");	// ��ѡ
		}else{
			ComputLayer<double> comptLayer("untitled");
			comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );
			MPI_Barrier(MPI_COMM_WORLD);
		}

		cout<<myRank<<" subWorkBR "<<subWorkBR.minIRow()<<" "<<subWorkBR.maxIRow()<<" "<<subWorkBR.nRows()<<endl;
		endtime = MPI_Wtime();
		if (myRank==0)
			cout<<"dcmp time is "<<endtime-starttime<<endl;

		////������ʵ����ָ������
		//starttime = MPI_Wtime();
		//ComputLayer<double> comptLayer("untitled");
		//char* comptfilename = "/data/aibb/PaRGO_develop/comp.tif";
		//comptLayer.readNeighborhood(compuNeighbor);
		//comptLayer.readFile(comptfilename);
		//CoordBR subWorkBR;
		//comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );
		//cout<<myRank<<" subWorkBR "<<subWorkBR.minIRow()<<" "<<subWorkBR.maxIRow()<<" "<<subWorkBR.nRows()<<endl;
		//endtime = MPI_Wtime();
		//if (myRank==0)
		//	cout<<"dcmp time is "<<endtime-starttime<<endl;

		//��������
		for(int i=0; i<vInputnames.size(); i++){
			vInputLayers[i]->readNeighborhood(dataNeighbor);
			vInputLayers[i]->readFile(vInputnames[i], subWorkBR, ROWWISE_DCMP);
		}
		fcmLayer.copyLayerInfo(*vInputLayers[0]); //�������ͼ��
		for(int i=0; i<clusterNum; i++){
			vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
		}
		//ִ�м���
		starttime = MPI_Wtime();
		FCMOperator fcmOper;
		fcmOper.initialization(vInputLayers.size(), clusterNum, maxIteration, tolerance, wm);
		fcmOper.inputLayer(vInputLayers);
		fcmOper.fcmLayer(fcmLayer);
		fcmOper.degLayer(vDegreeLayer);
		fcmOper.Run();
	}
	//starttime = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	if (myRank==0)
		cout<<"compute time is "<<endtime-starttime<<endl;
	fcmLayer.writeFile(outputfilename);
	//for( size_t i = 0; i < vDegreeLayer.size(); ++i ){
	//	vDegreeLayer[i]->writeFile(pDegLayerName[i]);
	//}
	cout<<"write done."<<endl;

	Application::END();
	//system("pause");
	return 0;
}
