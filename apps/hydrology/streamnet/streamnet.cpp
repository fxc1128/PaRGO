/***************************************************************************
* pSCA.cpp
*
* Project: GPRO_D8_SCA
* Purpose: Calculation of SCA based on D8 flow direction; Demonstration program for GPRO. 
*
* Author:  Ai Beibei
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
#include "application.h"
#include "SNoperator.h"
#include "communication.h"


using namespace std;
using namespace GPRO;

int main(int argc, char* argv[]) {
    /*  enum ProgramType{MPI_Type = 0,
                   MPI_OpenMP_Type,
                   CUDA_Type,
                   Serial_Type};*/
    Application::START(MPI_Type, argc, argv); //init

    char* inputsrcfile;
	char* inputdemfile;
	char* inputdirfile;
	char* inputareafile;
    char* neighborfile;
    char* outputordfile;
    char* treefile;
	char* coordfile;
	char* outputnetfile;
    char* outputwsfile;
    //int threadNUM;
    if (argc != 11) {
        cerr << "please input right parameter.";
        return 0;
    }
    else {
        inputsrcfile = argv[1];
		inputdemfile = argv[2];
		inputdirfile = argv[3];
		inputareafile = argv[4];
		neighborfile = argv[5];
		outputnetfile = argv[6];
		outputwsfile = argv[7];
		outputordfile = argv[8];
		treefile = argv[9];
		coordfile = argv[10];
		
        //threadNUM = atoi(argv[5]);
    }
    //omp_set_num_threads(threadNUM);
	RasterLayer<double> demLayer("demLayer");
    demLayer.readNeighborhood(neighborfile);
    demLayer.readFile(inputdemfile,ROWWISE_DCMP);//add rowwise_dcmp

	RasterLayer<double> srcLayer("srcLayer");
    srcLayer.readNeighborhood(neighborfile);
    srcLayer.readFile(inputsrcfile,ROWWISE_DCMP);//add rowwise_dcmp
	
	RasterLayer<double> areaLayer("areaLayer");
    areaLayer.readNeighborhood(neighborfile);
    areaLayer.readFile(inputareafile,ROWWISE_DCMP);//add rowwise_dcmp

    RasterLayer<double> dirLayer("dirLayer");
    dirLayer.readNeighborhood(neighborfile);
    dirLayer.readFile(inputdirfile,ROWWISE_DCMP);//add rowwise_dcmp

    RasterLayer<double> wsLayer("wsLayer");
    wsLayer.copyLayerInfo(dirLayer);

	RasterLayer<double> ordLayer("ordLayer");
    ordLayer.copyLayerInfo(dirLayer);

    double starttime;
    double endtime;
    MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();
	
	//tiffIO srcIO(inputsrcfile, SHORT_TYPE);  

    SNOperator snOper;
    snOper.demLayer(dirLayer);
    snOper.srcLayer(srcLayer);
	snOper.areaLayer(areaLayer);
    snOper.dirLayer(dirLayer);
	snOper.coordfile=coordfile;
	snOper.treefile=treefile;
	snOper.streamnetlyr="";
	snOper.streamnetsrc=outputnetfile;
	snOper.elevfile=inputdemfile;
	snOper.elevtype=(DATA_TYPE) FLOAT_TYPE;
	snOper.srcfile=inputsrcfile;
	snOper.srctype=(DATA_TYPE) FLOAT_TYPE;
    snOper.Run();

    MPI_Barrier(MPI_COMM_WORLD);
    endtime = MPI_Wtime();
    cout << "run time is " << endtime - starttime << endl;

    wsLayer.writeFile(outputwsfile);

    Application::END();
    return 0;
}
