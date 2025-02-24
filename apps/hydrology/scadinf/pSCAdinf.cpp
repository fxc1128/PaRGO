/***************************************************************************
* pSCADinf.cpp
*
* Project: GPRO_Dinf_SCA
* Purpose: Calculation of SCA based on Dinf flow direction ; Demonstration program for GPRO. 
*
* Author:  Fan Xingchen
****************************************************************************
* Copyright (c) 2022. Fan Xingchen
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
#include "dinfOperator.h"
#include "communication.h"

using namespace GPRO;

void Usage(const string& error_msg = "") {

    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }
    cout << " Usage: scadinf -dinf <flow direction D-inf file> -nbr <neighbor definition file> -out <output sca file>" << endl;
    cout << " Or use the Simple Usage: scadinf <d-inf file> <neighbor definition file> <output sca file>" << endl << endl;
    cout << "Example.1. scadinf -dinf /path/to/dinf.tif -nbr /path/to/moore.nbr -out /path/to/scadinf.tif" << endl;
    cout << "Example.4. scadinf /path/to/dinf.tif /path/to/moore.nbr /path/to/scadinf.tif" << endl;

    exit(1);
}

int main(int argc, char *argv[]) 
{

	/*!
	 * Parse input arguments.
	 * DO NOT start the application unless the required inputs are provided!
	 */
	if (argc < 4) {
        Usage("Too few arguments to run this program.");
	}
	// Input arguments
	char* inputfilename = nullptr;
	char* neighborfile = nullptr;
	char* outputfilename = nullptr;
	//int dirType;

    int i = 1;
    bool simpleusage = true;
	while (argc > i) {
		if (strcmp(argv[i], "-dinf") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                inputfilename = argv[i];
                i++;
            } else {
	            Usage("No argument followed '-dinf'!");
            }
		} else if (strcmp(argv[i], "-nbr") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                neighborfile = argv[i];
                i++;
			} else {
				Usage("No argument followed '-nbr'!");
			}
        } else if (strcmp(argv[i], "-out") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                outputfilename = argv[i];
                i++;
			} else {
                Usage("No argument followed '-out'!");
			}
        } else { // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputfilename = argv[1];
            neighborfile = argv[2];
            outputfilename = argv[3];
			//dirType = atoi(argv[4]);
			break;
        }
	}
	if (!FileExists(inputfilename)) {
        Usage("The input dir file not exists");
    }
    if (!FileExists(neighborfile)) {
        Usage("neighbor file not exists");
    }
	
	
	Application::START(MPI_Type, argc, argv);
	
	RasterLayer<double> dinfLayer("dinfLayer");
	dinfLayer.readNeighborhood(neighborfile);  
	dinfLayer.readFile(inputfilename, ROWWISE_DCMP);
	
	RasterLayer<double> scaLayer("scaLayer");
    scaLayer.copyLayerInfo(dinfLayer);
	
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();

	DinfOperator dinfOper;	
	dinfOper.dinfLayer(dinfLayer);
	dinfOper.scaLayer(scaLayer);
	dinfOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//if( myrank==0 )
	cout<<"run time is "<<endtime-starttime<<endl;
	scaLayer.writeFile(outputfilename);
	
	Application::END();
	//system("pause");
	return 0;
}