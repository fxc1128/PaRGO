#ifndef SNOPERATOR_H
#define SNOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <queue>
#include <stdint.h>
#include "ogr_api.h"
#include <gdal.h>
#include <cpl_conv.h>
#include <cpl_string.h>
#include <ogr_spatialref.h>
#include "mpi.h"

#define Eps 0.000001
#define MAX_STRING_LENGTH 255
#define MAXLN 4096
using namespace GPRO;
long LAST_ID = -1;
enum DATA_TYPE
{ SHORT_TYPE,
LONG_TYPE,
FLOAT_TYPE
};
class SNOperator : public RasterOperator<double> 
{
public:

	char* elevfile;
	DATA_TYPE elevtype;
	char* srcfile;
	DATA_TYPE srctype;
	//tiffIO tiffio;
	SNOperator()
		:RasterOperator<double>(),
		_pDEMLayer(0), _pareaLayer(0), _pDirLayer(0),_psrcLayer(0),num(0), _maxRow(0), _maxCol(0),Sup(0),Sdown(0),
		_contribs("contribs"),_wsgrid("wsgrid"),_idgrid("idgrid"),_lengths("lengths") {}

	~SNOperator() {}
	//char *elevfile;
	//DATA_TYPE elevtype;
	//tiffIO *elevIO(char *elevfile, DATA_TYPE elevtype);
	void demLayer(RasterLayer<double> &layerD);
	void srcLayer(RasterLayer<double> &layerD);
	void dirLayer(RasterLayer<double> &layerD);
	void areaLayer(RasterLayer<double> &layerD);

	//void tiffio(char &file, DATA_TYPE &elevtype);
	virtual bool isTermination();
	virtual bool Operator(const CellCoord &coord,bool operFlag);

	const char *treefile;
	const char *coordfile;
	char *streamnetsrc;
	char *streamnetlyr;

protected:
	double _cellSize;
	double _noData;
	int num;
	int _maxRow;
	int _maxCol;
	int _rank;
	long Sup,Sdown;
	RasterLayer<double> _contribs;
	RasterLayer<double> _idgrid;
	RasterLayer<double> _wsgrid;
	RasterLayer<double> _lengths;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pDirLayer;
	RasterLayer<double> *_psrcLayer;
	RasterLayer<double> *_pareaLayer;
	Neighborhood<double> *_pDEMNbrhood;





};

//  Parameters for WGS84 assumed for all geographic coordinates
const double elipa=6378137.000;
const double elipb=6356752.314;
const double boa=elipb/elipa;

class tiffIO{
private:
	GDALDatasetH fh;			//gdal file handle
	GDALDatasetH copyfh;
	int isFileInititialized;
	GDALDriverH hDriver;
	GDALRasterBandH bandh;    
	// gdal band
	int rank, size;			//MPI rank & size, rank=number for this process, size=number of processes
	uint32_t totalX;		//DGT	// unsigned long BT - ??width of entire grid in number of cells (all partitions)
	uint32_t totalY;		//DGT	// unsigned long BT - ??length of entire grid in number of cells (all partitions)
	//double dx;				//??width of each cell
	//double dy;				//??length of each cell
	double xllcenter;		//horizontal center point of lower left grid cell in grographic coordinates, not grid coordinates
	double yllcenter;		//vertical center point of lower left grid cell in geographic coordinates, not grid coordinates
	double xleftedge;		//horizontal coordinate of left edge of grid in geographic coordinates, not grid coordinates
	double ytopedge;		//vertical coordinate of top edge of grid in geographic coordinates, not grid coordinates
	DATA_TYPE datatype;		//datatype of the grid values and the nodata value: short, long, or float
	//void *nodata;    //pointer to the nodata value, the nodata value type is indicated by datatype
	double nodata;	// noDatarefactor 11/18/17		
	void *filenodata;       //pointer to no data value from the file.  This may be different from nodata because filedatatype and datatype are not equivalent 
	char filename[MAXLN];  //  Save filename for error or warning writes
	const char *valueUnit; //value units
	double *Yp;
	double *dxc,*dyc;
	double dxA,dyA,dlat,dlon,xleftedge_g,ytopedge_g,xllcenter_g,yllcenter_g;
	int IsGeographic;
	OGRSpatialReferenceH  hSRS;

	//  Mappings


public:
	tiffIO(char *fname, DATA_TYPE newtype);
	tiffIO(char *fname, DATA_TYPE newtype, double nodata, const tiffIO &copy); // noDatarefactor 11/18/17
	//tiffIO(char *fname, DATA_TYPE newtype, void* nd, const tiffIO &copy); 
	~tiffIO();

	//BT void read(unsigned long long xstart, unsigned long long ystart, unsigned long long numRows, unsigned long long numCols, void* dest);
	//BT void write(unsigned long long xstart, unsigned long long ystart, unsigned long long numRows, unsigned long long numCols, void* source);
	void read(long xstart, long ystart, long numRows, long numCols, void* dest);
	void write(long xstart, long ystart, long numRows, long numCols, void* source);

	bool compareTiff(const tiffIO &comp);

	//void geoToGlobalXY(double geoX, double geoY, unsigned long long &globalX, unsigned long long &globalY);
	//void globalXYToGeo(unsigned long long globalX, unsigned long long globalY, double &geoX, double &geoY);
	// geoToGlobalXY_real(double geoX, double geoY, int &globalX, int &globalY);
	void geoToGlobalXY(double geoX, double geoY, int &globalX, int &globalY);
	void globalXYToGeo(long globalX, long globalY, double &geoX, double &geoY);

	void geotoLength(double dlon,double dlat, double lat, double *xyc);


	uint32_t getTotalX(){return totalX;}  // DGT made uint32_t, was long
	uint32_t getTotalY(){return totalY;}  // DGT made uint32_t, was long


	double getdyc(int index) {
		if (index < 0 || index >= totalY)
			return -1;

		return dyc[index];
	}

	double getdxc(int index) {
		if (index < 0 || index >= totalY)
			return -1;

		return dxc[index];
	}

	double getdxA() { return fabs(dxc[totalY/2]); }
	double getdyA() { return fabs(dyc[totalY/2]); }
	double getdlon() {return dlon;}
	double getdlat() {return dlat;}
	int getproj() {return IsGeographic;}


	double getXllcenter(){return xllcenter;}
	double getYllcenter(){return yllcenter;}

	OGRSpatialReferenceH getspatialref(){return hSRS;} // return projection information for raster file
	double getXLeftEdge(){return xleftedge;}
	double getYTopEdge(){return ytopedge;}
	DATA_TYPE getDatatype(){
		return datatype;
	}
	double getNodata() { 
		return nodata; 
	}  // noDatarefactor 11/18/17
	//void* getNodata(){return nodata;}
};
struct point{
	long x;
	long y;
	float elev;
	float area;
	float length;
};
struct streamlink{
	int32_t Id;
	int32_t u1;
	int32_t u2;
	long d;
	long magnitude;
	long shapeId;//from outletshapefile
	double elevU;
	double elevD;
	double length;
	short order;
	queue<point> coord; //  We store the coordinates of a link in a queue
	long numCoords;
	bool terminated;
};
struct llnode{
	streamlink *data;
	llnode *next;
};
struct LinkedLL{
	llnode *head;
	int numLinks;
};

LinkedLL linkSet;
void makeLinkSet(){
	linkSet.head = NULL;
	linkSet.numLinks = 0;
	return;
}

point* initPoint(CellSpace<double>& elev, CellSpace<double>& areaD8,CellSpace<double>& lengths,long iRow,long iCol);
void setLinkInfo(long **LinkIdU1U2DMagShapeid,double **LinkElevUElevDLength,double **PointXY, float **PointElevArea, double cellarea, tiffIO *elevIO);
void getNumLinksAndPoints(long &myNumLinks,long &myNumPoints);
void SendAndReciveLinks(int nx,int ny,CellSpace<double>&idGrid, CellSpace<double>&contribs,CellSpace<double>& flowDir, CellSpace<double>& src);
long getMagnitude(int32_t Id);
void appendPoint(int32_t Id, point* addPoint);
void setDownLinkId(int32_t Id, long dId);
streamlink* createLink(int32_t u1, int32_t u2, long d, point* coord); //,long numCoords);  //, float dx, float dy);
void linkSetInsert(streamlink* linkToAdd);
long GetOrderOf(long ID);
long findLinkThatEndsAt(long x, long y, CellSpace<double>&);
bool recvLink(int src);
bool sendLink(int32_t Id, int dest);
void terminateLink(int32_t Id);
//long findLinkThatStartsAt(long x, long y);
streamlink* FindLink(int32_t Id);
streamlink* getFirstLink();
streamlink* takeOut(int32_t Id);
llnode* LSInsert(llnode *head, llnode *addNode);

llnode* LSInsert(llnode *head, llnode *addNode){
	addNode->next = head;
	head = addNode;
	return head;
}
streamlink* takeOut(int32_t Id){
	streamlink* linkToBeRemoved;
	llnode* nodeToBeRemoved;
	if(linkSet.head == NULL){
		return NULL;//linkSet is empty
	}
	if(linkSet.head->data->Id == Id){
		//remove llnode and return streamlink* to head data	
		linkToBeRemoved = linkSet.head->data;
		nodeToBeRemoved = linkSet.head;
		linkSet.head = linkSet.head->next;
		delete nodeToBeRemoved;
		linkSet.numLinks--;
		return linkToBeRemoved;
	}else{// not the first llnode in the list
		llnode *previous, *current;
		current = linkSet.head;
		previous = NULL;
		while(current->data->Id != Id && current->next != NULL){//search the linked list
			previous = current;
			current = current->next;
		}
		if(current->data->Id == Id){//found it
			previous->next = current->next;
			nodeToBeRemoved = current;
			linkToBeRemoved = current->data;
			delete nodeToBeRemoved;
			linkSet.numLinks--;
			return linkToBeRemoved;
		}
		//We didn't find the ID
	}
	return NULL;
}
streamlink* FindLink(int32_t Id){
	streamlink* linkToBeRemoved;
	llnode* nodeToBeRemoved;
	if(linkSet.head == NULL){
		return NULL;//linkSet is empty
	}
	if(linkSet.head->data->Id == Id){
		return linkSet.head->data;
	}else{// not the first llnode in the list
		llnode* previous, *current;
		current = linkSet.head;
		previous = NULL;
		while(current->data->Id != Id && current->next != NULL){//search the linked list
			previous = current;
			current = current->next;
		}
		if(current->data->Id == Id){//found it
			return current->data;
		}
		//We didn't find the ID
	}
	return NULL;
}
void terminateLink(int32_t Id){
	streamlink *linkToKill;
	linkToKill = FindLink(Id);
	linkToKill = linkToKill;
	if(linkToKill == NULL){
		MPI_Abort(MPI_COMM_WORLD,4);
	}
	else{
		linkToKill->terminated = true;
	}
	return;
}
bool sendLink(int32_t Id, int dest){
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);


	streamlink *toSend = takeOut(Id);
	if(toSend == NULL)
		return false;
	if(toSend->terminated == true){
		linkSetInsert(toSend);
		return false;
	}

	MPI_Request req;
	MPI_Send(&(toSend->Id)		,1,MPI_LONG		,dest,1,MPI_COMM_WORLD);//send Id
	MPI_Send(&(toSend->u1)		,1,MPI_LONG		,dest,2,MPI_COMM_WORLD);
	MPI_Send(&(toSend->u2)		,1,MPI_LONG		,dest,3,MPI_COMM_WORLD);
	MPI_Send(&(toSend->d)		,1,MPI_LONG		,dest,4,MPI_COMM_WORLD);
	MPI_Send(&(toSend->elevU)	,1,MPI_DOUBLE	,dest,5,MPI_COMM_WORLD);
	MPI_Send(&(toSend->elevD)	,1,MPI_DOUBLE	,dest,6,MPI_COMM_WORLD);
	MPI_Send(&(toSend->length)	,1,MPI_DOUBLE		,dest,7,MPI_COMM_WORLD);
	MPI_Send(&(toSend->order)	,1,MPI_SHORT		,dest,8,MPI_COMM_WORLD);
	MPI_Send(&(toSend->numCoords)	,1,MPI_LONG	,dest,9,MPI_COMM_WORLD);
	MPI_Send(&(toSend->magnitude)	,1,MPI_LONG	,dest,11,MPI_COMM_WORLD);
	MPI_Send(&(toSend->shapeId)	,1,MPI_LONG	,dest,12,MPI_COMM_WORLD);

	//More information on This at :https://computing.llnl.gov/tutorials/mpi/
	//create mpi Data types
	MPI_Datatype PointType, oldtypes[2];  
	int          blockcounts[2]; 

	MPI_Aint offsets[2], extent;
	MPI_Status stat; 
	//set up first blocks of storage
	offsets[0] = 0;
	oldtypes[0] = MPI_LONG;
	blockcounts[0]= 2;
	//set up second block of storage
	MPI_Type_extent(MPI_LONG, &extent);
	offsets[1] = 2 * extent;
	oldtypes[1] = MPI_FLOAT;
	blockcounts[1] = 3;
	//create define it as an MPI data type and comit it.
	MPI_Type_struct(2,blockcounts,offsets,oldtypes,&PointType);
	MPI_Type_commit(&PointType);

	MPI_Status status;

	char *ptr;
	int place;
	point *buf;
	point *sendArr;
	sendArr = new point[toSend->numCoords];
	for(int i = 0; i < toSend->numCoords ;i++)
	{
		sendArr[i].x = toSend->coord.front().x;//y elec area length
		sendArr[i].y = toSend->coord.front().y;
		sendArr[i].elev = toSend->coord.front().elev;
		sendArr[i].area = toSend->coord.front().area;
		sendArr[i].length = toSend->coord.front().length;
		toSend->coord.pop();
	}	


	int bsize = toSend->numCoords*sizeof(point)*2+MPI_BSEND_OVERHEAD;  // Experimentally this seems to need to have >47 added for overhead
	buf = new point[bsize];

	MPI_Buffer_attach(buf,bsize);
	MPI_Bsend(sendArr,toSend->numCoords,PointType,dest,10,MPI_COMM_WORLD);
	MPI_Buffer_detach(&ptr,&place);
	delete sendArr;
	delete toSend;

	return true;
}
bool recvLink(int src){
	streamlink* toRecv = new streamlink;
	MPI_Status stat;
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	MPI_Recv(&(toRecv->Id)		,1,MPI_LONG		,src,1,MPI_COMM_WORLD,&stat);//send Id
	MPI_Recv(&(toRecv->u1)		,1,MPI_LONG		,src,2,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->u2)		,1,MPI_LONG		,src,3,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->d)		,1,MPI_LONG		,src,4,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->elevU)	,1,MPI_DOUBLE		,src,5,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->elevD)	,1,MPI_DOUBLE		,src,6,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->length)	,1,MPI_DOUBLE		,src,7,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->order)	,1,MPI_SHORT		,src,8,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->numCoords)	,1,MPI_LONG		,src,9,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->magnitude)	,1,MPI_LONG		,src,11,MPI_COMM_WORLD,&stat);
	MPI_Recv(&(toRecv->shapeId)	,1,MPI_LONG		,src,12,MPI_COMM_WORLD,&stat);

	//More information on This at 
	//create mpi Data types
	MPI_Datatype PointType, oldtypes[2];  
	int          blockcounts[2]; 

	MPI_Aint offsets[2], extent;
	MPI_Status stat1; 
	//set up first blocks of storage
	offsets[0] = 0;
	oldtypes[0] = MPI_LONG;
	blockcounts[0]= 2;
	//set up second block of storage
	MPI_Type_extent(MPI_LONG, &extent);
	offsets[1] = 2 * extent;
	oldtypes[1] = MPI_FLOAT;
	blockcounts[1] = 3;
	//create define it as an MPI data type and comit it.
	MPI_Type_struct(2,blockcounts,offsets,oldtypes,&PointType);
	MPI_Type_commit(&PointType);
	int flag;
	//MPI_Request req;
	//MPI_Test(&req,&flag,&stat1);
	//if(!(stat1.MPI_SOURCE >=0 && stat1.MPI_SOURCE < size))
	//	return false;
	//	toRecv->coord = new point[toRecv->numCoords];

	point *recvArr;
	recvArr = new point[toRecv->numCoords];
	MPI_Recv(recvArr,toRecv->numCoords,PointType,src,10,MPI_COMM_WORLD,&stat);
	toRecv->terminated = false;

	point temp;
	for(int i=0;i<toRecv->numCoords > 0;i++)
	{
		temp.x = recvArr[i].x;//y elec area length
		temp.y = recvArr[i].y;//y elec area length
		temp.elev = recvArr[i].elev;//y elec area length
		temp.area = recvArr[i].area;//y elec area length
		temp.length = recvArr[i].length;//y elec area length
		toRecv->coord.push(temp);
	}	

	linkSetInsert(toRecv);
	delete recvArr;

	return true;
}
void linkSetInsert(streamlink* linkToAdd)
{
	llnode *newNode = new llnode;
	newNode->data = linkToAdd;
	newNode->next = NULL;

	linkSet.numLinks++;
	linkSet.head = LSInsert(linkSet.head, newNode); 	
	//cout << "New ID entered: " << linkToAdd->Id << endl;
	return;
}

point* initPoint(CellSpace<double>& elev, CellSpace<double>& areaD8,CellSpace<double>& lengths,long iRow,long iCol, long gi, long gj)
{
	point *thePoint;
	thePoint=new point[1];
	//int gi,gj;
	float tempFloat;
	//gi=LoctoGloRow(iRow);
	//localToGlobal((int)i,(int)j,gi,gj);  
	thePoint->x = gi;
	thePoint->y = gj;//ÓÐÊ²Ã´ÓÃ£¿
	thePoint->elev = elev[iRow][iCol];//->getData(i,j,tempFloat);
	thePoint->area = areaD8[iRow][iCol];//->getData(i,j,tempFloat);
	thePoint->length = lengths[iRow][iCol];//->getData(i,j,tempFloat);
	return(thePoint);
}

streamlink* createLink(int32_t u1, int32_t u2, long d, point* coord) //,long numCoords) //, double dx, double dy)
{
	int size ;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	int rank ;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	streamlink *newLink = new streamlink;

	if(LAST_ID == -1){
		newLink->Id = rank;
	}else{
		newLink->Id =LAST_ID + size;
	}
	LAST_ID = newLink->Id;

	newLink->u1 = u1;	//maybe NULL
	newLink->u2 = u2;	//maybe NULL
	newLink->d = d;  	//maybe NULL

	newLink->coord.push(*coord);

	//newLink->coord = new point[numCoords];
	//memcpy(newLink->coord,coord,numCoords*sizeof(point));

	newLink->numCoords = 1; //numCoords;
	newLink->elevU = coord->elev;  //newLink->coord[0].elev;
	newLink->elevD = coord->elev;  //newLink->coord[numCoords-1].elev;
	newLink->order = 1;//a link starts at order 1
	newLink->terminated = false;

	newLink->magnitude = 1;
	int i=0;
	newLink->length = 0;
	//  DGT commented code below for efficiency reasons - length not used
	//float diag = sqrt(dx*dx+dy*dy);
	//for(i = 0; i < numCoords-1; i++){
	//	if(newLink->coord[i].x != newLink->coord[i+1].x && newLink->coord[i].y != newLink->coord[i+1].y)
	//		newLink->length += diag;
	//	else if(newLink->coord[i].x == newLink->coord[i+1].x )
	//		newLink->length += dy;
	//	else
	//		newLink->length += dx;
	//}
	newLink->shapeId = -1;
	//cout << rank << " Putting link into link set..."<< newLink->Id << endl;
	linkSetInsert(newLink);
	//cout << rank << " link put in linkset." << endl;
	return newLink;
}
void setDownLinkId(int32_t Id, long dId){
	streamlink* myLink = FindLink(Id);
	myLink->d = dId;
	return;
}
void appendPoint(int32_t Id, point* addPoint){
	streamlink* myLink = FindLink(Id);
	//cout << "ID = " << Id << endl;
	if(myLink == NULL)
		MPI_Abort(MPI_COMM_WORLD,8);//link not found cannot proccess
	//point* points;
	myLink->numCoords = myLink->numCoords +1;//make room

	//points = new point[myLink->numCoords];
	//	if(points == NULL)
	//		MPI_Abort(MCW,48);//link not found cannot proccess
	//cout << "before memcopy " << Id << endl;
	long i;
	//for(i=0; i < (myLink->numCoords) -1; i++)
	//{
	//	points[i]=myLink->coord[i];
	//}
	//memcpy(points,myLink->coord,(myLink->numCoords-1)*sizeof(point));
	//cout << "after memcopy " << Id << endl;
	//points[myLink->numCoords-1] = *addPoint;
	//delete [] myLink->coord;
	//cout << "after delete " << Id << endl;
	myLink->coord.push(*addPoint);  
	myLink->elevD = addPoint->elev;
	//myLink->numCoords++;
	//myLink->length++;  // DGT this is not used 
	return;
}

long getMagnitude(int32_t Id){
	streamlink* myLink = FindLink(Id);
	if(myLink == NULL)
		MPI_Abort(MPI_COMM_WORLD,7);
	return myLink->magnitude;
}
void getNumLinksAndPoints(long &NumLinks,long &NumPoints){
long ID = -1;
	llnode* current = linkSet.head;
	NumLinks = linkSet.numLinks;
	NumPoints = 0;
	if(linkSet.numLinks == 0)
		return ;
	else
		while(current != NULL){
			NumPoints += current->data->numCoords;
			current = current->next;
		}	
	return;
}
void setLinkInfo(long **LinkIdU1U2DMagShapeid,double **LinkElevUElevDLength,double **PointXY, float **PointElevArea, double cellarea, tiffIO *elevIO)
{
	long counter = 0;
	long pointsSoFar = 0;
	llnode* current = linkSet.head;
	if(linkSet.numLinks == 0)
		return;
	else
	{
        long begcoord=0;
		while(current != NULL){
			LinkIdU1U2DMagShapeid[counter][0] = current->data->Id;
			LinkIdU1U2DMagShapeid[counter][1] = begcoord;
			LinkIdU1U2DMagShapeid[counter][2] = begcoord+current->data->numCoords-1;
			begcoord=LinkIdU1U2DMagShapeid[counter][2]+1;
			LinkIdU1U2DMagShapeid[counter][3] = current->data->d;
			LinkIdU1U2DMagShapeid[counter][4] = current->data->u1;
			LinkIdU1U2DMagShapeid[counter][5] = current->data->u2;
			LinkIdU1U2DMagShapeid[counter][6] = current->data->order;
			LinkIdU1U2DMagShapeid[counter][7] = current->data->shapeId;
			LinkIdU1U2DMagShapeid[counter][8] = current->data->magnitude;
/*
			LinkElevUElevDLength[counter][0] = current->data->elevU;
			LinkElevUElevDLength[counter][1] = current->data->elevD;
			LinkElevUElevDLength[counter][2] = current->data->length;
	*/		
			long i=0;
			//double cellarea = (elev->getdxA())*(elev->getdyA());
			for(i=0;i<current->data->numCoords;i++){
				//elevIO->globalXYToGeo(current->data->coord[i].x, current->data->coord[i].y,PointXY[pointsSoFar][0],PointXY[pointsSoFar][1]);
				elevIO->globalXYToGeo(current->data->coord.front().x, current->data->coord.front().y,PointXY[pointsSoFar][0],PointXY[pointsSoFar][1]);	
				//PointElevArea[pointsSoFar][0] = current->data->coord[i].length;
				PointElevArea[pointsSoFar][0] = current->data->coord.front().length;
				//PointElevArea[pointsSoFar][1] = current->data->coord[i].elev;
				PointElevArea[pointsSoFar][1] = current->data->coord.front().elev;
				//PointElevArea[pointsSoFar][2] = current->data->coord[i].area*cellarea;
				PointElevArea[pointsSoFar][2] = current->data->coord.front().area*cellarea;
				//  Pop off queue	
				current->data->coord.pop();
				pointsSoFar++;
			}
			current = current->next;
			counter++;
		}	
	}
	return;
}

const short TIFF = 42;
const short BIGTIFF = 43;

struct geotiff{
	long xresNum;              //??? unsigned long BT - numerator for horizontal resolution 
	long xresDen;              //??? unsigned long BT - denominator for horizontal resolution
	long yresNum;              //??? unsigned long BT - numerator for vertical resolution
	long yresDen;              //??? unsigned long BT - denominator for horizontal resolution
	short planarConfig;        //??? unsigned short BT - irrelevant when SamplesPerPixel=1

	uint32_t geoKeySize;      	//  DGT made uint32_t	//??? unsigned long long BT
	uint32_t geoDoubleSize;   	//  DGT made uint32_t	//??? unsigned long long BT
	uint32_t geoAsciiSize;    	//  DGT made uint32_t	//??? unsigned long long BT
	uint32_t gdalAsciiSize;   	//  DGT made uint32_t	//??? unsigned long long BT

	uint16_t *geoKeyDir;        //  DGT made uint16_t  	// unsigned short BT - pointer to extended geoKeyDir metadata
	double *geoDoubleParams;    //pointer to extended geoDoubleParams metadata
	char *geoAscii;             //pointer to extended geoAscii metadata
	char *gdalAscii;            //pointer to extended gdalAscii metadata

};

//Image File Directory (IFD) - TIFF Tag Structure
struct ifd {
	unsigned short tag;			//Tag ID#
	unsigned short type;		//Datatype of Values (TIFF datatypes, not C++)
	uint32_t count;	  //DGT	//Count of Values
	uint32_t offset;	//DGT	// unsigned long long BT - Values (if fits in 4 bytes for TIFF or 8 for BIGTIFF else Offset to Values)
};

#endif