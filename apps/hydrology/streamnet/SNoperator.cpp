#include"SNoperator.h"


OGRSFDriverH    driver;
OGRDataSourceH  hDS1;
OGRLayerH       hLayer1;
OGRFeatureDefnH hFDefn1;
OGRFieldDefnH   hFieldDefn1;
OGRFeatureH     hFeature1;
OGRGeometryH    geometry, line;
OGRSpatialReferenceH  hSRSraster,hSRSshapefile;
void SNOperator::
	demLayer(RasterLayer<double>& layerD) {
		_pDEMLayer = &layerD;
		_pDEMNbrhood = layerD.nbrhood();
		_cellSize = _pDEMLayer->_pMetaData->cellSize;
		_noData = _pDEMLayer->_pMetaData->noData;
		_maxRow = _pDEMLayer->_pMetaData->_localworkBR.maxIRow();
		_maxCol = _pDEMLayer->_pMetaData->_localworkBR.maxICol();
		_rank=_pDEMLayer->_pMetaData->myrank;
		_contribs.copyLayerInfo(*_pDEMLayer);
		_wsgrid.copyLayerInfo(*_pDEMLayer);
		_idgrid.copyLayerInfo(*_pDEMLayer);
		_lengths.copyLayerInfo(*_pDEMLayer);
		Configure(_pDEMLayer, false);
		Configure(&_contribs, true);
		Configure(&_wsgrid, true);
		Configure(&_idgrid, true);
		Configure(&_lengths, true);
}

void SNOperator::srcLayer(RasterLayer<double>& layerD) {
	_psrcLayer = &layerD;
	_pDEMNbrhood = layerD.nbrhood();
	Configure(_psrcLayer, true);
}
void SNOperator::dirLayer(RasterLayer<double>& layerD) {

	_pDirLayer = &layerD;
	_pDEMNbrhood = layerD.nbrhood();
	Configure(_pDirLayer, false);
}
void SNOperator::areaLayer(RasterLayer<double>& layerD) {

	_pareaLayer = &layerD;
	_pDEMNbrhood = layerD.nbrhood();
	Configure(_pareaLayer, false);
}


bool SNOperator::isTermination() {
	return true;
	//num--;
	//if(num > 0)
	//{
	//	return true;
	//}
	//else
	//{
	//	return false;
	//}
}
void tiffIO::geotoLength(double dlon,double dlat, double lat, double *xyc){
	double ds2,beta,dbeta;

	dlat=dlat*PI/180.;
	dlon=dlon*PI/180.;
	lat=lat*PI/180.;
	beta = atan(boa*tan(lat));
	dbeta=dlat*boa*(cos(beta)/cos(lat))*(cos(beta)/cos(lat));
	ds2=(pow(elipa*sin(beta),2)+pow(elipb*cos(beta),2))*pow(dbeta,2);
	xyc[0]=elipa*cos(beta)*abs(dlon);
	xyc[1]=double(sqrt(double(ds2)));    
}
void getLayername(char *inputogrfile, char *layername)
{  
	std::string filenamewithpath;
	filenamewithpath=inputogrfile;
	size_t found = filenamewithpath.find_last_of("/\\");
	std::string filenamewithoutpath;
	filenamewithoutpath=filenamewithpath.substr(found+1);
	const char *filename = filenamewithoutpath.c_str(); // convert string to char
	const char *ext; 
	ext = strrchr(filename, '.'); // getting extension
	//char layername[MAXLN];
	size_t len = strlen(filename);
	size_t len1 = strlen(ext);
	memcpy(layername, filename, len-len1);
	layername[len - len1] = 0; 
	printf("%s ", layername);
	return;
	//return layername;
}
const char *getOGRdrivername(char *datasrcnew){
	const char *ogrextension_list[5] = {".sqlite",".shp",".json",".kml",".geojson"};  // extension list --can add more 
	const char *ogrdriver_code[5] = {"SQLite","ESRI Shapefile","GeoJSON","KML","GeoJSON"};   //  code list -- can add more
	size_t extension_num=5;
	char *ext; 
	int index = 1; //default is ESRI shapefile
	ext = strrchr(datasrcnew, '.'); 
	if(!ext){

		index=1; //  if no extension then writing will be ESRI shapefile
	}
	else
	{

		//  convert to lower case for matching
		for(int i = 0; ext[i]; i++){
			ext[i] = tolower(ext[i]);
		}
		// if extension matches then set driver
		for (size_t i = 0; i < extension_num; i++) {
			if (strcmp(ext,ogrextension_list[i])==0) {
				index=i; //get the index where extension of the outputfile matches with the extensionlist 
				break;
			}
		}

	}

	const  char *drivername;
	drivername=ogrdriver_code[index];
	return drivername;
}
void createStreamNetShapefile(char *streamnetsrc,char *streamnetlyr,OGRSpatialReferenceH hSRSraster){

	/* Register all OGR drivers */
	OGRRegisterAll();
	//const char *pszDriverName = "ESRI Shapefile";
	const char *pszDriverName;
	pszDriverName=getOGRdrivername(streamnetsrc);

	//get driver by name
	driver = OGRGetDriverByName( pszDriverName );
	if( driver == NULL )
	{
		printf( "%s warning: driver not available.\n", pszDriverName );
		//exit( 1 );
	}

	// open datasource if the datasoruce exists 
	if(pszDriverName=="SQLite") hDS1 = OGROpen(streamnetsrc, TRUE, NULL );
	// create new data source if data source does not exist 
	if (hDS1 ==NULL){ 
		hDS1 = OGR_Dr_CreateDataSource(driver, streamnetsrc, NULL);}
	else { hDS1=hDS1 ;}

	if (hDS1 != NULL) {


		// layer name is file name without extension
		if(strlen(streamnetlyr)==0){
			// Chris George suggestion
			char streamnetlayername[MAXLN];
			getLayername(streamnetsrc, streamnetlayername); // get layer name if the layer name is not provided		  
			//char *streamnetlayername;
			//streamnetlayername=getLayername(streamnetsrc); // get layer name if the layer name is not provided
			hLayer1= OGR_DS_CreateLayer( hDS1,streamnetlayername,hSRSraster, wkbLineString, NULL );} 

		else {
			hLayer1= OGR_DS_CreateLayer( hDS1,streamnetlyr,hSRSraster, wkbLineString, NULL ); }// provide same spatial reference as raster in streamnetshp file
		if( hLayer1 == NULL )
		{
			printf( "warning: Layer creation failed.\n" );
			//exit( 1 );
		}



		/* Add a few fields to the layer defn */ //need some work for setfiled width
		// add fields, field width and precision
		hFieldDefn1 = OGR_Fld_Create( "LINKNO", OFTInteger ); // create new field definition
		OGR_Fld_SetWidth( hFieldDefn1, 6); // set field width
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0); // create field 

		hFieldDefn1 = OGR_Fld_Create( "DSLINKNO", OFTInteger );
		OGR_Fld_SetWidth( hFieldDefn1, 6);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "USLINKNO1", OFTInteger );
		OGR_Fld_SetWidth( hFieldDefn1, 6);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "USLINKNO2", OFTInteger );
		OGR_Fld_SetWidth( hFieldDefn1, 6);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "DSNODEID", OFTInteger );
		OGR_Fld_SetWidth( hFieldDefn1, 12);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "strmOrder", OFTInteger );
		OGR_Fld_SetWidth( hFieldDefn1, 6);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "Length", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1, 16);
		OGR_Fld_SetPrecision(hFieldDefn1, 1); // set precision
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "Magnitude", OFTInteger );
		OGR_Fld_SetWidth( hFieldDefn1, 6);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "DSContArea", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1, 16);
		OGR_Fld_SetPrecision(hFieldDefn1, 1);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "strmDrop", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1, 16);
		OGR_Fld_SetPrecision(hFieldDefn1, 2);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1= OGR_Fld_Create( "Slope", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1,16);
		OGR_Fld_SetPrecision(hFieldDefn1, 12);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "StraightL", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1,16);
		OGR_Fld_SetPrecision(hFieldDefn1, 1);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "USContArea", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1, 16);
		OGR_Fld_SetPrecision(hFieldDefn1, 1);
		OGR_L_CreateField(hLayer1,hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "WSNO", OFTInteger );
		OGR_Fld_SetWidth( hFieldDefn1, 6);
		OGR_Fld_SetPrecision(hFieldDefn1,0);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "DOUTEND", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1, 16);
		OGR_Fld_SetPrecision(hFieldDefn1, 1);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "DOUTSTART", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1, 16);
		OGR_Fld_SetPrecision(hFieldDefn1, 1);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);

		hFieldDefn1 = OGR_Fld_Create( "DOUTMID", OFTReal );
		OGR_Fld_SetWidth( hFieldDefn1, 16);
		OGR_Fld_SetPrecision(hFieldDefn1, 1);
		OGR_L_CreateField(hLayer1,  hFieldDefn1, 0);
	}
}

int reachshape(long *cnet,float *lengthd, float *elev, float *area, double *pointx, double *pointy, long np,tiffIO &obj)
{
	// Function to write stream network shapefile
	int nVertices;
	if (np < 2) {//singleton - will be duplicated
		nVertices = 2;
	}
	else {
		nVertices = np;
	}
	double *mypointx = new double[nVertices];
	double *mypointy = new double[nVertices];
	double x,y,length,glength,x1,y1,xlast,ylast,usarea,dsarea,dslast,dl,drop,slope;
	int istart,iend,j;

	istart=cnet[1];  //  start coord for first link
	iend=cnet[2];//  end coord for first link
	x1=pointx[0];
	y1=pointy[0];
	length=0.;
	xlast=x1;
	ylast=y1;
	usarea=area[0];
	dslast=usarea;
	dsarea=usarea;
	long prt = 0;

	for(j=0; j<np; j++)  //  loop over points
	{
		x=pointx[j];
		y=pointy[j];
		// we have to reverse order of pointx and pointy arrays to finish up with 
		// the point with point index 0 in the shape being the outlet point
		// (which is for backwards compatibility with previous versions of TauDEM)
		int i = np - (j + 1);
		mypointx[i] = x;
		mypointy[i] = y;
		if(obj.getproj()==1){
			double dlon1,dlat1;double dxdyc1[2];
			dlon1=x-xlast;dlat1=y-ylast;
			obj.geotoLength(dlon1,dlat1,y,dxdyc1);
			dl=sqrt(dxdyc1[0]*dxdyc1[0]+dxdyc1[1]*dxdyc1[1]);
		}
		else if(obj.getproj()==0){
			dl=sqrt((x-xlast)*(x-xlast)+(y-ylast)*(y-ylast));
		}

		if(dl > 0.)
		{
			length=length+dl;
			xlast=x;  ylast=y;
			dsarea=dslast;   // keeps track of last ds area
			dslast=area[j];
		}
	}
	drop=elev[0]-elev[np-1];
	slope=0.;
	float dsdist=lengthd[np-1];
	float usdist=lengthd[0];
	float middist=(dsdist+usdist)*0.5;
	if(length > 0.)slope=drop/length;

	if(obj.getproj()==1){
		double dlon2,dlat2;double dxdyc2[2];
		dlon2=x-x1;dlat2=y-y1;
		obj.geotoLength(dlon2,dlat2,y,dxdyc2);
		glength=sqrt(dxdyc2[0]*dxdyc2[0]+dxdyc2[1]*dxdyc2[1]);

	}
	else if(obj.getproj()==0){
		glength=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
	}

	// ensure at least two points (assuming have at least 1) by repeating singleton
	if (np < 2) {
		mypointx[1] = mypointx[0];
		mypointy[1] = mypointy[0];
	}

	//SHPObject *shape = SHPCreateSimpleObject(
	//	SHPT_ARC,						// type
	//	nVertices,						// number of vertices
	//	mypointx,						// X values
	//	mypointy,						// Y values
	//	NULL);							// Z values

	hFDefn1 = OGR_L_GetLayerDefn( hLayer1 );
	hFeature1 = OGR_F_Create( hFDefn1 );
	OGR_F_SetFieldInteger( hFeature1, 0, (int)cnet[0]); // set field value 
	OGR_F_SetFieldInteger( hFeature1, 1, (int)cnet[3]);
	OGR_F_SetFieldInteger( hFeature1, 2, (int)cnet[4]);
	OGR_F_SetFieldInteger( hFeature1, 3, (int)cnet[5]);
	OGR_F_SetFieldInteger( hFeature1, 4, (int)cnet[7]);
	OGR_F_SetFieldInteger( hFeature1, 5, (int)cnet[6]);
	OGR_F_SetFieldDouble( hFeature1, 6, length);
	OGR_F_SetFieldInteger( hFeature1, 7, (int)cnet[8]);
	OGR_F_SetFieldDouble( hFeature1, 8,  dsarea);
	OGR_F_SetFieldDouble( hFeature1, 9, drop);
	OGR_F_SetFieldDouble( hFeature1, 10, slope);
	OGR_F_SetFieldDouble( hFeature1, 11, glength);
	OGR_F_SetFieldDouble( hFeature1, 12, usarea);
	OGR_F_SetFieldInteger( hFeature1, 13, (int)cnet[0]);
	OGR_F_SetFieldDouble( hFeature1, 14, dsdist);
	OGR_F_SetFieldDouble( hFeature1, 15, usdist);
	OGR_F_SetFieldDouble( hFeature1, 16, middist);

	//creating geometry using OGR

	line = OGR_G_CreateGeometry( wkbLineString );
	for(j=0; j<np; j++) {
		OGR_G_AddPoint(line, mypointx[j], mypointy[j], 0);
	}
	OGR_F_SetGeometryDirectly(hFeature1,line); // set geometry to feature
	OGR_L_CreateFeature( hLayer1, hFeature1 ); //adding feature 

	delete[] mypointx;
	delete[] mypointy;

	return 0;
}
struct Slink{
	long id;
	short dest;
};



bool SNOperator::Operator(const CellCoord& coord, bool operFlag) {
	CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
	CellSpace<double>& src = *(_psrcLayer->cellSpace());
	CellSpace<double>& dir = *(_pDirLayer->cellSpace());
	CellSpace<double>& area = *(_pareaLayer->cellSpace());

	CellSpace<double>& lengths = *(_lengths.cellSpace());
	CellSpace<double>& contribs = *(_contribs.cellSpace());
	CellSpace<double>& wsGrid = *(_wsgrid.cellSpace());
	CellSpace<double>& idGrid = *(_idgrid.cellSpace());
	tiffIO elevIO(elevfile,elevtype);
	tiffIO srcIO(srcfile,srctype);
	//CellSpace<double>& d8L = *(_pD8Layer->cellSpace());
	//CellSpace<double>& scaL = *(_pSCALayer->cellSpace());
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	int minRow = _pDEMLayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pDEMLayer->_pMetaData->_localworkBR.minICol();

	queue <Slink> linkQ;
	long iRow = coord.iRow();
	long iCol = coord.iCol();
	long gi=_pDEMLayer->_pMetaData->_MBR.minIRow()+iRow;
	long gj=_pDEMLayer->_pMetaData->_MBR.minIRow()+iCol;
	//long Sup ;//= 0;
	//long Sdown;// = 0;
	short d;
	Slink temp;
	//node t;
	int p;
	int nexti,nextj;
	int d1[9]={-1,-1,-1,0,0,0,1,1,1};
	int d2[9]={-1,0,1,-1,0,1,-1,0,1};

	double lenx=iNeighborCells*_cellSize;
	double leny=iNeighborCells*_cellSize;
	long tRow,tCol;
	//int k;
	switch(num)
	{
	case 0:
	case 1:{


		if(num==0){
			cout<<"hehe"<<endl;
			//int d = 8;
			//if(iCol%25==0){
			//	cout<<"cc:"<<contribs[iRow][iCol];
			//	cout<<"row:"<<iRow<<"rank:"<<_rank<<endl;
			//}
			contribs[iRow][iCol]=-2;

			if(fabs(src[iRow][iCol] - _noData) > Eps&&src[iRow][iCol]==1&&!fabs(dir[iRow][iCol]-_noData)<Eps&&dir[iRow][iCol]>=0){	
				d=dir[iRow][iCol];
				tRow=iRow+d1[d];
				tCol=iCol+d2[d];

				if(fabs(src[tRow][tCol] - _noData) > Eps&&src[tRow][tCol]==1&&!fabs(dir[tRow][tCol]-_noData)<Eps&&dir[tRow][tCol]>=0){
					contribs[iRow][iCol]=1;
				}else
					contribs[iRow][iCol]=0;
			}

			if (iRow == _maxRow && iCol == _maxCol) {
				//cout<<"cccc:"<<contribs[iRow][iCol];
				MPI_Barrier(MPI_COMM_WORLD);
				num = 1;
				Termination = 0;

			}
			return true;
		}


		//if(cas==0){


		if (contribs[iRow][iCol] < 0 && !(iRow == _maxRow && iCol == _maxCol)) {
			return true;
		}
		if(num==1){
			if(!contribs[iRow][iCol]){

				d=dir[iRow][iCol];
				float llength=0.;
				tRow=iRow+d1[d];
				tCol=iCol+d2[d];
				if(!fabs(lengths[tRow][tCol]-_noData)<Eps)
				{
					llength += lengths[tRow][tCol];

					if(d==1||d==7)llength=llength+leny; // make sure that it is correct
					if(d==3||d==5)llength=llength+lenx;
					if(d%2==0)llength=llength+sqrt(lenx*lenx+leny*leny);
				}
				lengths[iRow][iCol]=llength;
				//totalP++;
				//  Find if neighbor points to me and reduce its dependency by 1


				contribs[iRow][iCol]=-1;
			}//else if(contribs[iRow][iCol] ==1){
				//d=dir[iRow][iCol];
				//tRow=d1[d]+iRow;
				//tCol=d2[d]+iCol;
				//if(contribs[tRow][tCol]==-1)
				//	contribs[iRow][iCol] =0;


			//}
			if (iRow == _maxRow && iCol == _maxCol) {
				cout<<"hehe"<<num<<endl;
				MPI_Barrier(MPI_COMM_WORLD);
				bool flag=true;
				num=2;
				for (int i = minRow; i <= _maxRow; ++i) {
					for (int j = minCol; j <= _maxCol; ++j) {
						if (contribs[i][j] >=0) {
							d=dir[i][j];
							tRow=d1[d]+i;
							tCol=d2[d]+j;
							if(contribs[tRow][tCol]==-1)
								contribs[i][j] =0;

							//cout<<"row:"<<i<<" col:"<<j<<" con:"<<contribs[i][j]<<endl;
							
							num=1;
							
							//num++;
							//cout<<num;
							//return true;
						}
						Termination = 0;
					}
				}
				MPI_Barrier(MPI_COMM_WORLD);
				cout<<"num:"<<num<<endl;
				//if(flag)//termination should not be changed in if,else
				//{
				//	num=2;
				//	Termination = 0;
				//}
			}
			return true;
		}
		   }
		   break;
	case 2:
		{
			MPI_Barrier(MPI_COMM_WORLD);
			cout<<"num:"<<num<<endl;
			num=5;
			if(num==7)
			{	
				contribs[iRow][iCol]=-2;

				if(fabs(src[iRow][iCol] - _noData) > Eps&&src[iRow][iCol]==1&&!fabs(dir[iRow][iCol]-_noData)<Eps&&dir[iRow][iCol]>=0){
					contribs[iRow][iCol]=0;
					for (tRow = iRow - 1; tRow <= iRow + 1; tRow++) {
						for (tCol = iCol - 1; tCol <= iCol + 1; tCol++) {
							if(!fabs(dir[tRow][tCol]-_noData)<Eps&&dir[tRow][tCol]>=0){

								d=dir[tRow][tCol];

								if (d1[d]+tRow==iRow && d2[d]+tCol==iCol&&d + dir[iRow][iCol] != 8&& d != 4 ) {
									//) 
									if(src[tRow][tCol]>01&&fabs(src[iRow][iCol] - _noData) > Eps)contribs[iRow][iCol]++;

								}
							}
						}
					}
				}
				Sup = 0;
				Sdown = 0;	
				if (iRow == _maxRow && iCol == _maxCol) {
					//cout<<"cccc:"<<contribs[iRow][iCol];
					MPI_Barrier(MPI_COMM_WORLD);
					num = 3;
					Termination = 0;

				}
				return true;
			}
			if (contribs[iRow][iCol] < 0 && !(iRow == _maxRow && iCol == _maxCol)) {
				return true;
			}
			if(num==3){
				bool linksToReceive=true;
				if(!contribs[iRow][iCol]){

					short nOrder[8];  // neighborOrders
					short inneighbors[8];  // inflowing neighbors
					bool junction; // junction set to true/false in newOrder
					int32_t Id;

					//  Count number of grid cells on stream that point to me.  
					// (There will be no more than 8, usually only 0, 1,2 or 3.  The queue and 
					// processing order ensures that the necessary results have already been 
					//  evaluated for these.  In this case it should be a partially complete link 
					//  object referenced by a value in a link object partition).  Number that flow 
					//  to cell is k.  Record the order and direction to each. 
					int k=0;
					//int d=8;
					//for(m=1;m<=8;++m){
					for (int tRow = iRow - 1; tRow <= iRow + 1; tRow++) {
						for (int tCol = iCol - 1; tCol <= iCol + 1; tCol++) {
							if(!fabs(dir[tRow][tCol]-_noData)<Eps&&dir[tRow][tCol]>=0){

								d=dir[tRow][tCol];

								if (d1[d]+tRow==iRow && d2[d]+tCol==iCol&&d + dir[iRow][iCol] != 8&& d != 4 ) {
									//) 
									if(src[tRow][tCol]>0&&fabs(src[tRow][tCol] - _noData) > Eps&&fabs(idGrid[tRow][tCol] - _noData) > Eps){
										inneighbors[k]=8-d;  
										nOrder[k]=(short)lengths[tRow][tCol];
										k++;
									}

								}
							}
						}
					}




					p=dir[iRow][iCol];
					nexti=iRow+d1[p];
					nextj=iCol+d2[p];
					// Determine downslope neighbor
					//long nextx,nexty;
					//p=flowDir->getData(i,j,tempShort);
					//nextx=t.x+d1[p];
					//nexty=t.y+d2[p];

					// Determine if mandatory junction indicated by a value in idGrid
					bool mandatoryJunction=false;
					int32_t shapeID;
					if(!fabs(idGrid[iRow][iCol] - _noData) < Eps)
					{
						mandatoryJunction=true;
						shapeID=idGrid[iRow][iCol];
					}

					//-	determine if terminal link by if downslope neighbor is out of domain
					bool terminal = false;
					if(fabs(src[nexti][nextj] - _noData) > Eps 
						|| src[nexti][nextj] != 1)
						terminal = true;

					//  instantiate point
					point* addPoint;
					addPoint = initPoint(dem,area,lengths,iRow,iCol,gi,gj);

					//  Case for start of flow path
					if(k==0){
						//  This is the start of a stream - instantiate new link
						streamlink *thisLink;
						thisLink = createLink(-1, -1, -1, addPoint); //,1); //, dx,dy);
						int32_t u1= thisLink->Id;
						idGrid[iRow][iCol]=u1;//->setData(i,j, u1);
						wsGrid[iRow][iCol]=u1;//->setData(i,j, u1);
						lengths[iRow][iCol]=(float)thisLink->order;//->setData(i,j,(float)thisLink->order);

						//  Handle mandatory junction
						if(mandatoryJunction){
							appendPoint(u1, addPoint);
							thisLink->shapeId = shapeID;
							if(!terminal){
								streamlink *newLink;
								newLink = createLink(u1,-1,-1, addPoint); //, 1); //, dx,dy);
								newLink->order =  1;	
								setDownLinkId(u1,newLink->Id);		
								newLink->magnitude = 1;
								terminateLink(u1);
								idGrid[iRow][iCol]=newLink->Id;
								//  Here idGrid and wsGrid intentionally out of sync
								//  idGrid to be used to retrieve upstream link for next cell down
								//  wsGrid used to hold the id for delineation of watersheds which 
								//  is the link that ends at this point
							}
						}else if(terminal){  // handle terminal that is not mandatory junction
							appendPoint(u1, addPoint);
						}
					}

					//  Case for ongoing flow path with single inflow
					else if(k==1){
						int32_t u1=idGrid[iRow+d1[inneighbors[0]]][iCol+d2[inneighbors[0]]];
						appendPoint(u1, addPoint);
						streamlink *thisLink;
						//TODO DGT Thinks
						// thisLink=u1;  for efficiency
						thisLink=FindLink(u1);
						//cout << "append done"  << endl;
						idGrid[iRow][iCol]=u1;//->setData(i,j,u1);
						wsGrid[iRow][iCol]=u1;//->setData(i,j,u1);
						lengths[iRow][iCol]=(float)nOrder[0];

						//  Handle mandatory junction case
						if(mandatoryJunction){
							thisLink->shapeId = shapeID;
							if(!terminal){
								streamlink *newLink;
								newLink = createLink(u1,-1,-1, addPoint); //, 1); //, dx,dy);
								newLink->order =  thisLink->order;	
								setDownLinkId(u1,newLink->Id);		
								newLink->magnitude = thisLink->magnitude;
								terminateLink(u1);
								idGrid[iRow][iCol]=newLink->Id;
								//  TODO Here need to assign WSNO = -1
							}
						}
					}

					// Case for multiple inflows
					else if(k>1){
						//	rank incoming flow paths from highest to lowest order and process in this order
						//cout << "k>1 junction: "  << k << endl;
						//Dans Code to sort the two largest to the bottom of the list
						short temp;
						for(int a=0;a<2;++a){
							for(int b=0;b<k-1;++b){
								if(nOrder[b]>nOrder[b+1]){
									temp=nOrder[b];
									nOrder[b]=nOrder[b+1];
									nOrder[b+1]=temp;

									temp = inneighbors[b];
									inneighbors[b] = inneighbors[b+1];
									inneighbors[b+1] = temp;//keep track of where they came from. 
								}
							}
						}	
						// Determine order out
						short oOut;//if the last two have the same order, there is a junction, so bump the order
						if(nOrder[k-2]==nOrder[k-1]){
							//junction=true;
							oOut=nOrder[k-1]+1;
						}else{
							oOut=nOrder[k-1];
						}
						int32_t u1, u2;
						u1 = idGrid[iRow+d1[inneighbors[k-1]]][iCol+d2[inneighbors[k-1]]];
						u2 = idGrid[iRow+d1[inneighbors[k-2]]][iCol+d2[inneighbors[k-2]]];
						appendPoint(u1, addPoint);
						appendPoint(u2, addPoint);
						streamlink *newLink;
						newLink = createLink(u1,u2,-1,addPoint); //,1); //,dx,dy);
						newLink->order = oOut;
						setDownLinkId(u1,newLink->Id);
						setDownLinkId(u2,newLink->Id);
						newLink->magnitude = getMagnitude(u1) + getMagnitude(u2);

						//  TODO add WSNO = -1

						terminateLink(u1);
						terminateLink(u2);

						//  Handle remaining incoming flow paths
						for(int m=2; m<k; m++)   //  Loop is only entered if k>=3
						{
							u1 = newLink->Id;
							u2 = idGrid[iRow+d1[inneighbors[k-m-1]]][iCol+d2[inneighbors[k-m-1]]];

							appendPoint(u1, addPoint);
							appendPoint(u2, addPoint);

							// create a new link
							newLink = createLink(u1,u2,-1,addPoint); //,1);  //,dx,dy);
							//  set order of new link as order out
							newLink->order = oOut;
							//  set u1 to previously created link
							setDownLinkId(u1, newLink->Id);
							setDownLinkId(u2, newLink->Id);
							newLink->magnitude = getMagnitude(u1) + getMagnitude(u2);

							//  TODO add WSNO = -1

							terminateLink(u1);
							terminateLink(u2);
						}
						idGrid[iRow][iCol]=newLink->Id;// Set the idGrid to the new link ID?
						wsGrid[iRow][iCol]=newLink->Id;
						lengths[iRow][iCol]=(float)oOut;
						if(mandatoryJunction){
							// create a zero length link to the mandatory junction
							u1 = newLink->Id;
							appendPoint(u1, addPoint);
							newLink->shapeId=shapeID;

							//  TODO assign u1 to WSNO
							if(!terminal){
								newLink = createLink(u1,-1,-1,addPoint); //,1); //,dx,dy);
								//  set order of new link as order out
								newLink->order = oOut;
								//  set u1 to previously created link
								setDownLinkId(u1, newLink->Id);
								newLink->magnitude = getMagnitude(u1);
								idGrid[iRow][iCol]=newLink->Id;// Set the idGrid to the new link ID?
								//  TODO assign -1 to WNSO
							}
						}else if(terminal){
							u1 = newLink->Id;
							appendPoint(u1, addPoint);
							//  TODO assign u1 to WSNO
						}
					}
					contribs[iRow][iCol]=-1;
					d=dir[iRow][iCol];
					nexti=iRow+d1[d];
					nextj=iCol+d2[d];
					//if(x>=0 && x<nx && y>=0 && y<ny)
					if(nexti<1||nexti>_maxRow){
						if(nexti<1&&_rank>0){
							temp.id = idGrid[iRow][iCol];//,j,tempLong);
							temp.dest = _rank-1;
							linkQ.push(temp);
							Sup++;
						}
						if(nexti>_maxRow&&_rank<size-1)
						{
							temp.id = idGrid[iRow][iCol];//,j,tempLong);
							temp.dest = _rank+1;
							linkQ.push(temp);
							Sdown++;
						}
					}


				}else if(contribs[iRow][iCol]>0){
					//cout<<"undo r:"<<iRow<<"c:"<<iCol<<"cc:"<<contribs[iRow][iCol]<<endl;
					bool g=true;
					for(tRow = iRow - iNeighborCells; tRow <= iRow + iNeighborCells; tRow++) {
						for (tCol = iCol - iNeighborCells; tCol <= iCol + iNeighborCells; tCol++){
							if(dir[tRow][tCol]>=0){
								d=dir[tRow][tCol];
								if(d1[d]+tRow ==iRow&& d2[d]+tCol==iCol &&d + dir[iRow][iCol] != 8&& d != 4 &&contribs[tRow][tCol]>-1)
									g=false;	
							}
						}
						if(g){
							contribs[iRow][iCol]=0;
						}
					}
				}
				if (iRow == _maxRow && iCol == _maxCol) {

					MPI_Barrier(MPI_COMM_WORLD);
					bool flag=true;
					for (int i = minRow; i <= _maxRow; ++i) {
						for (int j = minCol; j <= _maxCol; ++j) {
							if (contribs[i][j] >=0) {
								//cout<<"row:"<<i<<" col:"<<j<<" con:"<<contribs[i][j]<<endl;
								Termination = 0;
								///num++;
								//g=false;
								//cout<<num;
								//return true;
							}
						}
					}
					if(size > 1)
					{
						long numRecv;
						MPI_Status stat;
						if(rank < size-1){//only send down
							MPI_Send(&Sdown,1,MPI_LONG,rank+1,0,MPI_COMM_WORLD);//send message saying done
							for(int j =0; j<Sdown; j++){
								temp.id = linkQ.front().id;
								temp.dest = linkQ.front().dest;
								linkQ.pop();//all an the q should go to the same dest.
								while(temp.dest != rank+1){
									linkQ.push(temp);
									temp.id = linkQ.front().id;
									temp.dest = linkQ.front().dest;
									linkQ.pop();//all an the q should go to the same dest.
								}
								sendLink(temp.id,temp.dest);
							}
						}
						if(rank >0)
						{
							MPI_Recv(&numRecv,1,MPI_LONG,rank-1,0,MPI_COMM_WORLD,&stat);
							for(int j =0 ;j<numRecv;j++){
								recvLink(rank-1);//recv from up.
							}
						}
						if(rank > 0)
						{
							MPI_Send(&Sup,1,MPI_LONG,rank-1,0,MPI_COMM_WORLD);//send message saying up
							for(int j =0; j<Sup; j++){
								temp.id = linkQ.front().id;
								temp.dest = linkQ.front().dest;
								linkQ.pop();//all an the q should go to the same dest.
								sendLink(temp.id,temp.dest);
							}
						}
						if(rank < size -1)
						{
							MPI_Recv(&numRecv,1,MPI_LONG,rank+1,0,MPI_COMM_WORLD,&stat);
							for(int j = 0;j < numRecv; j++){
								recvLink(rank+1);//recv from below.
							}
						}

						//MPI_Barrier(MCW);
						Sup = 0;
						Sdown = 0;
					}

				}

				return true;
			}
		}break;

	case 4:
		{

			if(num==4)
			{
				contribs[iRow][iCol]=-2;
				if(!fabs(dir[iRow][iCol]-_noData)<Eps&&dir[iRow][iCol]>=0){
					if(fabs(src[iRow][iCol] - _noData) < Eps||src[iRow][iCol]<1){
						d=dir[iRow][iCol];
						tRow=iRow+d1[d];
						tCol=iCol+d2[d];
						if(fabs(src[tRow][tCol] - _noData) < Eps||src[tRow][tCol]<1){
							contribs[iRow][iCol]=(short)1;
							wsGrid[iRow][iCol]=_noData;
						}else{
							contribs[iRow][iCol]=(short)0;
						}
					}
				}

				if (iRow == _maxRow && iCol == _maxCol) {
					//cout<<"cccc:"<<contribs[iRow][iCol];
					MPI_Barrier(MPI_COMM_WORLD);
					num = 5;
					Termination = 0;

				}
				return true;
			}
			if (contribs[iRow][iCol] < 0 && !(iRow == _maxRow && iCol == _maxCol)) {
				return true;
			}
			if(num==5){
				if(iRow == _maxRow && iCol == _maxCol)
				{
					//tiffIO iii(treefile, FLOAT_TYPE);
					long myNumLinks;
					long myNumPoints;
					long totalNumLinks;
					long totalNumPoints;
					long relativePointStart;//if 9 then line 9 is first to be filled by me

					double **PointXY;
					float **PointElevArea;

					long **LinkIdU1U2DMagShapeidCoords;
					double **LinkElevUElevDLength;
					long *buf;
					long *ptr;
					int place;
					int bsize = sizeof(long)+ MPI_BSEND_OVERHEAD;  
					buf = new long[bsize];
					getNumLinksAndPoints(myNumLinks,myNumPoints);
					//MPI_Barrier(MCW);

					int i = 0;
					if(myNumLinks > 0){
						LinkIdU1U2DMagShapeidCoords = new long*[myNumLinks];
						for(i=0;i<myNumLinks;i++)
							LinkIdU1U2DMagShapeidCoords[i] = new long[9];

						LinkElevUElevDLength = new double*[myNumLinks];
						for(i = 0;i<myNumLinks;i++)	
							LinkElevUElevDLength[i] = new double[3];

						if(myNumPoints>0)
						{
							PointXY = new double*[myNumPoints];
							for(i=0;i<myNumPoints;i++)
								PointXY[i] = new double[2];

							PointElevArea = new float*[myNumPoints];
							for(i=0;i<myNumPoints;i++)
								PointElevArea[i] = new float[3];
							
							setLinkInfo(LinkIdU1U2DMagShapeidCoords,LinkElevUElevDLength,PointXY,PointElevArea,_cellSize*_cellSize,&elevIO);
						}
					}

					int ibuf;
					int procNumPoints;
					int totalPoints = 0;
					float pbuf[3];
					double pbuf1[2];
					MPI_Status mystatus;

					int numPointsPrinted =0;
					int procNumLinks = 0;
					long treeBuf[9];

					int ilink, ipoint;
					if(rank==0){//open output point file
						FILE *fTreeOut;
						fTreeOut = fopen(treefile,"w");// process 0 writes all of its stuff
						FILE *fout;
						fout = fopen(coordfile,"w");

						//  Open shapefile 
						//need spatial refeence information which is stored in the tiffIO object
						createStreamNetShapefile(streamnetsrc,streamnetlyr,hSRSraster); // need raster spatail information for creating spatial reference in the shapefile
						long ndots=100/size+4;  // number of dots to print per process
						long nextdot=0;
						long dotinc=myNumLinks/ndots;

						for(ilink=0;ilink<myNumLinks;ilink++){//only once per link
							fprintf(fTreeOut,"\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",LinkIdU1U2DMagShapeidCoords[ilink][0],LinkIdU1U2DMagShapeidCoords[ilink][1],LinkIdU1U2DMagShapeidCoords[ilink][2],LinkIdU1U2DMagShapeidCoords[ilink][3],LinkIdU1U2DMagShapeidCoords[ilink][4],LinkIdU1U2DMagShapeidCoords[ilink][5],LinkIdU1U2DMagShapeidCoords[ilink][6],LinkIdU1U2DMagShapeidCoords[ilink][7],LinkIdU1U2DMagShapeidCoords[ilink][8]);
							//fflush(fTreeOut);
							long i1=LinkIdU1U2DMagShapeidCoords[ilink][1];
							long i2=LinkIdU1U2DMagShapeidCoords[ilink][2];
							float *lengthd, *elev, *area;
							double *pointx, *pointy;
							lengthd = new float[i2-i1+1];
							elev = new float[i2-i1+1];
							area = new float[i2-i1+1];
							pointx = new double[i2-i1+1];
							pointy = new double[i2-i1+1];
							for(ipoint=i1;ipoint<=i2;ipoint++){
								fprintf(fout,"\t%f\t%f\t%f\t%f\t%f\n",PointXY[ipoint][0],PointXY[ipoint][1],PointElevArea[ipoint][0],PointElevArea[ipoint][1],PointElevArea[ipoint][2]);
								lengthd[ipoint-i1]=PointElevArea[ipoint][0];
								elev[ipoint-i1]=PointElevArea[ipoint][1];
								area[ipoint-i1]=PointElevArea[ipoint][2];
								pointx[ipoint-i1]=PointXY[ipoint][0];
								pointy[ipoint-i1]=PointXY[ipoint][1];
							}
							//  Write shape
							long cnet[9];
							for(int iref=0; iref<9; iref++)
								cnet[iref]=LinkIdU1U2DMagShapeidCoords[ilink][iref];
							reachshape(cnet,lengthd,elev,area,pointx,pointy,i2-i1+1,srcIO);

							delete lengthd; // DGT to free memory
							delete elev;
							delete area;
							delete pointx;
							delete pointy;
							//if(ilink > nextdot)  // Indicating progress
							//{
							//	fprintf(stderr,".");
							//	fflush(stderr);
							//	nextdot=nextdot+dotinc;
							//}


						} // process 0 recvs data from other processes, writes to file
						if(myNumLinks == 0)
							numPointsPrinted = 0;
						else
							numPointsPrinted = LinkIdU1U2DMagShapeidCoords[myNumLinks-1][2]+1;
						for(i=1;i<size;i++){// send message to next process to wake it up
							MPI_Send(&ibuf,1,MPI_INT,i,0,MPI_COMM_WORLD);// get numpoints from that process
							MPI_Recv(&procNumLinks,1,MPI_INT,i,0,MPI_COMM_WORLD,&mystatus);//get points one at a time and print them to file
							if(procNumLinks > 0)
							{
								dotinc=procNumLinks/ndots;  // For tracking progress
								nextdot=0;
								for(ilink=0;ilink<procNumLinks;++ilink){
									MPI_Recv(&treeBuf,9,MPI_LONG,i,1,MPI_COMM_WORLD,&mystatus);
									fprintf(fTreeOut,"\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",treeBuf[0],treeBuf[1]+numPointsPrinted,treeBuf[2]+numPointsPrinted,treeBuf[3],treeBuf[4],treeBuf[5],treeBuf[6],treeBuf[7],treeBuf[8]);

									MPI_Recv(&procNumPoints,1,MPI_INT,i,0,MPI_COMM_WORLD,&mystatus);//get points one at a time and print them to file
									// Variables for shape
									float *lengthd, *elev, *area;
									double *pointx, *pointy;
									lengthd = new float[procNumPoints];
									elev = new float[procNumPoints];
									area = new float[procNumPoints];
									pointx = new double[procNumPoints];
									pointy = new double[procNumPoints];

									for(ipoint=0;ipoint<procNumPoints;++ipoint){
										MPI_Recv(&pbuf1,2,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&mystatus);
										MPI_Recv(&pbuf,3,MPI_FLOAT,i,1,MPI_COMM_WORLD,&mystatus);
										fprintf(fout,"\t%f\t%f\t%f\t%f\t%f\n",pbuf1[0],pbuf1[1],pbuf[0],pbuf[1],pbuf[2]); 

										lengthd[ipoint]=pbuf[0];
										elev[ipoint]=pbuf[1];
										area[ipoint]=pbuf[2];
										pointx[ipoint]=pbuf1[0];
										pointy[ipoint]=pbuf1[1];
									}
									//  Write shape
									reachshape(treeBuf,lengthd,elev,area,pointx,pointy,procNumPoints,srcIO);
									delete lengthd; // DGT to free memory
									delete elev;
									delete area;
									delete pointx;
									delete pointy;	
									//if(ilink > nextdot)  // Indicating progress
									//{
									//	fprintf(stderr,".");
									//	fflush(stderr);
									//	nextdot=nextdot+dotinc;
									//}

								}
								//  DGT moved line below to out of the loop so this only increments on the last link - and only if there is one
								numPointsPrinted += treeBuf[2]+1;//might need adjustmetn JJN
							}
						} 
						fclose(fTreeOut);
						fclose(fout);
						//SHPClose(shp1);
						//DBFClose(dbf1);
						OGR_DS_Destroy( hDS1 ); // destroy data source
						/*
						shp1.close(streamnetshp);
						*/
					}else{//other processes send their stuff to process 0
						//first, wait to recieve word from process 0
						MPI_Recv(&ibuf,1,MPI_INT,0,0,MPI_COMM_WORLD,&mystatus);//then, send myNumPoints
						MPI_Send(&myNumLinks,1,MPI_INT,0,0,MPI_COMM_WORLD);//pack up each point into a buffer, and send it
						for(ilink=0;ilink<myNumLinks;++ilink){
							treeBuf[0] = LinkIdU1U2DMagShapeidCoords[ilink][0];
							treeBuf[1] = LinkIdU1U2DMagShapeidCoords[ilink][1];
							treeBuf[2] = LinkIdU1U2DMagShapeidCoords[ilink][2];
							treeBuf[3] = LinkIdU1U2DMagShapeidCoords[ilink][3];
							treeBuf[4] = LinkIdU1U2DMagShapeidCoords[ilink][4];
							treeBuf[5] = LinkIdU1U2DMagShapeidCoords[ilink][5];
							treeBuf[6] = LinkIdU1U2DMagShapeidCoords[ilink][6];
							treeBuf[7] = LinkIdU1U2DMagShapeidCoords[ilink][7];
							treeBuf[8] = LinkIdU1U2DMagShapeidCoords[ilink][8];
							MPI_Send(&treeBuf,9,MPI_LONG,0,1,MPI_COMM_WORLD);
							int ptsInLink=treeBuf[2]-treeBuf[1]+1;
							MPI_Send(&ptsInLink,1,MPI_INT,0,0,MPI_COMM_WORLD);//pack up each point into a buffer, and send it
							for(ipoint=treeBuf[1];ipoint<=treeBuf[2];++ipoint){
								pbuf1[0]=PointXY[ipoint][0];
								pbuf1[1]=PointXY[ipoint][1];
								pbuf[0]=PointElevArea[ipoint][0];
								pbuf[1]=PointElevArea[ipoint][1];
								pbuf[2]=PointElevArea[ipoint][2];
								MPI_Send(&pbuf1,2,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
								MPI_Send(&pbuf,3,MPI_FLOAT,0,1,MPI_COMM_WORLD);
								//		MPI_Recv(&ibuf,1,MPI_INT,0,0,MCW,&mystatus);
							}
						}
					}  // just be sure we're all here
					MPI_Barrier(MPI_COMM_WORLD); 
				}


			}

			if(num==5){
				if(!contribs[iRow][iCol])
				{
					d=dir[iRow][iCol];
					tRow=iRow+d1[d];
					tCol=iCol+d2[d];
					wsGrid[iRow][iCol]=wsGrid[tRow][tCol];
					contribs[iRow][iCol]=-1;
				}


				if (iRow == _maxRow && iCol == _maxCol) {

					MPI_Barrier(MPI_COMM_WORLD);
					bool flag=true;
					for (int i = minRow; i <= _maxRow; ++i) {
						for (int j = minCol; j <= _maxCol; ++j) {
							if (contribs[i][j] ==1) {
								d=dir[i][j];
								tRow=d1[d]+i;
								tCol=d2[d]+j;
								if(contribs[tRow][tCol]==-1)
									contribs[i][j] =0;
								//cout<<"row:"<<i<<" col:"<<j<<" con:"<<contribs[i][j]<<endl;
								Termination = 0;
								//num=1;
								flag=false;
								//num++;
								//cout<<num;
								//return true;
							}
						}
					}
					//if(flag)
					//	num=2;
				}
				return true;
			}
		}break;
	default:
		cout<<"overnum:"<<num<<endl;
		break;
	}





	return true;
}



//another method, only serial
//bool SCAOperator::Operator(const CellCoord &coord,bool operFlag)
//{
//	CellSpace<double> &d8L = *(_pD8Layer->cellSpace());
//	CellSpace<double> &scaL = *(_pSCALayer->cellSpace());
//	CellSpace<double> &degreeL = *(_degreeLayer.cellSpace());
//	Neighborhood<double>& nbrhoodD = *(_pD8Nbrhood);
//	int r = _pD8Layer->_pMetaData->myrank;
//	int iRow = coord.iRow();
//	int iCol = coord.iCol();
//	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
//
//	if( num ==0 ){
//		degreeL[iRow][iCol] = 0;	//init
//		scaL[iRow][iCol] = 1;
//
//		if( iRow == _maxRow && iCol == _maxCol ){
//			MPI_Barrier(MPI_COMM_WORLD);
//			num = 1;
//			Termination = 0;
//		}
//		return true;
//	}
//
//	if( num==1 ){
//		if( d8L[iRow][iCol]>=0 && d8L[iRow][iCol]<9 && d8L[iRow][iCol]!=4 ){
//			int objRow = iRow + ((int)d8L[iRow][iCol])/3 -1;
//			int objCol = iCol + ((int)d8L[iRow][iCol])%3 -1;
//			if( degreeL[objRow][objCol]!=-1 )
//				degreeL[objRow][objCol]++;
//		}else{
//			degreeL[iRow][iCol] = -1;	//done
//		}
//
//		if( iRow == _maxRow && iCol == _maxCol ){
//			MPI_Barrier(MPI_COMM_WORLD);
//			num = 2;
//			Termination = 0;
//		}
//		return true;
//	}
//
//	if( degreeL[iRow][iCol] == 0 ){
//		int objRow = iRow + ((int)d8L[iRow][iCol])/3 -1;
//		int objCol = iCol + ((int)d8L[iRow][iCol])%3 -1;
//		scaL[objRow][objCol] += scaL[iRow][iCol];
//		degreeL[objRow][objCol]--;
//		degreeL[iRow][iCol] = -1;
//
//		//if( fabs(scaL[iRow][iCol] - _noData) >Eps ){
//		//	scaL[iRow][iCol] = scaL[iRow][iCol] / _cellSize;
//		//}
//		//scaL[iRow][iCol] = degreeL[iRow][iCol];
//		Termination = 0;
//	}
//
//	return true;
//}
