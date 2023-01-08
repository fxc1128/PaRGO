#include"operator.h"

void operator::
demLayer(RasterLayer<double>& layerD) {
    _pDEMLayer = &layerD;
    cellSize = _pDEMLayer->_pMetaData->cellSize;
    noData = _pDEMLayer->_pMetaData->noData;
    Configure(_pDEMLayer, false);
}

void Operator::ssLayer(RasterLayer<double>& layerD) {
    _pssLayer = &layerD;
    Configure(_pssLayer, false);
}

void Operator::smoLayer(RasterLayer<double>& layerD) {
    _psmoLayer = &layerD;
    Configure(_psmoLayer, false);
}

bool Operator::isTermination() {
    num--;
    return num > 0;
}

bool Operator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& ss = *(_pssLayer->cellSpace());
	CellSpace<double>& smo_dem = *(_psmoLayer->cellSpace());
    //Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
    //int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	//coord=处理单元坐标
    int iRow = coord.iRow();
    int iCol = coord.iCol();

    double d=0.;

    int k = 0;
    int tag = 0;
	int dis = 0;
	float p[3]={0.4,0.1,0.05};
    for(int iRow1 = iRow - iNeighborCells; iRow1 <= iRow + iNeighborCells; iRow1++) {
        for (int iCol1 = iCol - iNeighborCells; iCol1 <= iCol + iNeighborCells; iCol1++) {
            d = dem[iRow1][iCol1];
			dis =fabs(iRow1-iRow)+fabs(iCol1-iCol);
            if ( dis> 1) {
                elevsum+=p[2]*d[k];
				wsum+=p[2];
            }
			else if(dis==1){
				elevsum+=p[1]*d[k];
				wsum+=p[1];
			}
			else{
				elevsum+=p[0]*d[k];
				wsum+=p[0];
			}
			if (fabs(d - noData) < Eps) {
                tag = 1;
            }
		}
	}
	if(!tag){
		elevsum=elevsum/wsum;
		smo_elev[iRow][iCol]=elevsum;
	}
	else{
		smo_elev[iRow][iCol]=dem[iRow][iCol];
	}

//////

	int k,ik,jk,jomax,iomax,bound;

	long x=iCol;
	long y=iRow-1;
	//-- Put smoothed elevations back in elevation grid--						
	float emax=smo_elev[x][y];
	int iomax=0;
	int jomax=0;   
	int bound=0; 
	//--Calculate Streams--
  	  //for(y=-1; y < elevny; y++){
		 // for(x=0; x < elevnx-1; x++)
		 // {
			   					/*  .false.  *		
									/*  --FIRST PASS FLAG MAX ELEVATION IN GROUP OF FOUR  */
	for(int ik=0; ik<2; ik++)
	{
		for(int jk=1-ik; jk < 2; jk++)
		{
			d=dem[x+jk][y+ik];
			if(fabs(d - noData) < Eps)
				bound= 1; 
			else if(d > emax)
			{
				emax=d;
				iomax=ik;
				jomax=jk;
			}
		}
	}
	ss[x+jomax][y+iomax]=(short)0;
	if(bound == 1)
	{
		for(ik=0; ik < 2; ik++)
			for(jk=0; jk< 2; jk++)
			{
				ss[x+jk][y+ik]=(short)0;
			}
		 			/*  c---Unflag max pixel */
			
	}else{/*  i.e. unflag flats.  */
		for(ik=0; ik < 2; ik++)
			for(jk=0; jk< 2; jk++)
			{
				if(dem[x+jk][y+ik] == emax)
				{
					ss[x+jk][y+ik]=(short)0;
				}
			}
	}
			
		
		

