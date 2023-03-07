#include"wsOperator.h"

void wsOperator::
netLayer(RasterLayer<double>& layerD) {
    _pnetLayer = &layerD;
	_pnetNbrhood = layerD.nbrhood();
	_noData = _pnetLayer->_pMetaData->noData;
	_maxrow = _pnetLayer->_pMetaData->_localworkBR.maxIRow();
    _maxcol = _pnetLayer->_pMetaData->_localworkBR.maxICol();
	_rank=_pnetLayer->_pMetaData->myrank;

	_contribs.copyLayerInfo(*_pnetLayer);

	_cellSize = _pnetLayer->_pMetaData->cellSize;
    Configure(_pnetLayer, false);
	Configure(&_contribs, true);

}

void wsOperator::dirLayer(RasterLayer<double>& layerD) {
    _pDirLayer = &layerD;
    Configure(_pDirLayer, false);
}

void wsOperator::wsLayer(RasterLayer<double>& layerD) {
    _pWatershed = &layerD;
    Configure(_pWatershed, true);
}
bool wsOperator::isTermination() {
    //num--;
    //return num > 0;
	return true;
}

bool wsOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& net = *(_pnetLayer->cellSpace());
	CellSpace<double>& contribs = *(_contribs.cellSpace());
	CellSpace<double>& ws = *(_pWatershed->cellSpace());
	CellSpace<double>& dir = *(_pDirLayer->cellSpace());
	Neighborhood<double>& nbrhoodD = *(_pnetNbrhood);
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	int iRow = coord.iRow();
    int iCol = coord.iCol();
	int d1[9]={-1,-1,-1,0,0,0,1,1,1};
	int d2[9]={-1,0,1,-1,0,1,-1,0,1};
	
	short d;
	int nexti,nextj;
	int minRow = _pDirLayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pDirLayer->_pMetaData->_localworkBR.minICol();
	double starttime;
	double endtime;

	if(num==0)
	{
		
		if(!fabs(net[iRow][iCol] - _noData) < Eps&&!fabs(dir[iRow][iCol]-_noData)<Eps&&dir[iRow][iCol]>=0){
			contribs[iRow][iCol]=0;
			ws[iRow][iCol]=0;
		}
		else{
			contribs[iRow][iCol]=-2;//-2 means ending cal.,
			ws[iRow][iCol]=_noData;
		}
		
		if (iRow == _maxrow && iCol == _maxcol) {
			MPI_Barrier(MPI_COMM_WORLD);
			num = 1;
			Termination = 0;
        }
		return true;
	}
	
	int tRow,tCol;
	if(num==1){
		if(net[iRow][iCol]>=1){	
			for (tRow = iRow - 1; tRow <= iRow + 1; tRow++) {
				for (tCol = iCol - 1; tCol <= iCol + 1; tCol++) {
					if(!fabs(dir[tRow][tCol]-_noData)<Eps&&dir[tRow][tCol]>=0){
						d=dir[tRow][tCol];
						if (d1[d]+tRow==iRow && d2[d]+tCol==iCol&&d + dir[iRow][iCol] != 8&& d != 4 ) {
							if(net[tRow][tCol]>=1)contribs[iRow][iCol]++;
						}
					}
				}
			}
		}
		d=dir[iRow][iCol];
		tRow=d1[d]+iRow;
		tCol=d2[d]+iCol;
		if(contribs[iRow][iCol]>0&&contribs[tRow][tCol]==-2)
			n2=n2+1;
		if(contribs[iRow][iCol]>1){
			for (tRow = iRow - 1; tRow <= iRow + 1; tRow++) {
				for (tCol = iCol - 1; tCol <= iCol + 1; tCol++) {
					if(!fabs(dir[tRow][tCol]-_noData)<Eps&&dir[tRow][tCol]>=0&&contribs[tRow][tCol]>0){
						d=dir[tRow][tCol];
						if (d1[d]+tRow==iRow && d2[d]+tCol==iCol&&d + dir[iRow][iCol] != 8&& d != 4 ) {
							n1=n1+1;
						}
					}
				}
			}
		}
		if (iRow == _maxrow && iCol == _maxcol) {
			n1=n1+n2;
            MPI_Barrier(MPI_COMM_WORLD);
			MPI_Request send_request, recv_request; 
			if(rank!=0){
				MPI_Recv(&n2,1,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				cout<<"- Process "<< rank <<" received " << n2 << " from Process " << rank-1 << endl; 
				n1=n2+n1;
			}
			MPI_Send(&n1,1,MPI_INT,(rank+1)%size,0,MPI_COMM_WORLD);
			cout << "-> Process " << rank << " sent " << n1 << " to Process"<< (rank+1)%size << endl; 	
			if(rank==0)
				n2=0;
			//cout<<"f_n1:"<<n1<<endl;
            num = 2;
            Termination = 0;	
        }
		return true;
	}
	
	if (contribs[iRow][iCol] < 0 && !(iRow == _maxrow && iCol == _maxcol)) {
        return true;
    }
	

	if(num>=2){
		d=dir[iRow][iCol];
		tRow=d1[d]+iRow;
		tCol=d2[d]+iCol;
		if(contribs[iRow][iCol]>0&&(contribs[tRow][tCol]==-2||contribs[tRow][tCol]==-3||contribs[tRow][tCol]>1)){
			ws[iRow][iCol]=n2+1;
			n2=n2+1;
			contribs[iRow][iCol]=-1;
		}else if(contribs[tRow][tCol]==-1||contribs[tRow][tCol]==-3){
			ws[iRow][iCol]=ws[tRow][tCol];
			if(contribs[iRow][iCol]<2)
				contribs[iRow][iCol]=-1;
			else
				contribs[iRow][iCol]=-3;
		}else if(contribs[iRow][iCol]==0&&contribs[tRow][tCol]==-2){
			ws[iRow][iCol]=ws[tRow][tCol];
			contribs[iRow][iCol]=-2;
		}
		
		if (iRow == _maxrow && iCol == _maxcol) {
			for (int i = minRow; i <= _maxrow; ++i) {
				for (int j = minCol; j <= _maxcol; ++j) {
					if (contribs[i][j] >=0) {
						//if(contribs[i][j] ==2)
							//contribs[i][j]=-1;
						Termination = 0;
						//num++;
					}		
				} 
            }		
		}       
	}

	return true;
}




			
		
		

