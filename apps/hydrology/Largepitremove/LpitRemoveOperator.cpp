#include"LpitRemoveOperator.h"

void PitRemoveOperator::
demLayer(RasterLayer<double>& layerD) {
    _pDEMLayer = &layerD;
    _pDEMNbrhood = layerD.nbrhood();
    cellSize = _pDEMLayer->_pMetaData->cellSize;
    noData = _pDEMLayer->_pMetaData->noData;
    _xSize = _pDEMLayer->_pMetaData->_localdims.nRows();
    _ySize = _pDEMLayer->_pMetaData->_localdims.nCols();
    //cout<<noData<<endl;

    Configure(_pDEMLayer, false);
}
void PitRemoveOperator::wslayer(RasterLayer<double>& layerD) {
    _pwsLayer = &layerD;
    Configure(_pwsLayer, false);
}

void PitRemoveOperator::outLayer(RasterLayer<double>& layerD) {
    _pOutLayer = &layerD;
    Configure(_pOutLayer, false);
}
void PitRemoveOperator::wdemLayer(RasterLayer<double>& layerD) {
    _pwDEMLayer = &layerD;

    Configure(_pwDEMLayer, true);
}

bool PitRemoveOperator::isTermination() {
    return flag;
}
void PitRemoveOperator::init(RasterLayer<double>& layerD){
	_pwDEMLayer=&layerD;
	CellSpace<double>& wdem = *(_pwDEMLayer->cellSpace());
	
	int maxRow = _pwDEMLayer->_pMetaData->_localworkBR.maxIRow();
    int maxCol = _pwDEMLayer->_pMetaData->_localworkBR.maxICol();
	int minRow = _pwDEMLayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pwDEMLayer->_pMetaData->_localworkBR.minICol();
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			wdem[i][j]=0;
		}	
	}
}
void PitRemoveOperator::writesca(RasterLayer<double>& layerD,CellCoord nw){
	_pwkLayer=&layerD;
	CellSpace<double>& wk = *(_pwkLayer->cellSpace());
	CellSpace<double>& wdem = *(_pwDEMLayer->cellSpace());
	//CellSpace<double>& dir = *(_pD8Layer->cellSpace());
	int maxRow = _pwkLayer->_pMetaData->_localworkBR.maxIRow();
    int maxCol = _pwkLayer->_pMetaData->_localworkBR.maxICol();
	int d1[9]={-1,-1,-1,0,0,0,1,1,1};
	int d2[9]={-1,0,1,-1,0,1,-1,0,1};
	int minRow = _pwkLayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pwkLayer->_pMetaData->_localworkBR.minICol();
	int gi = nw.iRow()-1;//_pwkLayer->_pMetaData->_MBR.minIRow();
	int gj = nw.iCol();//_pwkLayer->_pMetaData->_MBR.minICol();
	//cout<<"globali:"<<gi<<" globalj:"<<gj<<endl;
	//cout<<"scalayer"<<_pSCALayer->_pMetaData->_localworkBR.maxIRow()<<endl;
	int ti,tj;
	short d;
	double s=0.;
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			//cout<<"sca:"<<sca[i+gi][j+gj]<<endl;
			//cout<<"wk:"<<wk[i][j]<<endl;
			if(fabs(wdem[i+gi][j+gj]-noData)<=Eps||wdem[i+gi][j+gj]<wk[i][j]){
				wdem[i+gi][j+gj]=wk[i][j];
				//cout<<"wk:"<<wk[i][j]<<" sca:"<<sca[i+gi][j+gj]<<endl;
			}
		}
	}
}
void PitRemoveOperator::getarea(int &minrow,int &mincol,int &nrow,int &ncol, int _g,int buf,int id) {
	CellSpace<double>& ws = *(_pwsLayer->cellSpace());

	int pi,pj,qi,qj;
	int minRow = _pwsLayer->_pMetaData->_MBR.minIRow();
	int minCol = _pwsLayer->_pMetaData->_MBR.minICol();
	int maxRow = _pwsLayer->_pMetaData->_MBR.maxIRow();
	int maxCol = _pwsLayer->_pMetaData->_MBR.maxICol();

	pi=-2;
	pj=-2;
	qi=-2;
	qj=-2;
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			if(ws[i][j]==id){
				if(i<pi||pi==-2)
					pi=i;
				if(j<pj||pj==-2)
					pj=j;
				if(i>qi||qi==-2)
					qi=i;
				if(j>qj||qj==-2)
					qj=j;
			}
		}
	}
	
	if((pi-buf)>minRow)
		pi=pi-buf;
	else
		pi=minRow+1;
	if((pj-buf)>minCol)
		pj=pj-buf;
	else
		pj=minCol+1;
	if((qi+buf)<maxRow)
		qi=qi+buf;
	else
		qi=maxRow-1;
	if((qj+buf)<maxCol)
		qj=qj+buf;
	else
		qj=maxCol-1;

	minrow=pi*_g;
	mincol=pj*_g;
	nrow=(qi-pi+1)*_g;
	ncol=(qj-pj+1)*_g;
}
bool PitRemoveOperator::Operator(const CellCoord& coord, bool operFlag) {

    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& wdem = *(_pwDEMLayer->cellSpace());
	CellSpace<double>& out = *(_pOutLayer->cellSpace());
    int iRow = coord.iRow();
    int iCol = coord.iCol();

    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
    int iNeighborCells = (int)sqrt((double)nbrhoodD.size()) / 2;

    double gap = 0.0005;

    int i, j;
    if (num == 0) {
        if (fabs(dem[iRow][iCol] - noData) > Eps) {
            wdem[iRow][iCol] = 10000.0;
			if(out[iRow][iCol]>dem[iRow][iCol])
				dem[iRow][iCol]=out[iRow][iCol];
        }
        else {
            wdem[iRow][iCol] = noData;
        }
		
        if ((iRow == _xSize - 2) && (iCol == _ySize - 2)) {
            MPI_Barrier(MPI_COMM_WORLD);
            num = 1;
            Termination = 0;
        }
    }
    else {
        if (fabs(dem[iRow][iCol] - noData) <= Eps || wdem[iRow][iCol] <= dem[iRow][iCol]) {
            return true;
        }

        for (i = iRow - iNeighborCells; i <= iRow + iNeighborCells; i++) {
            for (j = iCol - iNeighborCells; j <= iCol + iNeighborCells; j++) {
                if ((dem[iRow][iCol] >= (wdem[i][j] + gap)) || fabs(wdem[i][j] - noData) < Eps) {
                    wdem[iRow][iCol] = dem[iRow][iCol];
                    Termination = 0;
                }else if (wdem[iRow][iCol] > (wdem[i][j] + gap)) {
                    wdem[iRow][iCol] = wdem[i][j] + gap;
                    Termination = 0;
                }

            }
        }
    }

    return true;
}
