#include <iostream>
#include <time.h>
using namespace std;

int main()
{

	clock_t start,end;
	start=clock();
	system("G: & cd G:\\egc\\logandemo\\logant\\test\\0324 & mpiexec -n 4 G:\\EGC\\PaRGO-dev\\build\\apps\\hydrology\\Debug\\scad8.exe ttd81.tif G:\\EGC\\PaRGO-dev\\neighbor\\moore.nbr prglscad82.tif ");
	end=clock();
	cout<<"\ntime: = "<<((double)end-start)/CLK_TCK;
	system("pause");
	return 0;
}