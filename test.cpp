#include <iostream>
#include <opencv2/opencv.hpp>

#include "image_process.hpp"
#include "random_star_catalog.hpp"
#include "coordinate_transfer.hpp"
#include "visual_position_calculate.hpp"

using namespace std;
using namespace cv;


int main()
{
    visual_position_calculate *vpc = new visual_position_calculate();
    char path[256];
    vpc->GetAppPath(path);
    printf("%s\n",path);

    //记录当前时间，用于程序耗时记录
    win_time_val_t wintv;  
    win_time_t wintime;
    win_gettimeofday(&wintv);  
	win_time(&wintv,&wintime);
    printf("%d\n",wintime.year);
    printf("%d\n",wintime.mon);
    printf("%d\n",wintime.day);
    printf("%d\n",wintime.wday);
    printf("%d\n",wintime.hour);
    printf("%d\n",wintime.min);
    printf("%d\n",wintime.sec);
    printf("%d\n",wintime.msec);


    RANDOMSTARCATALOG* testStarCatalog = new RANDOMSTARCATALOG();

    unsigned int catalogCapacity = 20;
    float magnitudeMAX_user = 6.0f;
    double a = 0.0; 
    double b = 10.0;
    double c = 70.0; 
    double d = 80.0;

    testStarCatalog->CapacitySet(catalogCapacity);
    testStarCatalog->MagnitudeMAXSet(magnitudeMAX_user);
    testStarCatalog->RangeSet(a, b, c, d);
    testStarCatalog->CatalogGenerate();

    STARINFO displayStarInfo;

    testStarCatalog->GetStarInfo(11, displayStarInfo);

    COORDINATETRANSFER* testCoordinateTrans = new COORDINATETRANSFER();
    testCoordinateTrans->SetPosTime(34.0, 113.0, wintime.year, wintime.mon, wintime.day, wintime.hour, wintime.min, wintime.sec, wintime.msec);

    


	
    
    return 0;
}