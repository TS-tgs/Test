//coordinate_transfer.hpp
//地平、赤道、图像空间坐标转换

#include <cmath>

#ifndef PI
#define PI (139755218526789.0 / 44485467702853.0)
#endif

struct SPHERICALCOORDINATE
{
    double longtitude;
    double latitude;
};  //球面坐标结构体，单位为弧度

struct IMAGECOORDINATE
{
    double X;
    double Y;
};  //图像坐标结构体，单位为像素

struct CAMERACOORDINATE
{
    double horizon;
    double vertical;
};  //相机（传感器）坐标结构图，单位为μm

class COORDINATETRANSFER
{
    public:
    
    //获得某一星点各项坐标的函数
    SPHERICALCOORDINATE GetEquatorialCoordinate(int num); 
    SPHERICALCOORDINATE GetHorizonalCoordinate(int num);
    IMAGECOORDINATE GetImageCoordinate(int num); 
    CAMERACOORDINATE GetCameraCoordinate(int num); 

    //获得各个成员变量的函数
    inline unsigned int GetResolutionX(){ return resolutionX; };
    inline unsigned int GetResolutionY(){ return resolutionY; };
    inline float GetSensorWith(){ return sensorWith; };
    inline float GetSensorHeight(){ return sensorHeight; };

    int InportEquatorialCoordinate(RANDOMSTARCATALOG* randomStarCatalog);
    int SetParameters(unsigned int resX, unsigned int resY, float sensorW, float sensorH, float flength);
    void SetPosTime(double la, double lo, int y, int m, int d, int h, int mi, int s, int ms);

    private:

    SPHERICALCOORDINATE equatorialCoordinate[MAX_CATALOG_CAPACITY];   //赤道坐标
    SPHERICALCOORDINATE horizonalCoordinate[MAX_CATALOG_CAPACITY];    //地平坐标
    IMAGECOORDINATE imageCoordinate[MAX_CATALOG_CAPACITY];            //图像坐标
    CAMERACOORDINATE cameraCoordinate[MAX_CATALOG_CAPACITY];          //相机坐标
    unsigned int resolutionX;                   //星图水平分辨率
    unsigned int resolutionY;                   //星图垂直分辨率
    float sensorWith;                           //传感器宽度，单位μm
    float sensorHeight;                         //传感器高度，单位μm
    float focalLength;                          //焦距，单位mm
    unsigned int capacity;                      //容量
    bool isReady = false;                       //数据就绪标识
    double lat; //纬度
    double lon; //经度
    int year;   //年
    int mon;    //月
    int day;    //日
    int hour;   //时
    int min;    //分
    int sec;    //秒
    int msec;   //毫秒

};

void COORDINATETRANSFER::SetPosTime(double la, double lo, int y, int m, int d, int h, int mi, int s, int ms)
{
    lat = la; lon = lo; year = y; mon = m; day = d; hour = h; min = mi; sec = s; msec = ms;
};

SPHERICALCOORDINATE COORDINATETRANSFER::GetEquatorialCoordinate(int num)
{
    if (num > MAX_CATALOG_CAPACITY)
    {
        num = MAX_CATALOG_CAPACITY;
    }

    return equatorialCoordinate[num]; 
}; 

SPHERICALCOORDINATE COORDINATETRANSFER::GetHorizonalCoordinate(int num)
{
    if (num > MAX_CATALOG_CAPACITY)
    {
        num = MAX_CATALOG_CAPACITY;
    }

    return horizonalCoordinate[num]; 
};

IMAGECOORDINATE COORDINATETRANSFER::GetImageCoordinate(int num)
{
    if (num > MAX_CATALOG_CAPACITY)
    {
        num = MAX_CATALOG_CAPACITY;
    }

    return imageCoordinate[num]; 
}; 

CAMERACOORDINATE COORDINATETRANSFER::GetCameraCoordinate(int num)
{
    if (num > MAX_CATALOG_CAPACITY)
    {
        num = MAX_CATALOG_CAPACITY;
    }

    return cameraCoordinate[num]; 
}; 

int COORDINATETRANSFER::InportEquatorialCoordinate(RANDOMSTARCATALOG* randomStarCatalog)
{
    capacity = randomStarCatalog->GetCapacity();

    int ii;
    STARINFO tmpStarInfo;
    for (ii = 0; ii < capacity; ii++)
    {
        randomStarCatalog->GetStarInfo(ii, tmpStarInfo);
        equatorialCoordinate[ii].latitude = tmpStarInfo.dec;
        equatorialCoordinate[ii].longtitude = tmpStarInfo.asc;
    }

    isReady = true;
    return 0;
};