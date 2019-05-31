//random_star_catalog.hpp
//生成随机星表的类

#include <stdlib.h>
#include <ctime>

#ifndef MAX_CATALOG_CAPACITY
#define MAX_CATALOG_CAPACITY 1024 //最大星数量
#endif

struct STARINFO //星点信息
{
    double asc; //赤经
    double dec;   //赤纬
    float magnitude;   //星等
};

class RANDOMSTARCATALOG //随机星表类
{
  public:
    
    int CapacitySet(unsigned int userNum);                                                                                //设定需要生成的星数量
    inline int MagnitudeMAXSet(float magnitudeMAX_user);                                                                         //设定星等上限
    int RangeSet(double ascMIN_user, double ascMAX_user, double decMIN_user, double decMAX_user); //设定天区范围
    inline int DataReset();                                                                                                      //数据重置
    inline unsigned int GetCapacity() { return catalogCapacity; };                                                        //获得星表容量（星数量）

    int CatalogGenerate();                                 //生成随机星表
    int CatalogClear();                                    //清除星表数据
    inline bool CatalogIsReady() { return catalogReady; }; //检查星表是否就绪（存在星表数据）

    int GetStarInfo(unsigned int starNum, STARINFO &starInfo); //获得某一颗星信息

  private:
    bool catalogReady = false;                  //星表是否就绪
    unsigned int catalogCapacity;               //星表容量（星数量）
    double ascMIN;                       //赤经下限
    double ascMAX;                       //赤经上限
    double decMIN;                         //赤纬下限
    double decMAX;                         //赤纬上限
    float magnitudeMAX;                         //星等上限
    STARINFO starCatalog[MAX_CATALOG_CAPACITY]; //生成的随机星表
};


inline int RANDOMSTARCATALOG::DataReset()
{
    //所有数据置零
    catalogCapacity = 0;
    ascMIN = 0.0;
    ascMAX = 0.0;
    decMIN = 0.0;
    decMAX = 0.0;
    magnitudeMAX = 0.0f;
    catalogReady = false;
    return 0;
}

int RANDOMSTARCATALOG::CapacitySet(unsigned int userNum)
{
    //设定容量超出最大限制，返回-1
    if (userNum > MAX_CATALOG_CAPACITY)
    {
        return -1;
    }

    //新设定容量小于原容量，将相差部分星表数据置零
    if (userNum < catalogCapacity)
    {
        unsigned int ii;

        for (ii = userNum; ii < catalogCapacity; ii++)
        {
            starCatalog[ii].dec = 0.0;
            starCatalog[ii].asc = 0.0;
            starCatalog[ii].magnitude = 0.0f;
        }
    }

    //设定新容量
    catalogCapacity = userNum;
    return 0;
}

inline int RANDOMSTARCATALOG::MagnitudeMAXSet(float magnitudeMAX_user)
{
    //设定星等上限
    magnitudeMAX = magnitudeMAX_user;
    return 0;
}

int RANDOMSTARCATALOG::RangeSet(double ascMIN_user, double ascMAX_user, double decMIN_user, double decMAX_user)
{
    // 赤经 最小值大于最大值，设定非法，返回-1
    if (ascMAX_user < ascMIN_user)
    {
        return -1;
    }

    // 赤纬 最小值大于最大值，设定非法，返回-2
    if (decMAX_user < decMIN_user)
    {
        return -2;
    }

    //设定赤经赤纬范围（上下限）
    ascMIN = ascMIN_user;
    ascMAX = ascMAX_user;
    decMIN = decMIN_user;
    decMAX = decMAX_user;
    return 0;
}

int RANDOMSTARCATALOG::CatalogClear()
{
    //清楚星表数据
    unsigned int ii;
    for (ii = 0; ii < catalogCapacity; ii++)
        ;
    {
        starCatalog[ii].dec = 0.0;
        starCatalog[ii].asc = 0.0;
        starCatalog[ii].magnitude = 0.0f;
    }

    //其他数据置零
    DataReset();
    return 0;
}

int RANDOMSTARCATALOG::CatalogGenerate()
{
    //数据非法，可能是相关数据未设定，返回-1
    if ((catalogCapacity == 0) || (ascMAX == 0.0) || (magnitudeMAX == 0.0f))
    {
        return -1;
    }

    unsigned int ii;

    //随机生成每颗星的 赤经 数据
    srand((unsigned)time(NULL));
    for (ii = 0; ii < catalogCapacity; ii++)
    {
        starCatalog[ii].asc = (ascMAX - ascMIN) * rand() / double(RAND_MAX) + ascMIN;
    }

    //随机生成每颗星的 赤纬 数据
    srand((unsigned)time(NULL));
    for (ii = 0; ii < catalogCapacity; ii++)
    {
        starCatalog[ii].dec = (decMAX - decMIN) * rand() / double(RAND_MAX) + decMIN;
    }

    //随机生成每颗星的 星等 数据，取小数点后两位
    srand((unsigned)time(NULL));
    for (ii = 0; ii < catalogCapacity; ii++)
    {
        starCatalog[ii].magnitude = magnitudeMAX * (rand() % 100) / 100.0f;
    }

    //星表就绪标识设为true
    catalogReady = true;
    return 0;
}

int RANDOMSTARCATALOG::GetStarInfo(unsigned int starNum, STARINFO &starInfo)
{
    //访问位置超限，返回-1
    if (starNum > catalogCapacity)
    {
        return -1;
    }

    //将星点信息传引用传递出来
    starInfo = starCatalog[starNum];
    return 0;
}