#include <windows.h>

//需要获得utc以及utc-tai
//恒星参数需要：赤经，赤纬，赤经自行，赤纬自行，视差，视向速度
#pragma once

  
#define SECS_TO_FT_MULT 10000000  
#ifndef PI
#define PI (139755218526789.0 / 44485467702853.0)
#endif
#define J2000 (2451545.0)
#define D2R (PI / 180.0)
#define MUC (2.0 * GAUSSK * GAUSSK / (CAUD * CAUD))
#define GAUSSK (0.01720209895)
#define CAUD (CKMS * 86400.0 * KM2AU)
#define CKMS (299792.458)
#define KM2AU (1.0 / 149597870.66)
#define R2A (648000.0 / PI)
#define JulCty (36525.0)
#define A2R (PI / 648000.0)
#define TWOPI  (2.0 * PI)
#define BIG_ENDIAN_TEST 0

using namespace std;


static LARGE_INTEGER base_time;  
  
typedef struct win_time_val  
{  
    /** The seconds part of the time. */  
    long    sec;  
  
    /** The miliseconds fraction of the time. */  
    long    msec;  
  
} win_time_val_t;  
  
typedef struct win_time  
{  
    /** This represents day of week where value zero means Sunday */  
    int wday;  
  
    /** This represents day of month: 1-31 */  
    int day;  
  
    /** This represents month, with the value is 0 - 11 (zero is January) */  
    int mon;  
  
    /** This represent the actual year (unlike in ANSI libc where 
     *  the value must be added by 1900). 
     */  
    int year;  
  
    /** This represents the second part, with the value is 0-59 */  
    int sec;  
  
    /** This represents the minute part, with the value is: 0-59 */  
    int min;  
  
    /** This represents the hour part, with the value is 0-23 */  
    int hour;  
  
    /** This represents the milisecond part, with the value is 0-999 */  
    int msec;  
  
}win_time_t; 

typedef double** DMatrix;

/* 1980 IAU nutation coefficients */
static const double IAU1980[106][9] = {
	{-171996., -174.2, 92025., 8.9,  0., 0., 0., 0., 1.},
	{   2062.,    0.2,  -895., 0.5,  0., 0., 0., 0., 2.},
	{     46.,     0.,   -24.,  0., -2., 0., 2., 0., 1.},
	{     11.,     0.,     0.,  0.,  2., 0.,-2., 0., 0.},
	{     -3.,     0.,     1.,  0., -2., 0., 2., 0., 2.},
	{     -3.,     0.,     0.,  0.,  1.,-1., 0.,-1., 0.},
	{     -2.,     0.,     1.,  0.,  0.,-2., 2.,-2., 1.},
	{      1.,     0.,     0.,  0.,  2., 0.,-2., 0., 1.},
	{ -13187.,   -1.6,  5736., -3.1, 0., 0., 2.,-2., 2.},
	{   1426.,   -3.4,    54., -0.1, 0., 1., 0., 0., 0.},
	{   -517.,    1.2,   224., -0.6, 0., 1., 2.,-2., 2.},
	{    217.,   -0.5,   -95.,  0.3, 0.,-1., 2.,-2., 2.},
	{    129.,    0.1,   -70.,  0.,  0., 0., 2.,-2., 1.},
	{     48.,     0.,     1.,  0.,  2., 0., 0.,-2., 0.},
	{    -22.,     0.,     0.,  0.,  0., 0., 2.,-2., 0.},
	{     17.,   -0.1,     0.,  0.,  0., 2., 0., 0., 0.},
	{    -15.,     0.,     9.,  0.,  0., 1., 0., 0., 1.},
	{    -16.,    0.1,     7.,  0.,  0., 2., 2.,-2., 2.},
	{    -12.,     0.,     6.,  0.,  0.,-1., 0., 0., 1.},
	{     -6.,     0.,     3.,  0., -2., 0., 0., 2., 1.},
	{     -5.,     0.,     3.,  0.,  0.,-1., 2.,-2., 1.},
	{      4.,     0.,    -2.,  0.,  2., 0., 0.,-2., 1.},
	{      4.,     0.,    -2.,  0.,  0., 1., 2.,-2., 1.},
	{     -4.,     0.,     0.,  0.,  1., 0., 0.,-1., 0.},
	{      1.,     0.,     0.,  0.,  2., 1., 0.,-2., 0.},
	{      1.,     0.,     0.,  0.,  0., 0.,-2., 2., 1.},
	{     -1.,     0.,     0.,  0.,  0., 1.,-2., 2., 0.},
	{      1.,     0.,     0.,  0.,  0., 1., 0., 0., 2.},
	{      1.,     0.,     0.,  0., -1., 0., 0., 1., 1.},
	{     -1.,     0.,     0.,  0.,  0., 1., 2.,-2., 0.},
	{  -2274.,   -0.2,   977., -0.5, 0., 0., 2., 0., 2.},
	{    712.,    0.1,    -7.,   0., 1., 0., 0., 0., 0.},
	{   -386.,   -0.4,   200.,   0., 0., 0., 2., 0., 1.},
	{   -301.,     0.,   129., -0.1, 1., 0., 2., 0., 2.},
	{   -158.,     0.,    -1.,   0., 1., 0., 0.,-2., 0.},
	{    123.,     0.,   -53.,   0.,-1., 0., 2., 0., 2.},
	{     63.,     0.,    -2.,   0., 0., 0., 0., 2., 0.},
	{     63.,    0.1,   -33.,   0., 1., 0., 0., 0., 1.},
	{    -58.,   -0.1,    32.,   0.,-1., 0., 0., 0., 1.},
	{    -59.,     0.,    26.,   0.,-1., 0., 2., 2., 2.},
	{    -51.,     0.,    27.,   0., 1., 0., 2., 0., 1.},
	{    -38.,     0.,    16.,   0., 0., 0., 2., 2., 2.},
	{     29.,     0.,    -1.,   0., 2., 0., 0., 0., 0.},
	{     29.,     0.,   -12.,   0., 1., 0., 2.,-2., 2.},
	{    -31.,     0.,    13.,   0., 2., 0., 2., 0., 2.},
	{     26.,     0.,    -1.,   0., 0., 0., 2., 0., 0.},
	{     21.,     0.,   -10.,   0.,-1., 0., 2., 0., 1.},
	{     16.,     0.,    -8.,   0.,-1., 0., 0., 2., 1.},
	{    -13.,     0.,     7.,   0., 1., 0., 0.,-2., 1.},
	{    -10.,     0.,     5.,   0.,-1., 0., 2., 2., 1.},
	{     -7.,     0.,     0.,   0., 1., 1., 0.,-2., 0.},
	{      7.,     0.,    -3.,   0., 0., 1., 2., 0., 2.},
	{     -7.,     0.,     3.,   0., 0.,-1., 2., 0., 2.},
	{     -8.,     0.,     3.,   0., 1., 0., 2., 2., 2.},
	{      6.,     0.,     0.,   0., 1., 0., 0., 2., 0.},
	{      6.,     0.,    -3.,   0., 2., 0., 2.,-2., 2.},
	{     -6.,     0.,     3.,   0., 0., 0., 0., 2., 1.},
	{     -7.,     0.,     3.,   0., 0., 0., 2., 2., 1.},
	{      6.,     0.,    -3.,   0., 1., 0., 2.,-2., 1.},
	{     -5.,     0.,     3.,   0., 0., 0., 0.,-2., 1.},
	{      5.,     0.,     0.,   0., 1.,-1., 0., 0., 0.},
	{     -5.,     0.,     3.,   0., 2., 0., 2., 0., 1.},
	{     -4.,     0.,     0.,   0., 0., 1., 0.,-2., 0.},
	{      4.,     0.,     0.,   0., 1., 0.,-2., 0., 0.},
	{     -4.,     0.,     0.,   0., 0., 0., 0., 1., 0.},
	{     -3.,     0.,     0.,   0., 1., 1., 0., 0., 0.},
	{      3.,     0.,     0.,   0., 1., 0., 2., 0., 0.},
	{     -3.,     0.,     1.,   0., 1.,-1., 2., 0., 2.},
	{     -3.,     0.,     1.,   0.,-1.,-1., 2., 2., 2.},
	{     -2.,     0.,     1.,   0.,-2., 0., 0., 0., 1.},
	{     -3.,     0.,     1.,   0., 3., 0., 2., 0., 2.},
	{     -3.,     0.,     1.,   0., 0.,-1., 2., 2., 2.},
	{      2.,     0.,    -1.,   0., 1., 1., 2., 0., 2.},
	{     -2.,     0.,     1.,   0.,-1., 0., 2.,-2., 1.},
	{      2.,     0.,    -1.,   0., 2., 0., 0., 0., 1.},
	{     -2.,     0.,     1.,   0., 1., 0., 0., 0., 2.},
	{      2.,     0.,     0.,   0., 3., 0., 0., 0., 0.},
	{      2.,     0.,    -1.,   0., 0., 0., 2., 1., 2.},
	{      1.,     0.,    -1.,   0.,-1., 0., 0., 0., 2.},
	{     -1.,     0.,     0.,   0., 1., 0., 0.,-4., 0.},
	{      1.,     0.,    -1.,   0.,-2., 0., 2., 2., 2.},
	{     -2.,     0.,     1.,   0.,-1., 0., 2., 4., 2.},
	{     -1.,     0.,     0.,   0., 2., 0., 0.,-4., 0.},
	{      1.,     0.,    -1.,   0., 1., 1., 2.,-2., 2.},
	{     -1.,     0.,     1.,   0., 1., 0., 2., 2., 1.},
	{     -1.,     0.,     1.,   0.,-2., 0., 2., 4., 2.},
	{      1.,     0.,     0.,   0.,-1., 0., 4., 0., 2.},
	{      1.,     0.,     0.,   0., 1.,-1., 0.,-2., 0.},
	{      1.,     0.,    -1.,   0., 2., 0., 2.,-2., 1.},
	{     -1.,     0.,     0.,   0., 2., 0., 2., 2., 2.},
	{     -1.,     0.,     0.,   0., 1., 0., 0., 2., 1.},
	{      1.,     0.,     0.,   0., 0., 0., 4.,-2., 2.},
	{      1.,     0.,     0.,   0., 3., 0., 2.,-2., 2.},
	{     -1.,     0.,     0.,   0., 1., 0., 2.,-2., 0.},
	{      1.,     0.,     0.,   0., 0., 1., 2., 0., 1.},
	{      1.,     0.,     0.,   0.,-1.,-1., 0., 2., 1.},
	{     -1.,     0.,     0.,   0., 0., 0.,-2., 0., 1.},
	{     -1.,     0.,     0.,   0., 0., 0., 2.,-1., 2.},
	{     -1.,     0.,     0.,   0., 0., 1., 0., 2., 0.},
	{     -1.,     0.,     0.,   0., 1., 0.,-2.,-2., 0.},
	{     -1.,     0.,     0.,   0., 0.,-1., 2., 0., 1.},
	{     -1.,     0.,     0.,   0., 1., 1., 0.,-2., 1.},
	{     -1.,     0.,     0.,   0., 1., 0.,-2., 2., 0.},
	{      1.,     0.,     0.,   0., 2., 0., 0., 2., 0.},
	{     -1.,     0.,     0.,   0., 0., 0., 2., 4., 2.},
	{      1.,     0.,     0.,   0., 0., 1., 0., 1., 0.}
};

// Find 1st Jan 1970 as a FILETIME   
static void get_base_time(LARGE_INTEGER *base_time)  
{  
    SYSTEMTIME st;  
    FILETIME ft;  
  
    memset(&st,0,sizeof(st));  
    st.wYear=1970;  
    st.wMonth=1;  
    st.wDay=1;  
    SystemTimeToFileTime(&st, &ft);  
      
    base_time->LowPart = ft.dwLowDateTime;  
    base_time->HighPart = ft.dwHighDateTime;  
    base_time->QuadPart /= SECS_TO_FT_MULT;  
}

int win_gettimeofday(win_time_val_t *tv)  
{  
    SYSTEMTIME st;  
    FILETIME ft;  
    LARGE_INTEGER li;  
    static char get_base_time_flag=0;  
  
    if (get_base_time_flag == 0)  
    {  
        get_base_time(&base_time);  
    }  
  
    /* Standard Win32 GetLocalTime */  
    GetLocalTime(&st);  
    SystemTimeToFileTime(&st, &ft);  
  
    li.LowPart = ft.dwLowDateTime;  
    li.HighPart = ft.dwHighDateTime;  
    li.QuadPart /= SECS_TO_FT_MULT;  
    li.QuadPart -= base_time.QuadPart;  
  
    tv->sec = li.LowPart;  
    tv->msec = st.wMilliseconds;  
  
    return 0;  
}
  
int win_time(const win_time_val_t *tv, win_time_t *time)  
{  
    LARGE_INTEGER li;  
    FILETIME ft;  
    SYSTEMTIME st;  
  
    li.QuadPart = tv->sec;  
    li.QuadPart += base_time.QuadPart;  
    li.QuadPart *= SECS_TO_FT_MULT;  
  
    ft.dwLowDateTime = li.LowPart;  
    ft.dwHighDateTime = li.HighPart;  
    FileTimeToSystemTime(&ft, &st);  
  
    time->year = st.wYear;  
    time->mon = st.wMonth;  
    time->day = st.wDay;  
    time->wday = st.wDayOfWeek;  
  
    time->hour = st.wHour;  
    time->min = st.wMinute;  
    time->sec = st.wSecond;  
    time->msec = tv->msec;  
  
    return 0;  
} 

class visual_position_calculate
{
	public:
	visual_position_calculate(int y, int m, int d, int hour, int minute, int second, int milisecond, double tai_utcData);


	int nYear;
	int nMonth;
	int nDay;
	int nHour;
	int nMinute;
	double nSecond;
	double m_utc;
	double m_tai_utc;
	bool timeReady;

    void matrix_mxmul(double a[], double b[], double c[],int m, int n, int k);//矩阵乘法函数
	int matrix_inv(double a[], int n);//矩阵求逆函数
	double matrix_sum(double a[], int m);//向量求和函数,a为向量，m为要求和的元素个数
	void matrix_mul(double a[], double b[], double c[], int m);//向量点乘函数
	bool GetAppPath(char *Path);//得到当前应用程序所在绝对路径


	void FinalCalculate(int starNum, double starRawData[][8], double starViualPos[][4]);


	/////////////////////用于恒星视位置计算的所有函数进行声明////////////////////////////

	private:

	double greenwichHourAngle;
	double jedz;
	int efirst;
	int bcoef, np, nv;
	double pc[18], vc[18], cbody[1200];
	double cbuf[15][3][8];
	double twot, dna, dt1, temp, ll, tc, vfac;
	double m_JD;        //用于行星和恒星
	double StateVector[15][15][2][6] ;

	bool m_bePlant;
	bool m_beLunar;
	bool m_beAlt;

	double   SS[3],cval[400],emrat,au;
	short int NUMDE,ipt[3][12],lpt[3];
	int  bary, bsav;
	double pvsun[6],db[1100];
	FILE * fpBinaryFile;
	char ttl[3][65];
	short tmpShort;
	int  ncon;
	char cnam[400][7];
	double tmpDouble;
	long LengthOfFile,ncoeff,BlockLength;
	int LengthOfHeader,NumBlocks;

	
	double Cal2JED(int m, int d, int y, double utc, double tai_utc);//将世界时转化为儒略日期

	bool flag_GetStateVector;//GetStateVector函数是否已执行标识
	bool GetStateVector(double jd, int targ, int cent, int recpol, double StateVector[15][15][2][6]);//得到星体的状态矢量 

	double hms2R(int a, int b, int c, int d );//将h，m，s，ms转化成一个实数——弧度
	void LightTime(double jed, int body, double earth_ssb[], double earth_hel[], double StarData[], double body_geo[], double body_hel[], double *Ltime);//将星体由几何太阳系质心位置转化到几何位置或日心位置
	void RayBend(double earth_hel[], double body_geo[], double body_hel[], double p1[]);//根据太阳的重力场来计算相对光偏差
	void SplitStateVector(double pv[], double p[], double v[]);
	void Aberrate(double p1[], double EBdot[], double p2[]);//计算光行差
	void GetPrecessParams(double jed1, double jed2, double *zeta, double *z, double *theta, double *zetadot, double *zdot, double *thetadot);//计算岁差参数
	
	bool flag_GetDpsiDeps;//GetDpsiDeps函数是否已执行标识
	void GetDpsiDeps(double jed, double *dpsi, double *deps, double *dpsidot, double *depsdot);//计算章动dpsi, deps, dpsidot, and depsdot
	
	void Obliquity(double jed1, double jed2, int m, double *obl, double *obldot);//计算平黄赤交角及它的导数
	void QRotate(double vin[], int axis, double phi, double phidot, int s, double vout[]);//矩阵旋转
	double deg(double x);//将 hh.mmssss 时间转化到 hh.hhhhh 时间
	double fix(double x);//得返回到实数的整数值
	
	bool flag_pleph;//pleph函数是否已执行标识
	bool pleph(double jd, int targ, int cent, double *rrd, int *inside);//读取JPL行星星表并给出星体相对于给定原点的位置和速度矢量
	
	void Rec2Pol(double a[], double b[]);//将直角坐标转变为极坐标
	double Vecmag(double a[]);//计算矢量的数值
	void Uvector(double a[], double unita[]);//计算单位矢量
	void Vdot(int n, double a[], double b[], double *adotb);//矢量乘法

	bool flag_FunArgIAU;//FunArgIAU函数是否已执行标识
	void FunArgIAU(double jed, double *funarg);//利用IAU表达式计算基本参数///计算章动的五个基本变量函数
	
	double amodulo(double a, double b);//计算一个在0 <= a < b范围内的a值
	DMatrix createDMatrix(int n);//创建一个n*n矩阵
	void GetQMatrix(double phi, double phidot, int axis, int s, DMatrix QMatrix);//对于给定的角度和时间倒数来计算矩阵的旋转
	void MatXVec(DMatrix a, double b[], double c[], int l, int m);//矩阵相乘

	bool flag_ephopn;//ephopn函数是否已执行标识
	FILE * ephopn();//读取自己生成的文件jeleph.405，初始化数据

	bool flag_state;//state函数是否已执行标识
	void state(double *jed, int LList[], double pv[][13], double *nut);

	void RotMat(int axis, double phi, DMatrix r, DMatrix dr);
	void convert_little_endian(char *ptr, int len);
	void split(double tt, double *ipart, double *fpart);//将给定的数值分解为小数和整数部分

	bool flag_interp;//interp函数是否已执行标识
	void interp(int buff, double *t, int ncf, int ncm, int na, int fl, double dumpv[][2]);//区分和内插切比雪夫多项式来给定位置和速度矢量
	
	void errprt(int i, char const *msg);
	void freeDMatrix(DMatrix m, int n);//释放矩阵m的所有内存

};

visual_position_calculate::visual_position_calculate(int y=-1, int m=-1, int d=-1, int hour=-1, int minute=-1, int second=-1, int milisecond=-1, double tai_utcData=0.)
{
	if( y == -1)
	{
		win_time_val_t tv;
		win_time_t t;
		win_gettimeofday(&tv);
		win_time(&tv, &t);
		y = t.year;
		m = t.mon;
		d = t.day;
		hour = t.hour;
		minute = t.min;
		second = t.sec;
		milisecond = t.msec;
	}

	nYear = y;
	nMonth = m;
	nDay = d;
	nHour = hour;
	nMinute = minute;
	nSecond = (double)second + (double)milisecond/1000.0;	//时间精度为毫秒

	m_utc = (double)nHour + (double)nMinute/100.0 + nSecond/10000.0;  //hh.mmssss格式时间
	m_tai_utc = tai_utcData;
	timeReady = true;

	/////////////////////用于恒星视位置计算的部分变量初始化////////////////////////////

	flag_GetStateVector = true;
	flag_pleph = true;
	flag_GetDpsiDeps = true;
	flag_FunArgIAU = true;
	flag_ephopn = true;
	flag_state = true;
	flag_interp = true;

	double cbuf[15][3][8] = {0.};

	m_bePlant = false;
	m_beLunar = false;
	m_beAlt = false;
}


void visual_position_calculate::matrix_mxmul(double a[], double b[], double c[], int m, int n, int k)
{
	int i,j,r;
	for(i=0; i<m; i++)
		for(j=0; j<k; j++)
		{
			c[i*k+j] = 0.0;
			for(r=0; r<n; r++)
				c[i*k+j] += a[r+i*n] * b[r*k+j];
		}
}

void visual_position_calculate::matrix_mul(double a[], double b[], double c[], int m)
{
	int i;
	for(i=1; i<=m; i++)
		c[i]=a[i]*b[i];
}

double visual_position_calculate::matrix_sum(double a[], int m)
{
	int i;
	double x=0;
	for(i=1;i<=m;i++)
		x=x+a[i];
	return x;
}

int visual_position_calculate::matrix_inv(double a[], int n)
{
	int *is,*js,i,j,k,l,u,v;
	double d,p;
	is=new int[n];
	js=new int[n];
	for (k=0; k<=n-1; k++)
	{ d=0.0;
	for (i=k; i<=n-1; i++)
		for (j=k; j<=n-1; j++)
		{ l=i*n+j; p=fabs(a[l]);
	if (p>d) { d=p; is[k]=i; js[k]=j;}
	}
	if (d+1.0==1.0)
	{ delete(is); delete(js); printf("err**not inv\n");
	return(0);
	}
	if (is[k]!=k)
		for (j=0; j<=n-1; j++)
		{ u=k*n+j; v=is[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
	if (js[k]!=k)
		for (i=0; i<=n-1; i++)
		{ u=i*n+k; v=i*n+js[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
	l=k*n+k;
	a[l]=1.0/a[l];
	for (j=0; j<=n-1; j++)
		if (j!=k)
		{ u=k*n+j; a[u]=a[u]*a[l];}
		for (i=0; i<=n-1; i++)
			if (i!=k)
				for (j=0; j<=n-1; j++)
					if (j!=k)
					{ u=i*n+j;
		a[u]=a[u]-a[i*n+k]*a[k*n+j];
		}
		for (i=0; i<=n-1; i++)
			if (i!=k)
			{ u=i*n+k; a[u]=-a[u]*a[l];}
	}
	for (k=n-1; k>=0; k--)
	{ if (js[k]!=k)
	for (j=0; j<=n-1; j++)
	{ u=k*n+j; v=js[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
	if (is[k]!=k)
		for (i=0; i<=n-1; i++)
		{ u=i*n+k; v=i*n+is[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
	}
	delete(is); delete(js);
	return(1);
}
/*矩阵乘法函数, a[]中存放着矩阵A的元素，b[]中存放着B的元素，c[]存放AXB,m为A的行数，
n为A的列数，同时也是B的行数，k是B的列数*/


///////////////////////////////////////恒星视位置计算函数////////////////////////////////////////////////////////


double visual_position_calculate::Cal2JED(int m, int d, int y, double utc, double tai_utc)
{
	double time, datnum, a, b, term1, corr, T,jed;

	/* Convert to decimal hours */
	time = deg(utc);

	if ((m == 1) || (m == 2)) 
	{
		y = y - 1;
		m = m + 12;
	}

	/*
	Test to see if date is in Gregorian calendar.
	Converts date to a number of the form YYYY.MMDD.
	*/
	datnum = (double) y + (double) m/100.0 + (double) d/10000.0;

	if (datnum >= 1582.1015) 
	{
		/* Gregorian calendar */
		a = fix(0.01 * (double) y);
		b = 2.0 - a + fix(0.25 * (double) a);
	} 
	else 
	{
		/* Julian calendar */
		a = 0.0;
		b = 0.0;
	}

	if (y < 0) {
		/* Handles negative years */
		term1 = fix(365.25 * (double) y - 0.75); /* change here */
	} else {
		term1 = fix(365.25 * (double) y);
	}

	jed = term1 + fix(30.6001 * ((double) m + 1.0)) +
		(double) d + time / 24.0 + 1720994.5 + (double) b;


	/* First convert to TDT */
	corr = (tai_utc + 32.184) / 86400.0;
	jed = jed + corr;
	T = (jed - 2451545.0) /36525;
	/* Now compute the new correction in microseconds */
	corr = 1656.675     * sin((35999.3729 * T + 357.5387) * D2R);
	corr +=  22.418     * sin((32964.467  * T + 246.199)  * D2R);
	corr +=  13.84      * sin((71998.746  * T + 355.057)  * D2R);
	corr +=   4.77      * sin(( 3034.906  * T +  25.463)  * D2R);
	corr +=   4.67      * sin((34777.259  * T + 230.394)  * D2R);
	corr +=  10.216 * T * sin((35999.373  * T + 243.451)  * D2R);
	corr +=   0.171 * T * sin((71998.746  * T + 240.98 )  * D2R);
	corr +=   0.027 * T * sin(( 1222.114  * T + 194.661)  * D2R);
	corr +=   0.027 * T * sin(( 3034.906  * T + 336.061)  * D2R);
	corr +=   0.026 * T * sin((  -20.186  * T +   9.382)  * D2R);
	corr +=   0.007 * T * sin((29929.562  * T + 264.911)  * D2R);
	corr +=   0.006 * T * sin((  150.678  * T +  59.775)  * D2R);
	corr +=   0.005 * T * sin(( 9037.513  * T + 256.025)  * D2R);
	corr +=   0.043 * T * sin((35999.373  * T + 151.121)  * D2R);
	/* Convert from microseconds to seconds */
	corr = corr * 0.000001;
	/* Convert to days */
	corr = corr / 86400.0;


	jed += corr;

	/* Return modified JED if requested 
	if (w == 1) *jed -= 2400000.5;*/
	return jed;
}
//////////////////得到星体的状态矢量/////////////////////////////////////////////////////////////
bool visual_position_calculate::GetStateVector(double jd, int targ, int cent, int recpol, double StateVector[15][15][2][6])  // changed 
{
	double b[6] = {0.};
	double rrd[6] = {0.};//函数通过内插计算得到的天体相对于坐标原点的状态矢量（包括位置矢量和速度矢量）
	double dpsi, deps, dpsidot, depsdot;
	int inside, i;

	if(flag_GetStateVector)
	{
		dpsi=0.0;
		deps=0.0;
		dpsidot=0.0;
		depsdot=0.0;
		flag_GetStateVector=false;
	}

	/* Get the state vector from JPL ephemeris */
	if( !pleph(jd, targ, cent, rrd, &inside))//读取jpleph文件并给出星体相对于给定原点的位置和速度矢量
		return false;

	if (!inside) {

		for (i=0; i<6; i++) {
			StateVector[targ-1][cent-1][recpol-1][i] = 999.;
		}
		return false;
	}


	if (recpol == 1) {
		for (i=0; i<6; i++) {
			StateVector[targ-1][cent-1][recpol-1][i] = rrd[i];
		}
	} else {
		/* Convert rectangular vector to polar vector */
		Rec2Pol(rrd, b);
		for (i=0; i<6; i++) {
			StateVector[targ-1][cent-1][recpol-1][i] = b[i];
		}
	}

	return true;
}

/////////////////将h，m，s，ms转化成一个实数——弧度///////////////////////
double visual_position_calculate::hms2R(int a, int b, int c, int d)
{
	double r=((double)a*15.0+((double)b/4.0)+((double)c+(double)d*0.001)/240.0)*D2R;
	return r;
}

/////////////////////将星体由几何太阳系质心位置转化到几何位置或日心位置//////////////////////////////////////
void visual_position_calculate::LightTime(double jed, int body, double earth_ssb[], double earth_hel[], double StarData[], double body_geo[], double body_hel[], double *Ltime)
{
	int inside=0, i=0;
	double E[3],Edot[3],RQB[6],RSB[6],RQ[6],RP[6],Q[3],Qdot[3],P[3],Pdot[3];
	double unit_body_hel[6], unit_body_geo[6];
	double magE=0.0,magQ=0.0,magP=0.0,trialtau=0.0,tau=0.0,xtime=0.0,RelTerm=0.0,
		trueHelDist=0,trueGeoDist=0;
	double AppHelDist=0.0,AppGeoDist=0.0,RA=0.0,DE=0.0,sinRA=0.0,cosRA=0.0,sinDE=0.0,cosDE=0.0;
	double plx=0.0,r=0.0,muRA=0.0,muDE=0.0,rdot=0.0,drda[3],drdd[3],DT=0.0;

	memset(E,0,sizeof(E));
	memset(Edot,0,sizeof(Edot));
	memset(RQB,0,sizeof(RQB));
	memset(RSB,0,sizeof(RSB));
	memset(RQ,0,sizeof(RQ));
	memset(RP,0,sizeof(RP));
	memset(Q,0,sizeof(Q));
	memset(Qdot,0,sizeof(Qdot));
	memset(P,0,sizeof(P));
	memset(Pdot,0,sizeof(Pdot));
	memset(unit_body_hel,0,sizeof(unit_body_hel));
	memset(unit_body_geo,0,sizeof(unit_body_geo));
	memset(drda,0,sizeof(drda));
	memset(drdd,0,sizeof(drdd));

	SplitStateVector(earth_hel, E, Edot);
	magE = Vecmag(E);

	if (body != 99)
	{
		trialtau = 0.0;
		do {
			xtime = jed - trialtau;

			/* body's barycentric state vector */
			pleph(xtime, body, 12, RQB, &inside);

			/* Sun's barycentric state vector */
			pleph(xtime, 11, 12, RSB, &inside);

			/* body's heliocentric state vector */
			for (i=0; i<6; i++) {
				RQ[i] = RQB[i] - RSB[i];
				/* body's geocentric state vector */
				RP[i] = RQB[i] - earth_ssb[i];
			}

			/* Calculate improved value of the light time */
			SplitStateVector(RQ, Q, Qdot);
			SplitStateVector(RP, P, Pdot);
			magQ = Vecmag(Q);
			magP = Vecmag(P);

			if (trialtau == 0)
			{
				trueHelDist = magQ;
				trueGeoDist = magP;
			}

			if (magQ == 0) 
			{
				RelTerm = 0.0;
			} 
			else 
			{
				RelTerm = MUC * log((magE + magP + magQ) / fabs(magE - magP +
					magQ));
			}
			tau = (magP + RelTerm) / CAUD;
			if (fabs(trialtau - tau) < 1e-9) 
			{
				break;
			}
			trialtau = tau;
		} while(1);
	} 
	else
	{

		RA = StarData[0];
		DE = StarData[1];
		sinRA = sin(RA);
		cosRA = cos(RA);
		sinDE = sin(DE);
		cosDE = cos(DE);
		if (StarData[2] == 0)
		{
			plx = 1e-7;
		} 
		else 
		{
			plx = StarData[2];
		}

		/* Star's distance in AU */
		r = R2A / plx;

		/* Star's barycentric position vector */
		RQB[0] = r*cosRA*cosDE;
		RQB[1] = r*sinRA*cosDE;
		RQB[2] = r      *sinDE;

		/* Convert proper motions and radial velocity */
		/* to units of AU/day */
		muRA = StarData[3] * 15.0 * cosDE / (JulCty * plx);
		muDE = StarData[4] / (JulCty * plx);
		rdot = StarData[5] * 86400.0 * KM2AU;

		/* Partial of unit vector wrt RA */
		drda[0] = - sinRA * cosDE;
		drda[1] =   cosRA * cosDE;
		drda[2] =   0.0;

		/* Partial of unit vector wrt DE */
		drdd[0] = - cosRA * sinDE;
		drdd[1] = - sinRA * sinDE;
		drdd[2] =           cosDE;

		/* Star's barycentric velocity vector */
		for (i=3; i<6; i++) {
			RQB[i] = muRA*drda[i-3] + muDE*drdd[i-3] + rdot*RQB[i-3]/r;
		}

		/* Correction for space motion */
		DT = jed - J2000;
		for (i=0; i<3; i++) {
			RQB[i] = RQB[i] + RQB[i+3]*DT;
		}

		/* Sun's barycentric state vector */
		pleph(jed, 11, 12, RSB, &inside);

		/* Star's heliocentric state vector */
		for (i=0; i<6; i++) {
			RQ[i] = RQB[i] - RSB[i];
		}

		/* Star's geo/topocentric state vector */
		/* This correction is annual parallax. */
		for (i=0; i<6; i++) {
			RP[i] = RQB[i] - earth_ssb[i];
		}

		SplitStateVector(RQ, Q, Qdot);
		SplitStateVector(RP, P, Pdot);
		magQ = Vecmag(Q);
		magP = Vecmag(P);

		trueGeoDist = magP;
		trueHelDist = magQ;

		tau = r / CAUD;
	}

	*Ltime = tau;
	AppHelDist = magQ;
	AppGeoDist = magP;

	for (i=0; i<6; i++) {
		body_hel[i] = RQ[i];
		body_geo[i] = RP[i];
	}

	Uvector(body_hel, unit_body_hel);
	Uvector(body_geo, unit_body_geo);

	/* I don't understand why this next loop is needed, but it is. */
	for (i=3; i<6; i++) {
		unit_body_hel[i] = RQ[i];
		unit_body_geo[i] = RP[i];
	}

	for (i=0; i<6; i++) {
		body_geo[i] = unit_body_geo[i] * trueGeoDist;
		body_hel[i] = unit_body_hel[i] * trueHelDist;
	}
}

/////////////根据太阳的重力场来计算相对光偏差///////////////////////////////
void visual_position_calculate::RayBend(double earth_hel[], double body_geo[], double body_hel[], double p1[])
{
	double ue[3], up[3], uq[3], E[3], Edot[3], P[3], Pdot[3], Q[3], Qdot[3];
	double magE, magP, pdotq, edotp, qdote;
	int i;

	/* extract the pos. portions of the state vectors */
	SplitStateVector(earth_hel, E, Edot);
	SplitStateVector( body_geo, P, Pdot);
	SplitStateVector( body_hel, Q, Qdot);

	/* form unit vectors */
	Uvector(E, ue);
	Uvector(P, up);
	Uvector(Q, uq);

	/* form dot products and other quantities */
	magE = Vecmag(E);
	magP = Vecmag(P);
	Vdot(3, up, uq, &pdotq);
	Vdot(3, ue, up, &edotp);
	Vdot(3, uq, ue, &qdote);

	for (i=0; i<3; i++) {
		p1[i] = up[i] + (MUC / magE) * (pdotq * ue[i] - edotp * uq[i])
			/ (1.0 + qdote);
		/* make p1[] a non-unit vector */
		p1[i] = magP * p1[i];
		p1[i+3] = body_geo[i+3];
	}
}

void visual_position_calculate::SplitStateVector(double pv[], double p[], double v[])
{
	int i;

	for (i=0; i<3; i++) {
		p[i] = pv[i];
		v[i] = pv[i+3];
	}
}
////////////////////计算光行差////////////////////////////
void visual_position_calculate::Aberrate(double p1[], double EBdot[], double p2[])
{
	double v[3], up[3], p[3], pdot[3], magP, magV, p1dotv, beta;
	int i;

	/* Extract the pos. portion of the state vector */
	SplitStateVector(p1, p, pdot);
	magP = Vecmag(p);

	/* Need to make Ppos() a unit vector */
	Uvector(p, up);

	for (i=0; i<3; i++) {
		v[i] = EBdot[i] / CAUD;
	}

	Vdot(3, up, v, &p1dotv);
	magV = Vecmag(v);

	beta = 1.0 / sqrt(1.0 - magV * magV);

	for (i=0; i<3; i++) {
		p2[i] = ((up[i] / beta) +
			(1.0 + p1dotv / (1.0 + (1.0 / beta))) * v[i]) / (1.0 + p1dotv);
		/* Make p2[] a non-unit vector */
		p2[i] = magP * p2[i];
		p2[i+3] = p1[i+3];
	}

	return;
}
///////////////////////////计算岁差参数////////////////////////////
void visual_position_calculate::GetPrecessParams(double jed1, double jed2, double *zeta, double *z, double *theta, double *zetadot, double *zdot, double *thetadot)
{
	double T1, T2, c1, c2, c3, c4, c5, c6, p1, p2, x, xdot;

	T1 = (jed1 - J2000) / JulCty;
	T2 = (jed2 - jed1) / JulCty;

	/* compute zeta, z, theta, zetadot, zdot, thetadot */
	c1 = 2306.2181;
	c2 =    1.39656;
	c3 =   -0.000139;
	c4 =    0.30188;
	c5 =   -0.000344;
	c6 =    0.017998;
	p1 = c1 + c2 * T1 + c3 * T1 * T1;
	p2 = c4 + c5 * T1;
	x = p1 * T2 + p2 * T2 * T2 + c6 * T2 * T2 * T2;
	xdot = p1 + 2.0 * p2 * T2 + 3.0 * c6 * T2 * T2;
	*zeta = x * A2R;
	*zetadot = xdot * A2R / JulCty;

	c1 = 2306.2181;
	c2 =    1.39656;
	c3 =   -0.000139;
	c4 =    1.09468;
	c5 =    0.000066;
	c6 =    0.018203;
	p1 = c1 + c2 * T1 + c3 * T1 * T1;
	p2 = c4 + c5 * T1;
	x = p1 * T2 + p2 * T2 * T2 + c6 * T2 * T2 * T2;
	xdot = p1 + 2.0 * p2 * T2 + 3.0 * c6 * T2 * T2;
	*z = x * A2R;
	*zdot = xdot * A2R / JulCty;

	c1 = 2004.3109;
	c2 =   -0.85330;
	c3 =   -0.000217;
	c4 =   -0.42665;
	c5 =   -0.000217;
	c6 =   -0.041833;
	p1 = c1 + c2 * T1 + c3 * T1 * T1;
	p2 = c4 + c5 * T1;
	x = p1 * T2 + p2 * T2 * T2 + c6 * T2 * T2 * T2;
	xdot = p1 + 2.0 * p2 * T2 + 3.0 * c6 * T2 * T2;
	*theta = x * A2R;
	*thetadot = xdot * A2R / JulCty;
}
///////////////////////计算章动dpsi, deps, dpsidot, and depsdot.////////////////////////////
void visual_position_calculate::GetDpsiDeps(double jed, double *dpsi, double *deps, double *dpsidot, double *depsdot)
{
	double L, Ldot, LP, LPdot, F, Fdot, D, Ddot, N, Ndot, LL, LLdot;
	double T, c1, c2, c3, c4, c5, lamp, oamp, ls, os, arg, argdot;
	double funarg[12];
	int j;
	double dpsiz, depsz, dpsidotz, depsdotz, jedz;

	if(flag_GetDpsiDeps)
	{
		dpsiz=0.0;
		depsz=0.0;
		dpsidotz=0.0;
		depsdotz=0.0;
		jedz=0.0;
		flag_GetDpsiDeps=false;
	}

	

	if (jedz != jed) {
		jedz = jed;

		T = (jed - J2000) / JulCty;

		//计算章动的五个基本变量函数月球的平近点角，太阳的平近点角，月球平升交角距，日月平角距，月球轨道对黄道平均升交点的黄经。
		FunArgIAU(jed, funarg);
		L = funarg[0];
		Ldot = funarg[6];
		LP = funarg[1];
		LPdot = funarg[7];
		F  = funarg[2];
		Fdot  = funarg[8];
		D  = funarg[3];
		Ddot  = funarg[9];
		N  = funarg[4];
		Ndot  = funarg[10];
		LL = funarg[5];
		LLdot = funarg[11];

		/* evaluate the series */
		dpsiz = 0.;  /* initialize to zero */
		depsz = 0.;
		dpsidotz = 0.;
		depsdotz = 0.;
		for (j=0; j<106; j++) {
			lamp = IAU1980[j][0];
			oamp = IAU1980[j][2];
			ls   = IAU1980[j][1];
			os   = IAU1980[j][3];
			c1   = IAU1980[j][4];
			c2   = IAU1980[j][5];
			c3   = IAU1980[j][6];
			c4   = IAU1980[j][7];
			c5   = IAU1980[j][8];
			arg  = c1 * L + c2 * LP + c3 * F + c4 * D + c5 * N;
			arg  = amodulo(arg, TWOPI);
			dpsiz += (lamp + ls * T) * sin(arg);
			depsz += (oamp + os * T) * cos(arg);
			argdot = c1 * Ldot + c2 * LPdot + c3 * Fdot + c4 * Ddot + c5 * Ndot;
			argdot = amodulo(argdot, TWOPI);
			dpsidotz += (lamp + ls * T) * argdot * cos(arg) +
				ls * sin(arg) / JulCty;
			depsdotz -= (oamp + os * T) * argdot * sin(arg) +
				os * cos(arg) / JulCty;
		}

		/* normalize and convert units */
		dpsiz *= 0.0001 * A2R;
		depsz *= 0.0001 * A2R;
		dpsidotz *= 0.0001 * A2R;
		depsdotz *= 0.0001 * A2R;
	}

	*dpsi = dpsiz;
	*deps = depsz;
	*dpsidot = dpsidotz;
	*depsdot = depsdotz;
}
//////////////////////////////计算平黄赤交角及它的导数//////////////////////////////
void visual_position_calculate::Obliquity(double jed1, double jed2, int m, double *obl, double *obldot)
{
	double t1, t2, e0, e1, e2, e3, e4, e5, e6, epsbar;
	double dpsi, deps, dpsidot, depsdot;

	t1 = (jed1 - J2000) / JulCty;
	t2 = (jed2 - jed1) / JulCty;

	e0 = 84381.448;
	e1 =   -46.815;
	e2 =    -0.00059;
	e3 =     0.001813;
	epsbar = e0 + e1 * t1 + e2 * t1 * t1 + e3 * t1 * t1 * t1;
	e1 = -46.815;
	e2 =  -0.00117;
	e3 =   0.005439;
	e4 =  -0.00059;
	e5 =   0.005439;
	e6 =   0.001813;
	*obl = epsbar + (e1 + t1 * (e2 + t1 * e3)) * t2
		+ (e4 + e5 * t1) * t2 * t2
		+ e6 * t2 * t2 * t2;
	*obldot = e1 + t1 * (e2 + t1 * e3)
		+ 2.0 * (e4 + e5 * t1) * t2
		+ 3.0 * e6 * t2 * t2;

	if (m == 1) {
		/* need true obliquity */
		GetDpsiDeps(jed2, &dpsi, &deps, &dpsidot, &depsdot);
		/* Unit conversion is needed because obl is  */
		/* in arc seconds, nutations are in radians. */
		*obl = *obl + deps * R2A;
		*obldot = *obldot + depsdot * R2A;
	}

	/* convert to radians and radians/day */
	*obl = *obl * A2R;
	*obldot = *obldot * A2R / JulCty;
}

////////////通过给定的角度来计算矩阵的旋转后的坐标矢量（右手法则）//////////////////////////////
//这个函数的功能基本上与Rotmat函数的功能一致，只不过它是6×6矩阵，并且遵循右手法则。
void visual_position_calculate::QRotate(double vin[], int axis, double phi, double phidot, int s, double vout[])
{
	int i;
	double temp[6];
	DMatrix QMatrix;

	QMatrix = createDMatrix(6);//创建一个n*n矩阵//

	GetQMatrix(phi, phidot, axis, s, QMatrix);//对于给定的角度和时间倒数来计算矩阵的旋转//

	for (i=0; i<6; i++) {
		temp[i] = vin[i];
	}

	MatXVec(QMatrix, temp, vout, 6, 6);//用来实现矩阵的乘法运算//

	freeDMatrix(QMatrix, 6);//释放矩阵m的所有内存//
}
///////////////将 hh.mmssss 时间转化到 hh.hhhhhh 时间///////////
double visual_position_calculate::deg(double x)
{
	double hh, h, fixhh, hhfixhh;

	if (x == 0.0) return (0.0);

	hh = fabs(x);
	fixhh = floor(hh);
	hhfixhh = hh - fixhh + 5.0e-10;  /* fudge factor */
	h = fixhh + floor(100.0 * hhfixhh) / 60.0;
	h = h + (10000.0 * hhfixhh - 100.0 * floor(100.0 * hhfixhh)) / 3600.0;

	return ((x / hh) * h);

}
///////////返回到实数的整数值////////////////////////
double visual_position_calculate::fix(double x)
{
	if (x < 0)
	{
		return ceil(x);
	} 
	else
	{
		return floor(x);
	}
}
//////////////////////得到当前应用程序所在绝对路径///////////////////
bool visual_position_calculate::GetAppPath(char *Path)
{
	char tmpStr[MAX_PATH] ;
	int i, j ;

	if (GetModuleFileName (NULL, (TCHAR*)tmpStr, MAX_PATH) == 0)
		return false ;
	_strrev (tmpStr) ;

	for (i =0 ; i < strlen(tmpStr); i++)
		if (tmpStr[i] == '\\')
			break ;

	_strrev (tmpStr) ;
	j = strlen(tmpStr) - i-1 ;
	tmpStr[j] = '\0' ;
	strcpy(Path, tmpStr) ;

	return true ;
}
////////////////读取JPL行星星表并给出星体相对于给定原点的位置和速度矢量//////////////
bool visual_position_calculate::pleph(double jd, int targ, int cent, double *rrd, int *inside)
{
	int i;
	double fac;
	int nemb, ipv, ncmp, lme;
	int pfirst;
	double jdtot;
	double pv[6][13]={0.};
	double embf[2], ve[2], jed[2];
	int LList[12], LLst[13];
	int L[2], tc[2];

	if(flag_pleph)
	{
		fac=0.0;
		nemb=0;
		ipv=0;
		ncmp=0;
		lme=0;
		pfirst=0;
		jdtot=0.0;
		memset(embf,0,sizeof(embf));
		memset(ve,0,sizeof(ve));
		memset(jed,0,sizeof(jed));
		memset(LList,0,sizeof(LList));
		memset(LLst,0,sizeof(LLst));
		memset(L,0,sizeof(L));
		memset(tc,0,sizeof(tc));
		flag_pleph=false;
	}

	/* necessary for zero-offset arrays */
	targ = targ - 1;
	cent = cent - 1;


	/* 1st time in, be sure ephemeris is initialized */
	if (!pfirst) {
		pfirst = true;
		ipv = 2; /* we want pos and vel */
		if ( ephopn() == NULL )
		{
			pfirst = false ;
			return false ;
		}
		ve[0] = 1.0/(1.0+emrat);
		ve[1] = emrat*ve[0];

		jed[0] = 0.0;
		jed[1] = 0.0;

		embf[0] = -1.0;
		embf[1] =  1.0;
		for (i=0; i<12; i++) {
			LList[i] = 0;
		}
		L[0] = 0;
		L[1] = 0;
		tc[0] = 0;
		tc[1] = 0;
		for (i=0; i<13; i++) {
			LLst[i] = i;
			if (i ==  2) LLst[i] = 9;
			if (i == 11) LLst[i] = 10;
			if (i == 12) LLst[i] = 2;
		}
		fac = 0.0;
		nemb = 1;
	}

	/* Initialize jed[] for state() and set up component count */
	jed[0] = jd;
	jed[1] = 0.0;

	jdtot = jed[0] + jed[1];

	if ((jdtot >= SS[0]) && (jdtot <= SS[1])) {
		*inside = true;
	} else {
		*inside = false;
		return false;
	}

	ncmp = 3*ipv; /* total number of components */
	if (targ == 13) {
		if (ipt[1][11] > 0) {
			LList[10] = ipv;
			state(jed, LList, pv, rrd);
			LList[10] = 0;
			return false;
		} 
		else 
		{
			//     LogMsg(stderr, "pleph: no nutations on the ephemeris file.");
			exit(1);
		}
	}	
	if (targ == 14) {
		if (lpt[1] > 0) {
			LList[11] = ipv;
			state(jed, LList, pv, rrd);
			LList[11] = 0;
			for (i=0; i<ncmp; i++) {
				rrd[i] = pv[i][10];
			}
			return false;
		} else {
			//     LogMsg(stderr, "pleph: no librations on the ephemeris file.");
			exit(1);
		}
	}

	/* check for targ = cent */
	if (targ == cent) {
		for (i=0; i<ncmp; i++) {
			rrd[i] = 0.0;
		}
		return false;
	}
	/* force barycentric output by state() */
	bsav = bary;
	bary = true;

	/* set up proper entries in LList[] array for state() call */
	tc[0] = targ;
	tc[1] = cent;
	lme = 0;

	//问题出在数组的下标不能为-1，因为传进来的参数是0，所以targ必须是-1了，要作修改
	for (i=0; i<2; i++) {
		L[i] = LLst[tc[i]];
		if (L[i] < 10) LList[L[i]] = ipv;
		if (tc[i] == 2) {
			lme = 2;
			fac = -ve[0];
		}
		else if (tc[i] == 9) {
			lme = 9;
			fac = ve[1];
		}
		else if (tc[i] == 12) {
			nemb = i;
		}
	}

	if ((LList[9] == ipv) && (L[0] != L[1])) LList[2] = ipv - LList[2];


	state(jed, LList, pv, rrd);

	for (i=0; i<ncmp; i++) {
		pv[i][10] = pvsun[i];
		pv[i][12] = pv[i][2];
		if (lme > 0) pv[i][lme] = pv[i][2]+fac*pv[i][9];
		rrd[i] = pv[i][targ] - pv[i][cent];
	}	
	/* clear state() body array and restore barycenter flag */
	LList[2] = 0;
	LList[L[0]] = 0;
	LList[L[1]] = 0;
	bary = bsav;
	return true;
}
///////////////将直角坐标转变为极坐标/////////////////////
void visual_position_calculate::Rec2Pol(double a[], double b[])
{double x, y, z, x_dot, y_dot, z_dot,
rho, r, lambda, beta, lambda_dot,
beta_dot, r_dot;

x = a[0];
y = a[1];
z = a[2];
x_dot = a[3];
y_dot = a[4];
z_dot = a[5];

rho = sqrt(x * x + y * y);
r = sqrt(rho * rho + z * z);

lambda = atan2(y, x);
if (lambda < 0.0)
lambda += TWOPI;

beta   = atan2(z, rho);
if (beta < 0.0)
beta += TWOPI;

if (z < 0) {
	beta = beta - TWOPI;
}

if (rho == 0) {
	lambda_dot = 0.0;
	beta_dot = 0.0;
} else {
	lambda_dot = (x*y_dot-y*x_dot) / (rho*rho);
	beta_dot = (z_dot*rho*rho-z*(x*x_dot+y*y_dot))/(r*r*rho);
}

r_dot = (x * x_dot + y * y_dot + z * z_dot) / r;

/* position vector components */
b[0] = lambda;
if (b[0] >= TWOPI)
b[0] = b[0] - TWOPI;

b[1] = beta;
b[2] = r;
/* total velocity vector components */
b[3] = r * lambda_dot * cos(beta);
b[4] = r * beta_dot;
b[5] = r_dot;

}
////////////计算矢量的数值///////////////
double visual_position_calculate::Vecmag(double a[])
{
	double x;

	x = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

	return (x);
}
/////////////////////计算单位矢量/////////////////
void visual_position_calculate::Uvector(double a[], double unita[])
{
	double maga;
	int i;

	maga = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	for (i=0; i<3; i++) {
		if (maga != 0) {
			unita[i] = a[i]/maga;
		} else {
			unita[i] = 0.0;
		}
	}
}
//////////////////////矢量乘法////////////////////////////////
void visual_position_calculate::Vdot(int n, double a[], double b[], double *adotb)
{
	int i;

	*adotb = 0.0;
	for (i = 0; i<n; i++) {
		*adotb += a[i] * b[i];
	}
}
//////////////利用IAU表达式计算基本参数///计算章动的五个基本变量函数
void visual_position_calculate::FunArgIAU(double jed, double *funarg)
{
	double T, T2, T3, L, Ldot, lp, lpdot, F, Fdot, D, Ddot;
	double N, Ndot, LL, LLdot;

	if(flag_FunArgIAU)
	{
		jedz=0.0;
		flag_FunArgIAU=false;
	}

	if (jed == jedz) return;
	jedz = jed;

	T = (jed - J2000) / JulCty;
	T2 = T * T;
	T3 = T2 * T;

	/* Compute fundamental arguments */
	L = 485866.733 + (1325. * 1296000. + 715922.633) * T
		+ 31.31 * T2 + 0.064 * T3;
	Ldot = (1325. * 1296000. + 715922.633) + 2. * 31.31 * T
		+ 3. * 0.064 * T2;
	L = amodulo(L * A2R, TWOPI);
	Ldot = amodulo(Ldot * A2R / 36525., TWOPI);
	lp = 1287099.804 + (99. * 1296000. + 1292581.224) * T
		- 0.577 * T2 - 0.012 * T3;
	lpdot = (99. * 1296000. + 1292581.224) - 2. * 0.577 * T
		- 3. * 0.012 * T2;
	lp = amodulo(lp * A2R, TWOPI);
	lpdot = amodulo(lpdot * A2R / 36525., TWOPI);
	F = 335778.877 + (1342. * 1296000. + 295263.137) * T
		- 13.257 * T2 + 0.011 * T3;
	Fdot = (1342. * 1296000. + 295263.137) - 2. * 13.257 * T
		+ 3. * 0.011 * T2;
	F = amodulo(F * A2R, TWOPI);
	Fdot = amodulo(Fdot * A2R / 36525., TWOPI);
	D = 1072261.307 + (1236. * 1296000. + 1105601.328) * T
		- 6.891 * T2 + 0.019 * T3;
	Ddot = (1236. * 1296000. + 1105601.328) - 2. * 6.891 * T
		+ 3. * 0.019 * T2;
	D = amodulo(D * A2R, TWOPI);
	Ddot = amodulo(Ddot * A2R / 36525., TWOPI);
	N = 450160.28 - (5. * 1296000. + 482890.539) * T
		+ 7.455 * T2 + 0.008 * T3;
	Ndot = (5. * 1296000. + 482890.539) + 2. * 7.455 * T
		+ 3. * 0.008 * T2;
	N = amodulo(N * A2R, TWOPI);
	Ndot = amodulo(Ndot * A2R / 36525., TWOPI);
	LL = 785939.157 + (1336. * 1296000. + 1108372.598) * T
		- 5.802 * T2 + 0.019 * T3;
	LLdot = (1336. * 1296000. + 1108372.598) - 2. * 5.802 * T
		+ 3. * 0.019 * T2;
	LL = amodulo(LL * A2R, TWOPI);
	LLdot = amodulo(LLdot * A2R / 36525., TWOPI);

	funarg[0]  = L;
	funarg[6]  = Ldot;
	funarg[1]  = lp;
	funarg[7]  = lpdot;
	funarg[2]  = F;
	funarg[8]  = Fdot;
	funarg[3]  = D;
	funarg[9]  = Ddot;
	funarg[4]  = N;
	funarg[10] = Ndot;
	funarg[5]  = LL;
	funarg[11] = LLdot;
}
/////////计算一个在0 <= a < b范围内的a值///////////
double visual_position_calculate::amodulo(double a, double b)
{
	double x;

	x = a - b * floor(a/b);
	return (x);
}
////////////////创建一个n*n矩阵//////////////////
DMatrix visual_position_calculate::createDMatrix(int n)
{
	DMatrix m;
	int i;

	m = (double**)calloc(n, sizeof(double*));
	if (m == NULL) {
		return NULL;
	}

	for (i=0; i<n; i++) {
		m[i] = (double*)calloc(n, sizeof(double));
		if (m[i] == NULL) {
			freeDMatrix(m, i); /* Avoids garbage */
			return NULL;
		}
	}

	return m;
}
////////////对于给定的角度和时间倒数来计算矩阵的旋转/////////////////////////////
void visual_position_calculate::GetQMatrix(double phi, double phidot, int axis, int s, DMatrix QMatrix)
{
	DMatrix r, dr;
	int i, j;

	r  = createDMatrix(3);
	dr = createDMatrix(3);

	/* form the 3X3 r[] and dr[] matrices */
	RotMat(axis, phi, r, dr);

	/* form the 6X6 Q-matrix */
	for (i=0; i<3; i++) 
	{
		for (j=0; j<3; j++) 
		{
			QMatrix[i][j] = r[i][j];
			QMatrix[i][j+3] = 0.0;
			QMatrix[i+3][j] = phidot * dr[i][j] * s;
			QMatrix[i+3][j+3] = r[i][j];
		}
	}

	freeDMatrix(r, 3);
	freeDMatrix(dr, 3);
}
//////////////////矩阵相乘//////////////用来实现矩阵的乘法运算/////
void visual_position_calculate::MatXVec(DMatrix a, double b[], double c[], int l, int m)
{
	int i, j;
	double s, temp[6];

	for (i=0; i<l; i++) {
		s = 0.0;
		for (j=0; j<m; j++) {
			s += a[i][j] * b[j];
		}
		temp[i] = s;
	}

	for (i=0; i<l; i++) {
		c[i] = temp[i];
	}
}
////////////////读取自己生成的文件jeleph.405，初始化数据/////////////////
FILE * visual_position_calculate::ephopn() 
{
	/* Make sure all variables read from the data file are global */
	char FileName[256];
	fpBinaryFile = NULL ;
	memset(FileName,0,sizeof(FileName));
	int i, j;
	long curpos;

	if(flag_ephopn)
	{
		efirst=0;
		flag_ephopn=false;
	}

	if (!efirst) {

		/* Default file name is JPLEPH */
		if (strlen(FileName) == 0) 
		{
			char path[256] ;
			GetAppPath(path) ;
			sprintf_s(FileName,"%s\\jpleph.405",path);
		}

		fpBinaryFile = fopen(FileName, "rb") ;
		if (fpBinaryFile == NULL)
		{
			return (NULL);
		}

		for (i = 0; i < 3; i++) {
			if (fread(ttl[i], sizeof(char), 65, fpBinaryFile) != 65) {
				return (NULL);
			}
			ttl[i][64] = '\0';
		}
		/* 2 ints */
		//读取两个字节，可是文件中没有
		if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
			return (NULL);
		}

		//    convert_little_endian((char *) &tmpShort, sizeof(short));
		ncon = (int) tmpShort;

		/* ncon*6 ints */
		for (j = 0; j < ncon; ++j) {
			/* char CNAM[NCON][7] */
			if (fread(&cnam[j], sizeof(char), 6, fpBinaryFile) != 6) {
				return (NULL);
			}
			cnam[j][6] = '\0';
		}

		/* 24 ints, 8 each */
		for (j = 0; j < 3; j++) 
		{
			/* double SS[3] */
			if (fread(&tmpDouble, sizeof(double), 1, fpBinaryFile) != 1) {
				return (NULL);
			}
			convert_little_endian((char *) &tmpDouble, sizeof(double));
			SS[j] = (double) tmpDouble;
		}

		/* 16 ints, 8 each */
		if (fread(&tmpDouble,    sizeof(double), 1, fpBinaryFile) != 1) {
			return (NULL);
		}
		convert_little_endian((char *) &tmpDouble, sizeof(double));
		au = (double) tmpDouble;

		if (fread(&tmpDouble, sizeof(double), 1, fpBinaryFile) != 1) {
			return (NULL);
		}
		convert_little_endian((char *) &tmpDouble, sizeof(double));
		emrat = (double) tmpDouble;

		/* 72 ints */
		for (i = 0; i < 3; i++) {
			/* short IPT[3][12] */
			for (j = 0; j < 12; j++) {
				if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
					return (NULL);
				}
				convert_little_endian((char *) &tmpShort, sizeof(short));
				ipt[i][j] = (short) tmpShort;
			}
		}

		/* 2 ints */
		if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
			return (NULL);
		}
		convert_little_endian((char *) &tmpShort, sizeof(short));
		NUMDE = (short) tmpShort;

		/* 6 ints, 2 each */
		for (i = 0; i < 3; i++) {
			if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
				return (NULL);
			}
			convert_little_endian((char *) &tmpShort, sizeof(short));
			lpt[i] = (short) tmpShort;
		}

		/* ncon*8 ints */
		for (j = 0; j < ncon; j++) {
			if (fread(&tmpDouble, sizeof(double), 1, fpBinaryFile) != 1) {
				return (NULL);
			}
			convert_little_endian((char *) &tmpDouble, sizeof(double));
			cval[j] = (double) tmpDouble;
		}
		/* END OF HEADER */

		if ((curpos = ftell(fpBinaryFile)) == -1L) {
			return (NULL);
		}

		if (fseek(fpBinaryFile, 0L, SEEK_END) != 0) {
			return (NULL);
		}

		LengthOfFile = ftell(fpBinaryFile);

		if (fseek(fpBinaryFile, curpos, SEEK_SET) != 0) {
			return (NULL);
		}
		/* end block */

		LengthOfHeader = 317L + 14L * (long) ncon;

		ncoeff = 0L;

		for (j = 0; j < 12; j++) {
			if (j == 11)
				ncoeff = ncoeff + (long) ipt[1][j] * (long) ipt[2][j] * 2L;
			else
				ncoeff = ncoeff + (long) ipt[1][j] * (long) ipt[2][j] * 3L;
		}

		ncoeff = ncoeff + (long) lpt[1] * (long) lpt[2] * 3L + 2L;
		BlockLength = ncoeff * 8L;
		NumBlocks = (LengthOfFile - (long) LengthOfHeader) / BlockLength;
		efirst = true;
	}

	return (fpBinaryFile);
}
/////////////读取和内插JPL星表文件///////////////////
void visual_position_calculate::state(double *jed, int LList[], double pv[][13], double *nut)
{
	double t[2], jd[4], temp[6];
	double dumpv[3][2]={{0.,0.},{0.,0.},{0.,0.}};
	int sfirst;
	int nrl;

	int buff, ncf, na, km = 0;
	double aufac, s, ipart, fpart;

	int i, j, k, m;
	int nr;

	if(flag_state)
	{
		memset(t,0,sizeof(t));
		memset(jd,0,sizeof(jd));
		memset(temp,0,sizeof(temp));
		sfirst=0;
		nrl=0;
		buff=0;
		ncf=0;
		na=0;
		aufac=0.0;
		s=0.0;
		ipart=0.0;
		fpart=0.0;
		flag_state=false;
	}

	/* 1st time in, get pointer data, etc., from ephemeris file */
	if (!sfirst) {
		sfirst = true;
		aufac = 1.0;
		nrl = 0L;
		//		ephopn("");
		if (km) {
			t[1] = SS[2]*86400.0;
		} else {
			t[1] = SS[2];
			aufac = 1.0 / au;
		}
	}

	/* main entry point -- check epoch and read right record */
	s = jed[0] - 0.5;
	split(s, &ipart, &fpart);
	jd[0] = ipart;
	jd[1] = fpart;
	split(jed[1], &ipart, &fpart);
	jd[2] = ipart;
	jd[3] = fpart;
	jd[0] = jd[0] + jd[2] + 0.5;
	jd[1] = jd[1] + jd[3];
	split(jd[1], &ipart, &fpart);
	jd[2] = ipart;
	jd[3] = fpart;
	jd[0] = jd[0] + jd[2];

	/* error return of epoch out of range */
	if ((jd[0] < SS[0]) || (jd[0]+jd[3]) > SS[1]) {
		//    LogMsg(stderr, "state: epoch out of range.\n");
		exit(1);
	}

	/* 'nr' is the int index of the first coefficient */
	nr = (long) (floor((jd[0]-SS[0])/SS[2]));
	nr = LengthOfHeader + nr * BlockLength;
	/* use previous block if necessary */
	if (jd[0] == SS[1])
		nr = nr - BlockLength;
	if (nr < 1L) {
		//   LogMsg(stderr, "state: block not present.\n");
		exit(1);
	}

	/* calculate relative time in interval (0 <= t[0] <= 1) */
	t[0] = ((jd[0]-(((double)nr - (double)(LengthOfHeader))/(double)BlockLength
		* SS[2] + SS[0])) + jd[3]) / SS[2];
	if (t[0] < 0.0)
		t[0] = t[0]+ 1.0;

	/* read correct record if not in core */
	if (nr != nrl) {
		nrl = nr;
		if (fseek(fpBinaryFile, (long) (nr), 0) != 0) {
			//   LogMsg(stderr, "state: fseek() failed.\n");
			exit(1);
		}
		k = 0;
		do {
			if (k < ncoeff) {
				fread(&tmpDouble, sizeof(double), 1, fpBinaryFile);
				//				convert_little_endian((char *) &tmpDouble, sizeof(double));
				db[k] = (double) tmpDouble;
			}
			k++;
		} while (!feof(fpBinaryFile) && k <= ncoeff);
	}

	/* interpolate SSBARY Sun */
	buff = ipt[0][10]-1; /* location of first coeff */
	ncf  = ipt[1][10]; /* number of coeffs per component */
	na   = ipt[2][10]; /* number of sets of coeffs per 32day interval */

	interp(buff, t, ncf, 3, na, 2, dumpv);

	k = 0;
	for (j=0; j<2; j++) {
		for (i=0; i<3; i++) {
			pvsun[k] = dumpv[i][j]*aufac;
			k++;
		}
	}

	/* check and interpolate whichever bodies are requested */
	for (i=0; i<10; i++) {
		if (LList[i] <= 0)
			continue;
		if (ipt[1][i] <= 0) {
			errprt(i ,"th body requested - not on file.\n");
		}
		buff = ipt[0][i]-1; /* location of first coeff */
		ncf  = ipt[1][i]; /* number of coeffs per component */
		na   = ipt[2][i]; /* number of sets of coeffs per 32day interval */

		interp(buff, t, ncf, 3, na, LList[i], dumpv);

		/* need to re-map dumpv[1..3][1..2] --> temp[1..6] */
		k = 0;
		for (j=0; j<2; j++) {
			for (m=0; m<3; m++) {
				temp[k] = dumpv[m][j];
				k++;
			}
		}

		for (j=0; j<(LList[i]*3); j++) {
			if ((i <= 8) && (!bary)) {
				pv[j][i] = temp[j]*aufac-pvsun[j];
			} else {
				pv[j][i] = temp[j]*aufac;
			}
		}
	}

	/* do nutations if requested and if on file */
	if ((LList[10] > 0) && (ipt[1][11] > 0)) {
		buff = ipt[0][11]-1; /* location of first coeff */
		ncf  = ipt[1][11]; /* number of coeffs per component */
		na   = ipt[2][11]; /* number of sets of coeffs per 32day interval */

		interp(buff, t, ncf, 2, na, LList[10], dumpv);

		/* need to re-map dumpv(1:3,1:2) --> temp(1:6) */
		k = 0;
		for (j=0; j<2; j++) {
			for (m=0; m<2; m++) {
				nut[k] = dumpv[m][j];
				k++;
			}
		}
		nut[4] = 0.0;
		nut[5] = 0.0;
	}

	/* get librations if requested and if on file */
	if ((lpt[1] > 0) && (LList[11] > 0)) {
		buff = lpt[0]-1; /* location of first coeff */
		ncf  = lpt[1]; /* number of coeffs per component */
		na   = lpt[2]; /* number of sets of coeffs per 32day interval */

		interp(buff, t, ncf, 3, na, LList[11], dumpv);

		pv[0][10] = dumpv[0][0];
		pv[1][10] = dumpv[1][0];
		pv[2][10] = dumpv[2][0];
		pv[3][10] = dumpv[0][1];
		pv[4][10] = dumpv[1][1];
		pv[5][10] = dumpv[2][1];
	}
}

////////计算3*3矩阵相对于给定第Ith坐标轴角度的旋转/////////////////
//这个函数用来计算一个3×3的旋转矩阵及它的导数绕第ith坐标轴旋转一个角度phi后所得到的矩阵。
//其中i=1时x轴；i=2时y轴；i＝3时，z轴。
void visual_position_calculate::RotMat(int axis, double phi, DMatrix r, DMatrix dr)
{
	double cosphi, sinphi;

	cosphi = cos(phi);
	sinphi = sin(phi);

	switch (axis) {
	case 1:  /* rotate about x-axis */
		r[0][0] = 1.0;
		r[0][1] = 0.0;
		r[0][2] = 0.0;
		r[1][0] = 0.0;
		r[1][1] =  cosphi;
		r[1][2] =  sinphi;
		r[2][0] = 0.0;
		r[2][1] = -sinphi;
		r[2][2] =  cosphi;
		dr[0][0] = 0.0;
		dr[0][1] = 0.0;
		dr[0][2] = 0.0;
		dr[1][0] = 0.0;
		dr[1][1] = -sinphi;
		dr[1][2] =  cosphi;
		dr[2][0] = 0.0;
		dr[2][1] = -cosphi;
		dr[2][2] = -sinphi;
		break;
	case 2:  /* rotate about y-axis */
		r[0][0] =  cosphi;
		r[0][1] = 0.0;
		r[0][2] = -sinphi;
		r[1][0] = 0.0;
		r[1][1] = 1.0;
		r[1][2] = 0.0;
		r[2][0] =  sinphi;
		r[2][1] = 0.0;
		r[2][2] =  cosphi;
		dr[0][0] = -sinphi;
		dr[0][1] = 0.0;
		dr[0][2] = -cosphi;
		dr[1][0] = 0.0;
		dr[1][1] = 0.0;
		dr[1][2] = 0.0;
		dr[2][0] =  cosphi;
		dr[2][1] = 0.0;
		dr[2][2] = -sinphi;
		break;
	case 3:  /* rotate about z-axis */
		r[0][0] =  cosphi;
		r[0][1] =  sinphi;
		r[0][2] = 0.0;
		r[1][0] = -sinphi;
		r[1][1] =  cosphi;
		r[1][2] = 0.0;
		r[2][0] = 0.0;
		r[2][1] = 0.0;
		r[2][2] = 1.0;
		dr[0][0] = -sinphi;
		dr[0][1] =  cosphi;
		dr[0][2] = 0.0;
		dr[1][0] = -cosphi;
		dr[1][1] = -sinphi;
		dr[1][2] = 0.0;
		dr[2][0] = 0.0;
		dr[2][1] = 0.0;
		dr[2][2] = 0.0;
		break;
	default:
		exit(1);
	}
}
void visual_position_calculate::convert_little_endian(char *ptr, int len) 
{
	if (BIG_ENDIAN_TEST) 
	{
		int i;
		char *tmp;

		if ((tmp = (char*)malloc(len)) == NULL) {
			exit(1);
		}

		for (i = 0; i < len; i++) {
			tmp[i] = ptr[(len-1)-i];
		}

		memcpy(ptr, tmp, len);
		free(tmp);
	}
}
/////////将给定的数值分解为小数和整数部分////////////////////
void visual_position_calculate::split(double tt, double *ipart, double *fpart)
{
	*ipart = floor(tt);
	*fpart = tt - *ipart;

	if ((tt < 0) && (*fpart != 0)) {
		*ipart = *ipart - 1.0;
		*fpart = *fpart + 1.0;
	}
}
///////////////////区分和内插切比雪夫多项式来给定位置和速度矢量/////////////////////////////////
void visual_position_calculate::interp(int buff, double *t, int ncf, int ncm, int na, int fl, double dumpv[][2])
{
	int i, j, l, n, m;
	

	if(flag_interp)
	{
		bcoef=0;
		np=0;
		nv=0;
		memset(pc,0,sizeof(pc));
		memset(vc,0,sizeof(vc));
		memset(cbody,0,sizeof(cbody));
		twot=0.0;
		dna=0.0;
		temp=0.0;
		ll=0.0;
		tc=0.0;
		vfac=0.0;
		flag_interp=false;
	}

	np = 2;
	nv = 3;
	twot = 0.0;
	pc[0] = 1.0;
	pc[1] = 0.0;
	vc[1] = 1.0;

	dna = (double) na;
	dt1 = floor(t[0]);
	temp = dna*t[0];
	ll = floor((temp - dt1) + 1.0);

	/* 'tc' is the normalized chebyshev time (-1 <= tc <= 1) */
	tc = 2.0*(fmod(temp, 1.0) + dt1) - 1.0;

	if (tc != pc[1]) {
		np = 2;
		nv = 3;
		pc[1] = tc;
		twot = tc + tc;
	}

	if (np < ncf) {
		for (i=np; i<ncf; i++) {
			pc[i] = twot * pc[i-1] - pc[i-2];
		}
		np = ncf;
	}

	/* interpolate to get position for each component */

	/* number of coefficients for body */
	bcoef= ncf*na*ncm;

	/* stored body's coefficients in an array */
	n = buff;
	for (m=0; m<bcoef; m++) {
		cbody[m] = db[n];
		n++;
	}

	/* fill the cbuf[][][] array */
	n = 0;
	/* loop for each sub-interval */
	for (l=0; l<na; l++) {
		/* loop for each component */
		for (i=0; i<ncm; i++) {
			/* loop for each set of coeffs */
			for (j=0; j<ncf; j++) {
				cbuf[j][i][l] = cbody[n];
				n++;
			}
		}
	}

	for (i=0; i<ncm; i++) {
		dumpv[i][0] = 0.0;
		for (j=ncf-1; j>=0; j--) {
			dumpv[i][0] = dumpv[i][0] + pc[j] * cbuf[j][i][((int)ll)-1];
		}
	}
	if (fl <= 1) return;

	/*
	If velocity interpolation is wanted, be sure enough
	derivative polynomials have been generated and stored.
	*/

	vfac = (dna+dna)/t[1];
	vc[2] = twot+twot;
	if (nv < ncf) {
		for (i=nv; i<ncf; i++) {
			vc[i]=twot * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2];
		}
		nv = ncf;
	}

	/* interpolate to get velocity for each component */

	for (i=0; i<ncm; i++) {
		dumpv[i][1] = 0.0;
		for (j=ncf-1; j>=1; j--) {
			dumpv[i][1] = dumpv[i][1] + vc[j] * cbuf[j][i][((int)ll)-1];
		}
		dumpv[i][1] = dumpv[i][1] * vfac;
	}
}

///////////报错/////////////////////////
void visual_position_calculate::errprt(int i, char const *msg)
{
	printf("%s\n",msg);
	//	exit(1);
}

////////////释放矩阵m的所有内存////////////////
void visual_position_calculate::freeDMatrix(DMatrix m, int n)
{
	int i;

	for (i=0; i<n; i++)
	{
		free(m[i]); /* storage for doubles in row i */
	}
	free(m); /* storage for pointers to rows */

	return;
}



void visual_position_calculate::FinalCalculate(int starNum, double starRawData[][8], double starViualPos[][4])
{
	if(timeReady)
	{
		double jed1 = 0.0, julianCenturies = 0.0, S0 = 0.0, S1 = 0.0;
		double jed = 0.0;
		double Ltime = 0.0;

		//double inside;

		double P1[3], P2[3], P3[3], eb[3], ebdot[3], earth_ssb[6], earth_hel[6], body_geo[6], body_hel[6];

		//double E[3],RQB[6],Q[3],RSB[6],Qdot[3],Pdot[3],unit_body_geo[6],unit_body_hel[6],Edot[3],magP,magQ,RP[6],RQ[6];
		//double RA,DE,sinRA,cosRA,sinDE,cosDE;
		//double plx,r,muRA,muDE,rdot,drda[3],drdd[3],DT,trueGeoDist;

		memset(P1, 0, sizeof(P1));
		memset(P2, 0, sizeof(P2));
		memset(P3, 0, sizeof(P3));
		memset(eb, 0, sizeof(eb));
		memset(ebdot, 0, sizeof(ebdot));
		memset(earth_ssb, 0, sizeof(earth_ssb));
		memset(earth_hel, 0, sizeof(earth_hel));
		memset(body_geo, 0, sizeof(body_geo));
		memset(body_hel, 0, sizeof(body_hel));


		jed1 = Cal2JED(nMonth, nDay, nYear, 0., m_tai_utc);																																		 //世界时为0时的儒略时
		julianCenturies = (jed1 - J2000) / JulCty;																																	 //从J2000.0起算的儒略世纪数
		S0 = 6.0 + 41.0 / 60.0 + 50.54841 / 3600.0 + 8640184.812866 / 3600.0 * julianCenturies + 0.093104 / 3600.0 * julianCenturies * julianCenturies - 0.0000062 / 3600.0 * julianCenturies * julianCenturies * julianCenturies; //世界时0时的格林尼治平恒星时
		//S1=S0+hour+0.002737909350795*hour;//世界时等于hour时的格林尼治恒星时

		S1 = S0 + (double)nHour + (double)nMinute / 60.0 + (double)nSecond / 3600.0 + 0.002737909350795 * ((double)nHour + (double)nMinute / 60.0 + (double)nSecond / 3600.0);

		greenwichHourAngle = (S1 * 15.0 / 360.0 - int(S1 * 15.0 / 360.0)) * 360.0; //春分点格林时角为姿态计算准备

		jed = Cal2JED(nMonth, nDay, nYear, m_utc, m_tai_utc); //将世界时转化为儒略日期
		m_JD = jed;

		GetStateVector(jed, 3, 12, 1, StateVector); //得到星体的状态矢量
		for (int i = 0; i < 6; i++)
		{
			earth_ssb[i] = StateVector[3 - 1][12 - 1][1 - 1][i];
		}

		GetStateVector(jed, 3, 11, 1, StateVector); //得到星体的状态矢量
		//上式3为星体编号,对恒星视位置而言,3为地球
		for (int i = 0; i < 6; i++)
		{
			earth_hel[i] = StateVector[3 - 1][11 - 1][1 - 1][i];
		}

		//岁差章动修正,多星体共用//

		double zeta = 0.0, z = 0.0, theta = 0.0, zetadot = 0.0, zdot = 0.0, thetadot = 0.0;
		double dpsi = 0.0, deps = 0.0, dpsidot = 0.0, depsdot = 0.0;
		double trueEps = 0.0, MeanEps = 0.0, trueEpsDot = 0.0, MeanEpsDot = 0.0;

		GetPrecessParams(J2000, jed, &zeta, &z, &theta, &zetadot, &zdot, &thetadot); //计算岁差参数
		GetDpsiDeps(jed, &dpsi, &deps, &dpsidot, &depsdot); //计算章动dpsi, deps, dpsidot, and depsdot.
		Obliquity(J2000, jed, 0, &MeanEps, &MeanEpsDot);	//计算平黄赤交角及它的导数
		Obliquity(J2000, jed, 1, &trueEps, &trueEpsDot);
		SplitStateVector(earth_ssb, eb, ebdot);

		double StarData[6];		//用于计算的星点相关信息
		double X = 0.0, Y = 0.0, Z = 0.0;	//直角坐标
		double alpha = 0.0; //恒星赤经
		double delta = 0.0;  //恒星赤纬

		for (int i = 0; i < starNum; i++) //循环计算识别星的视位置
		{
			alpha = 0.0;
			delta = 0.0;
			jed = m_JD;
						
			StarData[0] = starRawData[i][2]; //赤经(弧度)
			StarData[1] = starRawData[i][3]; //赤纬(弧度)
			StarData[2] = starRawData[i][7]; ///视差
			StarData[3] = starRawData[i][4]; //赤经自行
			StarData[4] = starRawData[i][5]; //赤纬自行
			StarData[5] = starRawData[i][6]; //视向速度

			LightTime(jed, 99, earth_ssb, earth_hel, StarData, body_geo,
					  body_hel, &Ltime); //将星体由几何太阳系质心位置转化到几何位置或日心位置
			//LightTime函数将天体的质心坐标转化为地心坐标函数,修正周年视差和自行

			RayBend(earth_hel, body_geo, body_hel, P1); /////根据太阳的重力场来计算相对光偏差
			//RayBend函数为相对论的光偏差修正函数,用来计算由太阳系引力所引起的光偏差

			Aberrate(P1, ebdot, P2); //光行差的修正函数用来计算光行差


			//通过给定的角度来计算矩阵的旋转后的坐标矢量（右手法则）
			//QRotate这个函数的功能基本上与Rotmat函数的功能一致，只不过它是6×6矩阵，并且遵循右手法则。
			
			//岁差章动修正
			QRotate(P2, 3, -zeta, -zetadot, 1, P3); 
			QRotate(P3, 2, theta, thetadot, 1, P3);
			QRotate(P3, 3, -z, -zdot, 1, P3);
			QRotate(P3, 1, MeanEps, MeanEpsDot, 1, P3);
			QRotate(P3, 3, -dpsi, -dpsidot, 1, P3);
			QRotate(P3, 1, -trueEps, -trueEpsDot, 1, P3);

			//直角坐标
			X = P3[0];
			Y = P3[1];
			Z = P3[2];

			//转换为球面坐标——赤经赤纬(弧度)
			alpha = (atan2(Y, X));
			if (alpha < 0)
				alpha = 2.0 * PI + alpha;
			delta = atan(Z / (sqrt(X * X + Y * Y)));


			starViualPos[i][0] = starRawData[i][0]; 	//通道号
			starViualPos[i][1] = starRawData[i][1];		//星号
			starViualPos[i][2] = alpha;   //赤经(弧度)
			starViualPos[i][3] = delta;	//赤纬(弧度)
		}
	}
	else
	{
		/* code */
	}
	
}