#include <cv.h>
#include <highgui.h>

#define XBLOCK 30	//x方向分块个数
#define YBLOCK 20	//y方向分块个数
#define MAX_RADIUS 15	//星象最大半径
#define MIN_STAR_RECOGNIZED 2 //识别成功可进行计算的最少恒星数量
#define MAX_STAR_RECOGNIZED 32 //识别成功可进行计算的最大恒星数量
#define ALPHA 4.5	//阈值系数（方差倍数）

class image_process
{
    image_process();

public:

    int starsInImage[MAX_STAR_RECOGNIZED][6];
    inline int getStarNum()		{	return N_star;	};
    void ImageProcess(IplImage *rawImage, IplImage *starImage, double &thresh, double alpha); //阈值分割
    void Centroid(IplImage *rawImage, int cameraNum, bool withThresh); //求取质心

private:

	int imageWidth;
	int imageHeight;
    int N_star;

    void SplitImage(const IplImage *imgSrc, double &thresh, double alfa); //分块法求取阈值
    void PointProcess(IplImage *plate); //滤除单点
};


void image_process::ImageProcess(IplImage* rawImage, IplImage* starImage, double &thresh, double alpha)
{
	cvZero(starImage);	//图像清零

	SplitImage(rawImage,thresh, alpha);	//确定阈值

	cvThreshold(rawImage,starImage,thresh,255.0,CV_THRESH_TOZERO);	//阈值分割

	PointProcess(starImage);	//滤除单点

	return;
}


image_process::image_process()
{
    imageWidth = 0;
	imageHeight = 0;
    N_star = 0;
}


void image_process::SplitImage(const IplImage *imgSrc, double &thresh, double alpha=ALPHA)
{
	double mean = 0.0, delta = 0.0;
	IplImage *img = cvCloneImage(imgSrc);
	int subWidth  = cvFloor(img->width/XBLOCK);		//块宽度
	int subHeight = cvFloor(img->height/YBLOCK);	//快高度

	for (int y=0; y<YBLOCK; y++)
	{
		for (int x=0; x<XBLOCK; x++)
		{
			cvSetImageROI(img, cvRect(x*subWidth, y*subHeight, subWidth, subHeight));
			IplImage *roiImg = cvCreateImage(cvSize(subWidth,subHeight),img->depth,1);
			cvCopy(img, roiImg);

			CvScalar mean2;
			CvScalar delta2;

			cvAvgSdv(roiImg, &mean2, &delta2);	//计算小块的均值标准差

			double l = 0.0;
			l = mean2.val[0] + delta2.val[0]*alpha;		//计算阈值

			if (l < thresh || thresh == 0.0)
			{
				thresh = l;				//选择最小阈值为全局阈值
			}
			cvResetImageROI(img);	
		}
	}

	cvReleaseImage(&img);

	return;
}

void image_process::PointProcess(IplImage* plate)
{
	char* plateData = (char*)plate->imageData;
	int widthStep = plate->widthStep;

	int tempSum = 0;

	for (int j = 1; j < plate->height - 1; j++)
	{
		for (int i = 1; i < plate->width - 1; i++)
		{		
			if (plateData[j * widthStep + i] == 0)
				continue;
			else
			{
				tempSum = plateData[(j-1) * widthStep + (i-1)] +
					plateData[(j-1) * widthStep + i] +
					plateData[(j-1) * widthStep + (i+1)] +
					plateData[j * widthStep + (i-1)] +
					plateData[j * widthStep + (i+1)] +
					plateData[(j+1) * widthStep + (i-1)] +
					plateData[(j+1) * widthStep + i] +
					plateData[(j+1) * widthStep + (i+1)];		//八邻域所有点灰度值之和

				if (tempSum  == 0)			//和为零该点为孤立单个像素,滤除
				{
					plateData[j * widthStep + i] = 0;		
				}
			}
		}
	}
	return;
}

void image_process::Centroid(IplImage* rawImage, int cameraNum, bool withThresh)
{
	double alpha = ALPHA;
	double thresh = 0.0;

	IplImage * threshImage = cvCloneImage(rawImage);

	imageWidth = rawImage->width;
	imageHeight = rawImage->height;

	cvZero(threshImage);
 	ImageProcess(rawImage, threshImage, thresh, alpha);

	/*cvNamedWindow("1");
	cvShowImage("1",threshImage);
	cvWaitKey(0);
	cvDestroyWindow("1");*/
	//cvSaveImage("thresh.bmp",threshImage);

	int width = threshImage->width;
	int height = threshImage->height;
	int widthStep = threshImage->widthStep;

	IplImage* dst = cvCloneImage(threshImage);
	IplImage* starNumPlate = cvCreateImage(cvSize(width,height), IPL_DEPTH_64F, 1);

	cvZero(starNumPlate);

	char* plateData = (char*)threshImage->imageData;
	double* starNumPlateData = (double*)starNumPlate->imageData;
	char* dstData = (char*)dst->imageData;

	int count = 0;

	//获得图像所有连续区域,每个区域一个编号
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			if (plateData[j*widthStep + i] == 0)
				continue;
			else
			{
				if(j!=0)
				{
					if(i!=0)
					{
						if ((starNumPlateData[(j-1) * widthStep + (i-1)] > 0))
						{
							starNumPlateData[j * widthStep + i] = starNumPlateData[(j-1) * widthStep + (i-1)];
							continue;
						}
						if ((starNumPlateData[j * widthStep + (i-1)] > 0))
						{
							starNumPlateData[j * widthStep + i] = starNumPlateData[j * widthStep + (i-1)];
							continue;
						}
					}
					if ((starNumPlateData[(j-1) * widthStep + i] > 0))
					{
						starNumPlateData[j * widthStep + i] = starNumPlateData[(j-1) * widthStep + i];
						continue;
					}
					if(i!=width-1)
					{
						if ((starNumPlateData[(j-1) * widthStep + (i+1)] > 0))
						{
							starNumPlateData[j * widthStep + i] = starNumPlateData[(j-1) * widthStep + (i+1)];
							continue;
						}
						if ((starNumPlateData[j * widthStep + (i+1)] > 0))
						{
							starNumPlateData[j * widthStep + i] = starNumPlateData[j * widthStep + (i+1)];
							continue;
						}
					}
				}

				if(i!=0)
				{
					if ((starNumPlateData[j * widthStep + (i-1)] > 0))
					{
						starNumPlateData[j * widthStep + i] = starNumPlateData[j * widthStep + (i-1)];
						continue;
					}
				}
				if(i!=width-1)
				{
					if ((starNumPlateData[j * widthStep + (i+1)] > 0))
					{
						starNumPlateData[j * widthStep + i] = starNumPlateData[j * widthStep + (i+1)];
						continue;
					}
				}

				if(j!=height-1)
				{
					if(i!=0)
					{
						if ((starNumPlateData[(j+1) * widthStep + (i-1)] > 0))
						{
							starNumPlateData[j * widthStep + i] = starNumPlateData[(j+1) * widthStep + (i-1)];
							continue;
						}
					}
					if ((starNumPlateData[(j+1) * widthStep + i] > 0))
					{
						starNumPlateData[j * widthStep + i] = starNumPlateData[(j+1) * widthStep + i];
						continue;
					}
					if(i!=width-1)
					{
						if ((starNumPlateData[(j+1) * widthStep + (i+1)] > 0))
						{
							starNumPlateData[j * widthStep + i] = starNumPlateData[(j+1) * widthStep + (i+1)];
							continue;
						}
					}
				}
				
				if(starNumPlateData[j * widthStep + i] == 0)
				{
					count = count+1;
					starNumPlateData[j * widthStep + i] = count;
				}
			}
		}
	}

	double bright;

	uchar l = 0;
	int x=0, y=0;

	for (int k = 0; k < MAX_STAR_RECOGNIZED; k++)        //最大星数量
	{
		//StarInfo si;
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				if ((uchar)(plateData[j * widthStep + i]) > l)
				{
					l = (uchar)(plateData[j * widthStep + i]);
					x = i;
					y = j;
				}
			}
		}

 		if (l == 0.0)
			goto line1;


		bright = starNumPlateData[y * widthStep + x];

		starsInImage[N_star][0] = cameraNum;	//通道号,取值范围{1,2,3}
		starsInImage[N_star][1] = N_star;			//通道内序号
		starsInImage[N_star][4] = l;			//星等,即灰度值


		l = 0.0;
		double r = 0.0;	
		double x_sum = 0.0, y_sum = 0.0, gray_sum = 0.0, tmp = 0.0;

		for (int j = 0; j < height; j++)
		{
			if(j<0 || j>height-1)
				continue;
			for (int i = 0; i < width; i++)
			{
				if(i<0 || i>width-1)
					continue;
				if (starNumPlateData[j * widthStep + i] == bright)
				{
					plateData[j * widthStep + i] = 0;

					if (r < pow((i-x),2.0) + pow((j-y),2.0))
						r = pow((i-x),2.0) + pow((j-y),2.0);
				}
			}
		}

		int _r = cvFloor(sqrt(r));

		starsInImage[N_star][5]=_r;

		for (int j = y-_r; j < y+_r+1; j++)
		{
			if(j<0 || j>height-1)
				continue;
			for (int i = x-_r; i < x+_r+1; i++)
			{
				if(i<0 || i>width-1)
					continue;
				plateData[j * widthStep + i] = 0;
			}
		}

		
		

		for (int j = y-_r; j < y+_r+1; j++)
		{
			if(j<0 || j>height-1)
				continue;
			for (int i = x-_r; i < x+_r+1; i++)
			{
				if(i<0 || i>width-1)
					continue;
				tmp = (uchar)dstData[j * widthStep + i];
				if(withThresh)
				{
					tmp = tmp - thresh;
					if(tmp<0.0)
						tmp = 0.0;
				}

				x_sum += i*tmp;
				y_sum += j*tmp;

				gray_sum += tmp;

			}
		}

		starsInImage[N_star][2] = x_sum/gray_sum;		//质心横坐标
		starsInImage[N_star][3] = y_sum/gray_sum;		//质心纵坐标

		N_star++;			//星总数+1

	}

line1:

	/*FILE *pFile;
	pFile = fopen("g2.txt","w");
	for(int u = 0; u < N_star; u++)
	{
		fprintf(pFile,"double star%d = {%d, %d, %f, %f, %f, %f};\n", u+1, (int)starsInImage[u][0], (int)starsInImage[u][1], starsInImage[u][2], starsInImage[u][3], starsInImage[u][4], starsInImage[u][5]);
	}
	fclose(pFile);*/

	cvReleaseImage(&threshImage);
	cvReleaseImage(&dst);
	cvReleaseImage(&starNumPlate);
	return;
}