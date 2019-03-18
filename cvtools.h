#include "io.h"
#include "ds.h"


ClImage* BGR2Gray(ClImage* img);

ClImage* resize_bilinear(ClImage* src, Size dsize, double fx, double fy);

ClImage* resize_nearest(ClImage* src, Size dsize, double fx, double fy);


void getClImageInfo(ClImage* img, int* height, int* width, int* channel);
void getFtImageInfo(FtImage* img, int* height, int* width, int* channel);

void printClImageInfo(ClImage* img);
void printFtImageInfo(FtImage* img);

Mat* getGaussianKernel(int n, double sigma);
ClImage* GaussianBlur(ClImage* src, Size ksize, double sigma1, double sigma2);

ClImage* combine(ClImage* img1, ClImage* img2);
void drawLine(ClImage* img, Point2d p1, Point2d p2, uchar* color);

void write(KeyPoint* keypoints, Descriptor* descriptor, char* pt_path, char* descp_path);
void draw(ClImage* img, Point2d point, uchar* color);
float L2dis(float* p1, float* p2, int n);