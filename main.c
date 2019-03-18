#include "io.h"
#include "cvtools.h"
#include<string.h>
#include "ds.h"
#include "sift.h"
#include <stdbool.h>
#include "fast-edge.h"
#include <math.h>
#include "main.h"
#include "time.h"
#define MAX_NUM 500

ClImage* Sobel(ClImage *img, int *Window, float threshold)
{
	int Height = Window[1] - Window[0];
	int Width = Window[3] - Window[2];
	int width = img->width;
	ClImage* dst = (ClImage*)malloc(sizeof(ClImage));
	dst->imageData = (uchar*)malloc(sizeof(uchar)*Height*Width);
	dst->height = Height;
	dst->width = Width;
	dst->channels = 1;
	int Sobel_X[9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
	int Sobel_Y[9] = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };
	int temp1,temp2,temp3;
	for (int i = Window[0]; i < Window[1]; i++)
	{
		for (int j = Window[2]; j < Window[3]; j++)
		{
			dst->imageData[(i - Window[0])*Width + j - Window[2]] = 0;
			temp1 = Sobel_X[0] *  img->imageData[(i - 1)*width+j - 1] + Sobel_X[1] * img->imageData[(i - 1)*width + j ] + Sobel_X[2] * img->imageData[(i - 1)*width + j + 1] + Sobel_X[3] * img->imageData[i*width + j - 1] + Sobel_X[4] * img->imageData[i*width + j] + Sobel_X[5] * img->imageData[i*width + j + 1] + Sobel_X[6] * img->imageData[(i + 1)*width + j - 1] + Sobel_X[7] * img->imageData[(i + 1)*width + j ] + Sobel_X[8] * img->imageData[(i + 1)*width + j + 1];
			temp2 = Sobel_Y[0] * (int)img->imageData[(i - 1)*width + j - 1] + Sobel_Y[1] * (int)img->imageData[(i - 1)*width + j] + Sobel_Y[2] * (int)img->imageData[(i - 1)*width + j + 1] + Sobel_Y[3] * (int)img->imageData[i*width + j - 1] + Sobel_Y[4] * (int)img->imageData[i*width + j] + Sobel_Y[5] * (int)img->imageData[i*width + j + 1] + Sobel_Y[6] * (int)img->imageData[(i + 1)*width + j - 1] + Sobel_Y[7] * (int)img->imageData[(i + 1)*width + j] + Sobel_Y[8] * (int)img->imageData[(i + 1)*width + j + 1];
			temp3 = sqrt(temp2*temp2+temp2*temp2);
			if (temp3 > threshold)
			{
				dst->imageData[(i - Window[0])*Width + j - Window[2]] = 255;
			}
			else
			{
				dst->imageData[(i - Window[0])*Width + j - Window[2]] = 0;
			}
		}
	}
	return dst;
}
ClImage* Image_cut(ClImage *img, int *Window)
{
	int Height = Window[1] - Window[0];
	int Width = Window[3] - Window[2];
	int width = img->width;
	ClImage* dst = (ClImage*)malloc(sizeof(ClImage));
	dst->imageData = (uchar*)malloc(sizeof(uchar)*Height*Width);
	dst->height = Height;
	dst->width = Width;
	dst->channels = 1;
	for (int i = Window[0]; i < Window[1]; i++)
	{
		for (int j = Window[2]; j < Window[3]; j++)
		{
			dst->imageData[(i - Window[0])*Width + j - Window[2]] = img->imageData[i*width + j];
		}
	}
	return dst;
}
int* Fix_line(ClImage * img)
{
	int Height = img->height;
	int Width = img->width;
	int * Two_line = (int *)malloc((Width-6) * 2 * sizeof(int));
	int index_min, index_max;
	int start_hl, start_hh;
	bool flag1 = false, flag2 = false;
	//find the start point
	for (int i = 3; i < Width-3; i++)
	{
		for (int j = 3; j < Height-3; j++)
		{
			if (img->imageData[j*Width + i] == 255 && j<Height/2 && flag1==0)
			{
				Two_line[0] = j;
				flag1 = 1;
			}
			if (img->imageData[j*Width + i] == 255 && j>Height / 2 && flag2 == 0)
			{
				Two_line[1] = j;
				flag2 = 1;
			}
		}
		if (flag1 && flag2)
		{
			break;
		}
	}
	//fix the line
	for (int i = 4; i < Width-3; i++)
	{
		flag2 = false;
		index_min = Height / 2;
		index_max = Height / 2;
		for (int j = 3; j < Height-3; j++)
		{
			if (img->imageData[j*Width + i] == 255 && j<Height / 2)
			{
				index_min = j;
			}
			if (img->imageData[j*Width + i] == 255 && j>Height / 2 && flag2==0)
			{
				index_max = j;
				flag2 = 1;
			}
		}
		if(index_min == Height/2)
		{
			Two_line[2 * (i-3)] = Two_line[2 * (i - 4)];
		}
		else
		{
			Two_line[2 * (i - 3)] = index_min;
		}
		if (index_max == Height / 2)
		{
			Two_line[2 * (i - 3) + 1] = Two_line[2 * (i - 4) + 1];
		}
		else
		{
			Two_line[2 * (i - 3) + 1] = index_max;
		}
	}
	free(img);
	return Two_line;
}
ClImage * Mean_curvature(ClImage *img)
{
	int Height = img->height;
	int Width = img->width;
	ClImage* img_out = (ClImage*)malloc(sizeof(ClImage));
	img_out->imageData = (uchar*)malloc(sizeof(uchar)*Height*Width);
	img_out->height = Height;
	img_out->width = Width;
	img_out->channels = 1;
	int fx, fy, fxx, fyy,fxy;
	//for(int i=1;i<)
	return img_out;
}
int main(int argc, char* argv[])
{
	int a = 3;
	int b = 3;
	int c = sqrt(0.3);
	/*******前期裁剪*******/
	/*******预处理，高斯模糊函数******/
	ClImage* img = clLoadImage("C:\\Users\\daiguomeng\\Desktop\\Finger\\Finger Vein Database\\Finger Vein Database\\1\\index_1.bmp");
	Size size = { 0,0 };
	double sig_diff = 0.8;
	//int Window[4] = { 0,240,0,320 };
	int Window[4] = { 40,225,15,301 };
	ClImage* img_cut = Image_cut(img, Window);
	//free(img->imageData);
	//free(img);
	ClImage* Blur_img = GaussianBlur(img_cut, size, sig_diff, sig_diff);
	free(img_cut->imageData);
	free(img_cut);;
	clSaveImage("index_1_blur.bmp", Blur_img);
	//ClImage *Blur_img_canny = canny_edge_detect(Blur_img);
	//int *two_line = Fix_line(Blur_img_canny);

	//ClImage *Img_vein = Mean_curvature(Blur_img);
	/********预处理，裁剪*************/
	//clSaveImage("detectpt1.bmp", Blur_img_canny);
	free(Blur_img->imageData);
	free(Blur_img);
	//free(Blur_img_canny->imageData);
	//free(Blur_img_canny);
	uchar color[7][3] = { { 255, 0, 0 }, { 255, 0, 255 }, { 0, 127, 255 },
	{ 11, 134, 184 }, { 0, 255, 140 }, { 237, 149, 100 }, { 255, 255, 0 } };
	uchar pt_color[3] = { 0, 255, 255 };
	clock_t t1 = clock();
	//ClImage* img1 = clLoadImage("C:\\Users\\daiguomeng\\Desktop\\Finger\\Finger Vein Database\\Finger Vein Database\\2\\index_1.bmp");
	KeyPoint* keypoints1 = (KeyPoint*)malloc(sizeof(KeyPoint)*MAX_NUM);
	Descriptor* descp1 = detect(img, keypoints1, false);
	clock_t t2 = clock();
	printf("%d", t2 - t1);
	ClImage* img2 = clLoadImage("C:\\Users\\daiguomeng\\Desktop\\Finger\\Finger Vein Database\\Finger Vein Database\\10\\index_2.bmp");
	KeyPoint* keypoints2 = (KeyPoint*)malloc(sizeof(KeyPoint)*MAX_NUM);
	Descriptor* descp2 = detect(img2, keypoints2, false);
	
	write(keypoints1, descp1, "pt1.txt", "descp1.txt");
	write(keypoints2, descp2, "pt2.txt", "descp2.txt");

	ClImage* cmb_img = combine(img, img2);
	int matchp[5000][2];	
	int len = match(img, img2, descp1, descp2, keypoints1, keypoints2, matchp);
	printf("Get %d matched points.\n", len);

	int i;
	for (i = 0; i < len; i++)
	{
		Point2d p1 = { (int)(keypoints1[matchp[i][0]].pt.x), (int)(keypoints1[matchp[i][0]].pt.y) };
		draw(cmb_img, p1, pt_color);
		Point2d p2 = { (int)(keypoints2[matchp[i][1]].pt.x+img->width), (int)(keypoints2[matchp[i][1]].pt.y) };
		draw(cmb_img, p2, pt_color);
		drawLine(cmb_img, p1, p2, color[i % 7]);
	}
	clSaveImage("match.bmp", cmb_img);

	for (i = 0; i < descp1->height; i++)
	{
		Point2d p = { (int)(keypoints1[i].pt.x), (int)(keypoints1[i].pt.y) };
		draw(img, p, pt_color);
	}
	clSaveImage("detectpt1.bmp", img);

	for (i = 0; i < descp2->height; i++)
	{
		Point2d p = { (int)(keypoints2[i].pt.x), (int)(keypoints2[i].pt.y) };
		draw(img2, p, pt_color);
	}
	clSaveImage("detectpt2.bmp", img2);
	system("pause");
}