
#include "cvtools.h"
#include "math.h"
#include "stdio.h"

void getClImageInfo(ClImage* img, int* height, int* width, int* channel)
{
	*height = img->height;
	*width = img->width;
	*channel = img->channels;
}

void getFtImageInfo(FtImage* img, int* height, int* width, int* channel)
{
	*height = img->height;
	*width = img->width;
	*channel = img->channels;
}

void printClImageInfo(ClImage* img)
{
	int height, width, channel;
	getClImageInfo(img, &height, &width, &channel);
	printf("Image height: %d\n", height);
	printf("Image width: %d\n", width);
	printf("Image channel: %d\n", channel);
	printf("Image data(from top-left to bottom-right):\n");
	int i;
	if (channel == 3)
	{
		// print the image data in [H, W, C] format, BGR order
		for (i = 0; i < height*width*channel; i++)
		{
			if (i>0 && i%channel == 0)
				printf("\n");
			printf("%d ", (int)(img->imageData[i]));
		}
	}
	else if (channel == 1)
	{
		for (i = 0; i < height*width; i++)
		{
			if (i>0 && i%width == 0)
				printf("\n");
			printf("%d ", (int)(img->imageData[i]));
		}
	}
	printf("\n");
	
}

void printFtImageInfo(FtImage* img)
{
	int height, width, channel;
	getFtImageInfo(img, &height, &width, &channel);
	printf("Image height: %d\n", height);
	printf("Image width: %d\n", width);
	printf("Image channel: %d\n", channel);
	printf("Image data(from top-left to bottom-right):\n");
	int i;
	if (channel == 3)
	{
		// print the image data in [H, W, C] format, BGR order
		for (i = 0; i < height*width*channel; i++)
		{
			if (i>0 && i%channel == 0)
				printf("\n");
			printf("%.7f ", (img->imageData[i]));
		}
	}
	else if (channel == 1)
	{
		for (i = 0; i < height*width; i++)
		{
			if (i>0 && i%width == 0)
				printf("\n");
			printf("%.7f ", (img->imageData[i]));
		}
	}
	printf("\n");

}

ClImage* BGR2Gray(ClImage* src)
{
	int i;
	int width, height, channel;
	getClImageInfo(src, &height, &width, &channel);
	uchar* in_img = src->imageData;
	ClImage* bmpImg = (ClImage*)malloc(sizeof(ClImage));
	bmpImg->imageData = (uchar*)malloc(sizeof(uchar)*width*height);
	bmpImg->width = width;
	bmpImg->height = height;
	bmpImg->channels = 1;
	for (i = 0; i < width*height; i++)
	{
		float gray = 0.114 * in_img[i*channel] + 0.587 * in_img[i*channel + 1] + 0.299 * in_img[i*channel + 2];
		if (gray - (int)(gray) < 0.5)
			bmpImg->imageData[i] = floor(gray);
		else
			bmpImg->imageData[i] = ceil(gray);
	}

	return bmpImg;
}

// resize using bilinear interpolation
ClImage* resize_bilinear(ClImage* src, Size dsize, double fx, double fy)
{
	int width, height, channel;
	int i, j, k;
	getClImageInfo(src, &height, &width, &channel);
	uchar* srcData = src->imageData; // read only
	if (dsize.height != 0 && dsize.width != 0)
	{
		fx = (double)dsize.width / (double)width;
		fy = (double)dsize.height / (double)height;
	}
	else if (fx != 0 && fy != 0)
	{
		dsize.width = (int)(width * fx);
		dsize.height = (int)(height * fy);
	}
	else
	{
		//
	}

	ClImage* dst = (ClImage*)malloc(sizeof(ClImage));
	dst->width = dsize.width;
	dst->height = dsize.height;
	dst->channels = channel;
	dst->imageData = (uchar*)malloc(sizeof(uchar)*dsize.width*dsize.height*channel);
	
	int dst_step = dsize.width * channel;
	int src_step = width * channel;
	for (i = 0; i < dsize.height; i++)
	{
		//uchar* dstData = resize_dst + i * dst_step;
		float srcy = (float)((i + 0.5) / fy - 0.5);
		int y = floor(srcy);
		float v = srcy - y;
		y = min(y, height - 2);
		y = max(0, y);
		/*if (y < 0)
		{
			y = 0;
			v = 0;
		}
		if (y >= height - 1)
		{
			y = height - 2;
			v = 0;
		}*/
		short s_minus_v = (short)((1. - v) * 2048); // 1-v
		short s_v = 2048 - s_minus_v;

		for (j = 0; j < dsize.width; j++)
		{
			float srcx = ((j + 0.5) / fx - 0.5);
			int x = floor(srcx);
			float u = srcx - x;
			if (x < 0)
			{
				x = 0;
				u = 0;
			}
			if (x >= width - 1)
			{
				x = width - 2;
				u = 0;
			}
			short s_minus_u = (short)((1. - u) * 2048); // 1-u
			short s_u = 2048 - s_minus_u;
			for (k = 0; k < channel; k++)
			{
				dst->imageData[i*dst_step + j*channel + k] = (srcData[y*src_step + x*channel + k] * s_minus_u * s_minus_v + \
					srcData[(y + 1)*src_step + x*channel + k] * s_minus_u * s_v + \
					srcData[y*src_step + (x + 1)*channel + k] * s_u * s_minus_v + \
					srcData[(y + 1) * src_step + (x + 1)*channel + k] * s_u * s_v) >> 22;
			}
		}

	}
	return dst;
}

ClImage* resize_nearest(ClImage* src, Size dsize, double fx, double fy)
{
	int width, height, channel;
	int i, j, k;
	getClImageInfo(src, &height, &width, &channel);
	if (dsize.height != 0 && dsize.width != 0)
	{
		fx = (double)dsize.width / (double)width;
		fy = (double)dsize.height / (double)height;
	}
	else if (fx != 0 && fy != 0)
	{
		dsize.width = (int)(width * fx);
		dsize.height = (int)(height * fy);
	}
	else
	{
		//
	}

	ClImage* dst = (ClImage*)malloc(sizeof(ClImage));
	dst->width = dsize.width;
	dst->height = dsize.height;
	dst->channels = channel;
	dst->imageData = (uchar*)malloc(sizeof(uchar)*dsize.width*dsize.height*channel);

	int dst_step = dsize.width * channel;
	int src_step = width * channel;

	for (i = 0; i < dsize.height; i++)
	{
		int srcy = floor(i / fy);
		srcy = min(srcy, height - 1);

		for (j = 0; j < dsize.width; j++)
		{
			int srcx = floor(j / fx);
			srcx = min(srcx, width - 1);
			for (k = 0; k < channel; k++)
			{
				dst->imageData[i*dst_step + j*channel + k] = src->imageData[srcy*src_step + srcx*channel + k];
			}
		}

	}

	return dst;
}


Mat* getGaussianKernel(int n, double sigma)
{
	static const float small_gaussian_tab[][7] =
	{
		{ 1.f },
		{ 0.25f, 0.5f, 0.25f },
		{ 0.0625f, 0.25f, 0.375f, 0.25f, 0.0625f },
		{ 0.03125f, 0.109375f, 0.21875f, 0.28125f, 0.21875f, 0.109375f, 0.03125f }
	};

	/*如果sigma小于0，且n为不大于7的奇整数，则核的滤波系数固定了，其固定在数组

	small_gaussian_tab中，根据其n的长度来选择具体的值 ，如果不满足上面的，则固定核为0
	固定核为0表示自己计算其核*/

	const float* fixed_kernel = n % 2 == 1 && n <= 7 && sigma <= 0 ?
		small_gaussian_tab[n >> 1] : 0;

	Mat* kernel = (Mat*)malloc(sizeof(Mat));
	kernel->height = n;
	kernel->width = n;
	kernel->channels = 1;
	kernel->imageData = (float*)malloc(sizeof(float) * n * n);
	float* cf = (float*)malloc(sizeof(float)*n);

	double sigmaX = sigma > 0 ? sigma : ((n - 1) * 0.5 - 1) * 0.3 + 0.8;//当sigma小于0时，采用公式得到sigma(只与n有关)
	double scale2X = -0.5 / (sigmaX*sigmaX);//高斯表达式后面要用到
	double sum = 0;

	int i;
	for (i = 0; i < n; i++)
	{
		double x = i - (n - 1) * 0.5;
		//如果自己算其核的话，就常用公式exp(scale2X*x*x)计算，否则就用固定系数的核
		double t = fixed_kernel ? (double)fixed_kernel[i] : exp(scale2X*x*x);
		cf[i] = (float)t;//单精度要求时存入cf数组中
		sum += cf[i];//进行归一化时要用到
	}

	sum = 1. / sum;//归一化时核中各元素之和为1
	for (i = 0; i < n; i++)
	{
		cf[i] = (float)(cf[i] * sum);//归一化后的单精度核元素
	}
	int j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			kernel->imageData[i*n + j] = cf[i] * cf[j];
		}
	}


	return kernel;//返回n*1的数组，其元素或是单精度或是双精度，且符合高斯分布
}


ClImage* GaussianBlur(ClImage* src, Size ksize, double sigma1, double sigma2)
{
	int height, width, channel;
	int i, j, k, x, y;
	getClImageInfo(src, &height, &width, &channel);

	if (ksize.width <= 0 && sigma1 > 0)
	{
		ksize.width = (int)round(sigma1 * 3 * 2 + 1);
		ksize.width = ksize.width | 1;
	}
		
	if (ksize.height <= 0 && sigma1 > 0)
	{
		ksize.height = (int)round(sigma2 * 3 * 2 + 1);
		ksize.height = ksize.height | 1;
	}
	sigma1 = max(sigma1, 0);
	sigma2 = max(sigma2, 0);

	Mat* kernel = getGaussianKernel(ksize.width, sigma2);
	int n = kernel->height;
	int mid = (n - 1) / 2;

	// apply gaussian blur to the image
	ClImage* dst = (ClImage*)malloc(sizeof(ClImage));
	dst->imageData = (uchar*)malloc(sizeof(uchar)*height*width*channel);
	dst->height = height;
	dst->width = width;
	dst->channels = channel;
	int step = width*channel;
	int kstep = n;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			for (k = 0; k < channel; k++)
			{
				float value = 0.;
				for (y = 0; y < n; y++)
				{
					for (x = 0; x < n; x++)
					{
						int p = i + (y - mid);
						int q = j + (x - mid);
						if (p < 0 || p >= height || q < 0 || q >= width)
							continue;
						value += (float)src->imageData[p*step + q*channel + k] *
							kernel->imageData[y*kstep + x];
					}
				}
				dst->imageData[i*step + j*channel + k] = (uchar)value;
			}
		}
	}
	return dst;
	
}


ClImage* combine(ClImage* img1, ClImage* img2)
{
	ClImage* img = (ClImage*)malloc(sizeof(ClImage));
	int h1, w1, c1, h2, w2, c2;
	getClImageInfo(img1, &h1, &w1, &c1); 
	getClImageInfo(img2, &h2, &w2, &c2);
	img->height = h1;
	img->width = (w1 + w2);
	img->channels = c1;
	img->imageData = (uchar*)malloc(sizeof(uchar)*h1*(w1 + w2)*c1);
	int i, j, k;
	int stepcomb = img->width*img->channels;
	int step1 = w1*c1;
	int step2 = w2*c2;

	for (i = 0; i < h1; i++)
	{
		for (j = 0; j < w1 + w2; j++)
		{
			for (k = 0; k < c1; k++)
			{
				if (j < w1)
				{
					img->imageData[i*stepcomb + j*c1 + k] = img1->imageData[i*step1 + j*c1 + k];
				}
				else
				{
					img->imageData[i*stepcomb + j*c1 + k] = img2->imageData[i*step2 + (j - w1)*c2 + k];
				}
			}
		}
	}
	return img;
}


void drawLine(ClImage* img, Point2d p1, Point2d p2, uchar* color)
{
	int step = img->width*img->channels;
	int x1 = p1.x, y1 = p1.y;
	int x2 = p2.x, y2 = p2.y;

	int dx = x2 - x1;
	int dy = y2 - y1;
	int ux = ((dx > 0) << 1) - 1;
	int uy = ((dy > 0) << 1) - 1;
	int x = x1, y = y1, eps;

	eps = 0;
	dx = abs(dx);
	dy = abs(dy);
	if (dx > dy)
	{
		for (x = x1; x != x2; x += ux)
		{
			img->imageData[y*step + x*img->channels] = color[0];
			img->imageData[y*step + x*img->channels+1] = color[1];
			img->imageData[y*step + x*img->channels+2] = color[2];
			eps += dy;
			if ((eps << 1) >= dx)
			{
				y += uy;
				eps -= dx;
			}
		}
		
	}
	else
	{
		for (y = y1; y != y2; y += uy)
		{
			img->imageData[y*step + x*img->channels] = color[0];
			img->imageData[y*step + x*img->channels + 1] = color[1];
			img->imageData[y*step + x*img->channels + 2] = color[2];
			eps += dx;
			if ((eps << 1) >= dy)
			{
				x += ux;
				eps -= dy;
			}
		}
		
	}
}

void write(KeyPoint* keypoints, Descriptor* descriptor, char* pt_path, char* descp_path)
{
	int num = descriptor->height, dim = descriptor->width;
	int i, j;
	FILE* pfile = fopen(pt_path, "w");
	for (i = 0; i < num; i++)
	{
		fprintf(pfile, "%d ", (int)(keypoints[i].pt.x));
		fprintf(pfile, "%d\n", (int)(keypoints[i].pt.y));
	}
	fclose(pfile);

	FILE* file = fopen(descp_path, "w");
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < dim; j++)
		{
			fprintf(file, "%.5f ", descriptor->imageData[i*dim + j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void draw(ClImage* img, Point2d point, uchar* color)
{
	int step = img->width*img->channels;
	int r = 2;
	int i;

	int x = point.x;
	int y = point.y;
	int p, q;
	for (p = y - r; p <= y + r; p++)
	{
		for (q = x - r; q <= x + r; q++)
		{
			img->imageData[p*step + q*img->channels] = color[0];
			img->imageData[p*step + q*img->channels + 1] = color[1];
			img->imageData[p*step + q*img->channels + 2] = color[2];
		}
	}

}

float L2dis(float* p1, float* p2, int n)
{
	int i;
	float dis = 0.f;
	for (i = 0; i < n; i++)
	{
		dis += (p2[i] - p1[i])*(p2[i] - p1[i]);
	}
	dis = sqrtf(dis);
	return dis;
}