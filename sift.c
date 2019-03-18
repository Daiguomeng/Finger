#include "ds.h"
#include "io.h"
#include "cvtools.h"
#include "math.h"
#include "string.h"

# define M_PI 3.14159265358979323846

// default number of sampled intervals per octave
static const int SIFT_INTVLS = 3;

// default sigma for initial gaussian smoothing
static const float SIFT_SIGMA = 1.6f;

// default threshold on keypoint contrast |D(x)|
static const float SIFT_CONTR_THR = 0.04f;

// default threshold on keypoint ratio of principle curvatures
static const float SIFT_CURV_THR = 10.f;

// double image size before pyramid construction?
static const bool SIFT_IMG_DBL = true;

// default width of descriptor histogram array
static const int SIFT_DESCR_WIDTH = 4;

// default number of bins per histogram in descriptor array
static const int SIFT_DESCR_HIST_BINS = 8;

// assumed gaussian blur for input image
static const float SIFT_INIT_SIGMA = 0.5f;

// width of border in which to ignore keypoints
static const int SIFT_IMG_BORDER = 5;

// maximum steps of keypoint interpolation before failure
static const int SIFT_MAX_INTERP_STEPS = 5;

// default number of bins in histogram for orientation assignment
static const int SIFT_ORI_HIST_BINS = 36;

// determines gaussian sigma for orientation assignment
static const float SIFT_ORI_SIG_FCTR = 1.5f;

// determines the radius of the region used in orientation assignment
static const float SIFT_ORI_RADIUS = 3 * 1.5f;

// orientation magnitude relative to max that results in new feature
static const float SIFT_ORI_PEAK_RATIO = 0.8f;

// determines the size of a single descriptor orientation histogram
static const float SIFT_DESCR_SCL_FCTR = 3.f;

// threshold on magnitude of elements of descriptor vector
static const float SIFT_DESCR_MAG_THR = 0.2f;

// factor used to convert floating-point descriptor to unsigned char
static const float SIFT_INT_DESCR_FCTR = 512.f;

static const int SIFT_FIXPT_SCALE = 1;

static const float FLT_EPSILON = 1e-7;

// SIFT parameters
int nfeatures = 0;
int nOctaveLayers = 3;
double contrastThreshold = 0.07;//0.04;
double edgeThreshold = 10;
double sigma = 1.6;

ClImage* createInitialImage(ClImage* img, bool doubleImageSize, float sigma)
{
	int height, width, channel;
	getClImageInfo(img, &height, &width, &channel);
	ClImage* gray = NULL;
	if (channel == 3)
		gray = BGR2Gray(img);
	else
		gray = img;
	float sig_diff;
	Size size = { 0, 0 };
	ClImage* init_img;
	if (doubleImageSize)
	{
		ClImage* dbl;
		Size target_size = { height * 2, width * 2 };
		sig_diff = sqrtf(max(sigma * sigma - SIFT_INIT_SIGMA * SIFT_INIT_SIGMA * 4, 0.01f));
		dbl = resize_bilinear(gray, target_size, 0, 0);
		init_img = GaussianBlur(dbl, size, sig_diff, sig_diff);
		return init_img;
	}
	else
	{
		sig_diff = sqrtf(max(sigma * sigma - SIFT_INIT_SIGMA * SIFT_INIT_SIGMA, 0.01f));
		init_img = GaussianBlur(gray, size, sig_diff, sig_diff);
		return init_img;
	}
}

ClImage** buildGaussianPyramid(ClImage* base, int nOctaves)
{
	double* sig = (double*)malloc(sizeof(double)*(nOctaveLayers + 3));
	sig[0] = sigma;
	double k = pow(2., 1. / nOctaveLayers);
	int i;
	for (i = 1; i < nOctaveLayers + 3; i++)
	{
		double sig_prev = pow(k, (double)(i - 1))*sigma;
		double sig_total = sig_prev*k;
		sig[i] = sqrt(sig_total*sig_total - sig_prev*sig_prev);
	}
	int o;
	ClImage** pyr = (ClImage**)malloc(sizeof(ClImage*)*nOctaves*(nOctaveLayers+3));
	for (i = 0; i < nOctaves*(nOctaveLayers + 3); i++)
		pyr[i] = (ClImage*)malloc(sizeof(ClImage));
	for (o = 0; o < nOctaves; o++)
	{
		for (i = 0; i < nOctaveLayers + 3; i++)
		{
			if (o == 0 && i == 0)
				pyr[o*(nOctaveLayers + 3) + i] = base;
			else if (i == 0)
			{
				ClImage* src = pyr[(o - 1)*(nOctaveLayers + 3) + nOctaveLayers];
				
				Size size = { (int)(src->height / 2), (int)(src->width / 2) };
				pyr[o*(nOctaveLayers + 3) + i] = resize_nearest(src, size, 0, 0);

			}
			else
			{
				ClImage* src = pyr[o*(nOctaveLayers + 3) + i - 1];
				Size size = { 0, 0 };
				pyr[o*(nOctaveLayers + 3) + i] = GaussianBlur(src, size, sig[i], sig[i]);
			}
				
		}
	}
	return pyr;
}

FtImage** buildDoGPyramid(ClImage** gpyr, int nOctaves)
{
	FtImage** dogpyr = (FtImage**)malloc(sizeof(FtImage*)*nOctaves*(nOctaveLayers + 2));
	int i, o, j;
	for (i = 0; i < nOctaves*(nOctaveLayers + 2); i++)
		dogpyr[i] = (FtImage*)malloc(sizeof(FtImage));
	for (o = 0; o < nOctaves; o++)
	{
		for (i = 0; i < nOctaveLayers+2; i++)
		{
			ClImage* src1 = gpyr[o*(nOctaveLayers + 3) + i];
			ClImage* src2 = gpyr[o*(nOctaveLayers + 3) + i + 1];
			dogpyr[o*(nOctaveLayers + 2) + i]->height = src1->height;
			dogpyr[o*(nOctaveLayers + 2) + i]->width = src1->width;
			dogpyr[o*(nOctaveLayers + 2) + i]->channels = src1->channels;
			dogpyr[o*(nOctaveLayers + 2) + i]->imageData = (float*)malloc(sizeof(float)*src1->height*src1->width*src1->channels);
			for (j = 0; j < src1->height*src1->width*src1->channels; j++)
			{
				/*int v = (int)(src2->imageData[j]) - (int)(src1->imageData[j]);
				if (v < 0)
					v = 0;*/
				dogpyr[o*(nOctaveLayers + 2) + i]->imageData[j] = (float)(src2->imageData[j]) - (float)(src1->imageData[j]);
			}
		}

	}
	return dogpyr;
}

int solve(float** H, float* b, float* X_new, int n)
{
	int r = n, c = n;
	int* loc = (int*)malloc(sizeof(int)*c);
	float* X = (float*)malloc(sizeof(float)*c);

	int i, j, k;
	for (i = 0; i < c; i++)
		loc[i] = i;

	for (i = 0; i < r - 1; i++)
	{
		float MAX = -INT_MAX;
		int p, q, mi, mj;
		for (p = i; p < n; p++)
		{
			for (q = i; q < n; q++)
			{
				if (MAX < fabs(*((float*)H+n*p+q)))
				{
					MAX = fabs(*((float*)H + n*p + q));
					mi = p;
					mj = q;
				}
			}
		}
		float tmp;
		// swap the row i and row mi
		for (j = 0; j < n; j++)
		{
			tmp = *((float*)H + n*i + j);
			*((float*)H + n*i + j) = *((float*)H + n*mi + j);
			*((float*)H + n*mi + j) = tmp;
		}
		tmp = b[i];
		b[i] = b[mi];
		b[mi] = tmp;

		// swap the column i and column mj
		for (j = 0; j < n; j++)
		{
			tmp = *((float*)H + n*j + i);
			*((float*)H + n*j + i) = *((float*)H + n*j + mj);
			*((float*)H + n*j + mj) = tmp;
		}
		tmp = loc[mj];
		loc[mj] = loc[i];
		loc[i] = tmp;

		if (*((float*)H + n*i + i) == 0)
			return -1;

		for (j = i + 1; j < n; j++)
		{
			float div = *((float*)H + n*j + i) / *((float*)H + n*i + i);
			for (k = i; k < n; k++)
				*((float*)H + n*j + k) -= *((float*)H + n*i + k) * div;
			b[j] -= b[i] * div;
		}
		
	}
	if (*((float*)H + n*(n-1) + (n-1)) == 0 && b[n - 1] != 0)
		return -1; // 无解
	else if (*((float*)H + n*(n - 1) + (n - 1)) == 0 && b[n - 1] == 0)
		return 0; // 无穷解
	else
	{
		for (i = n - 1; i >= 0; i--)
		{
			float sum = 0.;
			for (j = i + 1; j < n; j++)
				sum += X[j] * *((float*)H + n*i + j);
			X[i] = (b[i] - sum) / *((float*)H + n*i + i);
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (loc[j] == i)
				break;
		}
	
		X_new[i] = X[j];
	}
	return 1;
}

bool adjustLocalExtrema(FtImage** dog_pyr, KeyPoint* kpt, int octv,
	int* layer, int* r, int* c, int nOctaveLayers,
	float contrastThreshold, float edgeThreshold, float sigma)
{
	const float img_scale = 1.f / (255 * SIFT_FIXPT_SCALE);
	const float deriv_scale = img_scale*0.5f;
	const float second_deriv_scale = img_scale;
	const float cross_deriv_scale = img_scale*0.25f;

	float xi = 0, xr = 0, xc = 0, contr = 0;
	int i = 0;

	for (; i < SIFT_MAX_INTERP_STEPS; i++)
	{
		int idx = octv*(nOctaveLayers + 2) + (*layer);
		FtImage* img = dog_pyr[idx];
		FtImage* prev = dog_pyr[idx - 1];
		FtImage* next = dog_pyr[idx + 1]; 
		int step = img->width*img->channels;
		float dD[3] = { ((float)(img->imageData[(*r)*step + (*c) + 1]) - (float)(img->imageData[(*r)*step + (*c) - 1])) * deriv_scale,
			((float)(img->imageData[(*r + 1)*step + (*c)]) - (float)(img->imageData[(*r - 1)*step + (*c)])) * deriv_scale,
			((float)(next->imageData[(*r)*step + (*c)]) - (float)(prev->imageData[(*r)*step + (*c)])) * deriv_scale };
		float v2 = (float)img->imageData[(*r)*step + (*c)] * 2;
		float dxx = ((float)(img->imageData[(*r)*step + (*c) + 1]) + (float)(img->imageData[(*r)*step + (*c) - 1]) - v2) * second_deriv_scale;
		float dyy = ((float)(img->imageData[(*r + 1)*step + (*c)]) + (float)(img->imageData[(*r - 1)*step + (*c)]) - v2) * second_deriv_scale;
		float dss = ((float)(next->imageData[(*r)*step + (*c)]) + (float)(prev->imageData[(*r)*step + (*c)]) - v2) * second_deriv_scale;
		float dxy = ((float)(img->imageData[(*r + 1)*step + (*c) + 1]) - (float)(img->imageData[(*r + 1)*step + (*c) - 1]) -
			(float)(img->imageData[(*r - 1)*step + (*c) + 1]) + (float)(img->imageData[(*r - 1)*step + (*c) - 1])) * cross_deriv_scale;
		
		float dxs = ((float)(next->imageData[(*r)*step + (*c) + 1]) - (float)(next->imageData[(*r)*step + (*c) - 1]) -
			(float)(prev->imageData[(*r)*step + (*c) + 1]) + (float)(prev->imageData[(*r)*step + (*c) - 1])) * cross_deriv_scale;
	
		float dys = ((float)(next->imageData[(*r + 1)*step + (*c)]) - (float)(next->imageData[(*r - 1)*step + (*c)]) -
			(float)(prev->imageData[(*r + 1)*step + (*c)]) + (float)(prev->imageData[(*r - 1)*step + (*c)])) * cross_deriv_scale;
		
		float H[3][3] = { { dxx, dxy, dxs }, { dxy, dyy, dys }, { dxs, dys, dss } };
		// solve X;
		float* X = (float*)malloc(sizeof(float)*3);
		if (solve(H, dD, X, 3) > 0)
		{
			xi = -X[2];
			xr = -X[1];
			xc = -X[0];
		}
		else
		{
			xi = 0;
			xr = 0;
			xc = 0;
		}
		

		if (fabs(xi) < 0.5f && fabs(xr) < 0.5f && fabs(xc) < 0.5f)
			break;

		if (fabs(xi) > (float)(INT_MAX / 3) ||
			fabs(xr) > (float)(INT_MAX / 3) ||
			fabs(xc) > (float)(INT_MAX / 3))
			return false;

		(*c) += round(xc);
		(*r) += round(xr);
		(*layer) += round(xi);

		if (*layer < 1 || *layer > nOctaveLayers ||
			*c < SIFT_IMG_BORDER || *c >= img->width - SIFT_IMG_BORDER ||
			*r < SIFT_IMG_BORDER || *r >= img->height - SIFT_IMG_BORDER)
			return false;
	}
	if (i >= SIFT_MAX_INTERP_STEPS)
		return false;
	{
		int idx = octv*(nOctaveLayers + 2) + (*layer);
		FtImage* img = dog_pyr[idx];
		FtImage* prev = dog_pyr[idx - 1];
		FtImage* next = dog_pyr[idx + 1];
		int step = img->width*img->channels;
		float t = ((float)(img->imageData[(*r)*step + (*c) + 1]) - (float)(img->imageData[(*r)*step + (*c) - 1])) * deriv_scale * xc +
			((float)(img->imageData[(*r + 1)*step + (*c)]) - (float)(img->imageData[(*r - 1)*step + (*c)])) * deriv_scale * xr +
			((float)(next->imageData[(*r)*step + (*c)]) - (float)(prev->imageData[(*r)*step + (*c)])) * deriv_scale * xi;
		contr = img->imageData[(*r)*step + (*c)] * img_scale + t * 0.5f;
		if (fabs(contr) * nOctaveLayers < contrastThreshold)
			return false;

		float v2 = (float)img->imageData[(*r)*step + (*c)] * 2;
		float dxx = ((float)(img->imageData[(*r)*step + (*c) + 1]) + (float)(img->imageData[(*r)*step + (*c) - 1]) - v2) * second_deriv_scale;
		float dyy = ((float)(img->imageData[(*r + 1)*step + (*c)]) + (float)(img->imageData[(*r - 1)*step + (*c)]) - v2) * second_deriv_scale;
		float dxy = ((float)(img->imageData[(*r + 1)*step + (*c) + 1]) - (float)(img->imageData[(*r + 1)*step + (*c) - 1]) -
			(float)(img->imageData[(*r - 1)*step + (*c) + 1]) + (float)(img->imageData[(*r - 1)*step + (*c) - 1])) * cross_deriv_scale;
		float tr = dxx + dyy;
		float det = dxx * dyy - dxy * dxy;
		if (det <= 0 || tr*tr*edgeThreshold >= (edgeThreshold + 1)*(edgeThreshold + 1)*det)
			return false;

	}
	kpt->pt.x = (*c + xc) * (1 << octv);
	kpt->pt.y = (*r + xr)*(1 << octv);
	kpt->octave = octv + ((*layer) << 8) + ((int)(round((xi + 0.5) * 255)) << 16);
	kpt->size = sigma*powf(2.f, (*layer + xi) / nOctaveLayers)*(1 << octv) * 2;
	kpt->response = fabs(contr);
	return true;
}

float calcOrientationHist(ClImage* img, Point2d pt, int radius,
	float sigma, float* hist, int n)
{
	int i, j, k, len = (radius * 2 + 1)*(radius * 2 + 1);
	float expf_scale = -1.f / (2.f*sigma*sigma);
	float* buf = (float*)malloc(sizeof(float)*(len * 4 + n + 4));
	float *X = buf, *Y = X + len, *Mag = X, *Ori = Y + len, *W = Ori + len;
	float* temphist = W + len + 2;
	for (i = 0; i < n; i++)
		temphist[i] = 0.f;

	int step = img->width*img->channels;

	for (i = -radius, k = 0; i <= radius; i++)
	{
		int y = pt.y + i;
		if (y <= 0 || y >= img->height - 1)
			continue;
		for (j = -radius; j <= radius; j++)
		{
			int x = pt.x + j;
			if (x <= 0 || x >= img->width - 1)
				continue;
			float dx = (float)img->imageData[y*step + x + 1] - (float)img->imageData[y*step + x - 1];
			float dy = (float)img->imageData[(y - 1)*step + x] - (float)img->imageData[(y + 1)*step + x];
			X[k] = dx;
			Y[k] = dy;
			W[k] = (i*i + j*j)*expf_scale;
			k++;
		}
	}
	len = k;

	for (k = 0; k < len; k++)
		W[k] = exp(W[k]);
	for (k = 0; k < len; k++)
		Ori[k] = atan2f(Y[k], X[k]) * 180 / M_PI;
	for (k = 0; k < len; k++)
		Mag[k] = sqrtf(X[k] * X[k] + Y[k] * Y[k]);

	for (k = 0; k < len; k++)
	{
		int bin = round((n / 360.f)*Ori[k]);
		if (bin >= n)
			bin -= n;
		if (bin < 0)
			bin += n;
		temphist[bin] += W[k] * Mag[k];
	}
	temphist[-1] = temphist[n - 1];
	temphist[-2] = temphist[n - 2];
	temphist[n] = temphist[0];
	temphist[n + 1] = temphist[1];
	for (i = 0; i < n; i++)
	{
		hist[i] = (temphist[i - 2] + temphist[i + 2])*(1.f / 16.f) +
			(temphist[i - 1] + temphist[i + 1]) * (4.f / 16.f) +
			temphist[i] * (6.f / 16.f);

	}
	free(buf);
	float maxval = hist[0];
	for (i = 1; i < n; i++)
		maxval = max(maxval, hist[i]);
	return maxval;

}

void findScaleSpaceExtrema(ClImage** gauss_pyr, FtImage** dog_pyr, KeyPoint* keypoints, int* kpnum, int nOctaves)
{
	int threshold = floor(0.5*contrastThreshold / nOctaveLayers * 255 * SIFT_FIXPT_SCALE);
	const int n = SIFT_ORI_HIST_BINS;
	float* hist = (float*)malloc(sizeof(float)*n);
	KeyPoint kpt;
	int o, i, r, c;
	(*kpnum) = 0;
	for (o = 0; o < nOctaves; o++)
	{
		for (i = 1; i <= nOctaveLayers; i++)
		{
			int idx = o * (nOctaveLayers + 2) + i;
			FtImage* img = dog_pyr[idx];
			FtImage* prev = dog_pyr[idx - 1];
			FtImage* next = dog_pyr[idx + 1];
			int step = img->width*img->channels;
			int rows = img->height, cols = img->width;
			for (r = SIFT_IMG_BORDER; r < rows - SIFT_IMG_BORDER; r++)
			{
				for (c = SIFT_IMG_BORDER; c < cols - SIFT_IMG_BORDER; c++)
				{
					int cur = r*step + c;
					int l0 = cur - step - 1;
					int l1 = cur - step;
					int l2 = cur - step + 1;
					int l3 = cur - 1;
					int l4 = cur + 1;
					int l5 = cur + step - 1;
					int l6 = cur + step;
					int l7 = cur + step + 1;
					float val = img->imageData[cur];
					if (fabs(val) > threshold &&
						((val > 0 && val >= img->imageData[l0] && val >= img->imageData[l1] &&
						val >= img->imageData[l2] && val >= img->imageData[l3] && val >= img->imageData[l4] &&
						val >= img->imageData[l5] && val >= img->imageData[l6] && val >= img->imageData[l7] &&
						val >= prev->imageData[cur] && val >= prev->imageData[l0] && val >= prev->imageData[l1] &&
						val >= prev->imageData[l2] && val >= prev->imageData[l3] && val >= prev->imageData[l4] &&
						val >= prev->imageData[l5] && val >= prev->imageData[l6] && val >= prev->imageData[l7] &&
						val >= next->imageData[cur] && val >= next->imageData[l0] && val >= next->imageData[l1] &&
						val >= next->imageData[l2] && val >= next->imageData[l3] && val >= next->imageData[l4] &&
						val >= next->imageData[l5] && val >= next->imageData[l6] && val >= next->imageData[l7]) ||
						(val < 0 && val <= img->imageData[l0] && val <= img->imageData[l1] &&
						val <= img->imageData[l2] && val <= img->imageData[l3] && val <= img->imageData[l4] &&
						val <= img->imageData[l5] && val <= img->imageData[l6] && val <= img->imageData[l7] &&
						val <= prev->imageData[cur] && val <= prev->imageData[l0] && val <= prev->imageData[l1] &&
						val <= prev->imageData[l2] && val <= prev->imageData[l3] && val <= prev->imageData[l4] &&
						val <= prev->imageData[l5] && val <= prev->imageData[l6] && val <= prev->imageData[l7] &&
						val <= next->imageData[cur] && val <= next->imageData[l0] && val <= next->imageData[l1] &&
						val <= next->imageData[l2] && val <= next->imageData[l3] && val <= next->imageData[l4] &&
						val <= next->imageData[l5] && val <= next->imageData[l6] && val <= next->imageData[l7])))
					{
						int r1 = r, c1 = c, layer = i;
						if (!adjustLocalExtrema(dog_pyr, &kpt, o, &layer, &r1, &c1, nOctaveLayers,
							(float)contrastThreshold, (float)edgeThreshold, (float)sigma))
							continue;
						float scl_octv = kpt.size*0.5f / (1 << o);
						Point2d p = { c1, r1 };
						float omax = calcOrientationHist(gauss_pyr[o*(nOctaveLayers + 3) + layer],
							p, round(SIFT_ORI_RADIUS * scl_octv),
							SIFT_ORI_SIG_FCTR * scl_octv,
							hist, n);
						float mag_thr = (float)(omax * SIFT_ORI_PEAK_RATIO);
						int j;
						for (j = 0; j < n; j++)
						{
							int l = j >0 ? j - 1 : n - 1;
							int r2 = j < n - 1 ? j + 1 : 0;
							
							if (hist[j] >hist[l] && hist[j] > hist[r2] && hist[j] >= mag_thr)
							{
								float bin = j + 0.5f * (hist[l] - hist[r2]) / (hist[l] - 2 * hist[j] + hist[r2]);
								bin = bin < 0 ? n + bin : bin >= n ? bin - n : bin;
								kpt.angle = 360.f - (float)((360.f / n) * bin);
								if (fabs(kpt.angle - 360.f) < FLT_EPSILON)
									kpt.angle = 0.f;
								keypoints[*kpnum] = kpt;
								(*kpnum)++;
							}
						}
					}
				}
			}
		}
	}
}


void removeDuplicated(KeyPoint* keypoints, int* kpnum)
{
	int i, j;
	float thresh = 5;
	uchar* mark = (uchar*)malloc(sizeof(uchar)*(*kpnum));
	for (i = 0; i < (*kpnum); i++)
		mark[i] = (uchar)1;
	for (i = 0; i < *kpnum; i++)
	{
		if (mark[i])
		{
			for (j = i + 1; j < *kpnum; j++)
			{
				if (keypoints[i].pt.x == keypoints[j].pt.x && keypoints[i].pt.y == keypoints[j].pt.y &&
					keypoints[i].size == keypoints[j].size && keypoints[i].angle == keypoints[j].angle)
					mark[j] = (uchar)0;
			}
		}
	}
	int len = 0;
	for (i = 0; i < *kpnum; i++)
	{
		if (mark[i])
		{
			keypoints[len] = keypoints[i];
			len++;
		}
	}
	*kpnum = len;
}

void unpackOctave(KeyPoint kpt, int* octave, int* layer, float* scale)
{
	*octave = kpt.octave & 255;
	*layer = (kpt.octave >> 8) & 255;
	*octave = (*octave) < 128 ? (*octave) : (-128 | (*octave));
	*scale = (*octave)>0 ? 1.f / (1 << (*octave)) : (float)(1 << -(*octave));
}

void calcSIFTDescriptor(ClImage* img, Point2f ptf, float ori, float scl, int d, int n, float* dst, int ptid)
{
	Point2d pt = { round(ptf.x), round(ptf.y) };
	float cos_t = cosf(ori*(float)(M_PI / 180));
	float sin_t = sinf(ori*(float)(M_PI / 180));
	float bins_per_rad = n / 360.f;
	float exp_scale = -1.f / (d * d * 0.5f);
	float hist_width = SIFT_DESCR_SCL_FCTR * scl;
	int radius = round(hist_width * 1.4142135623730951f * (d + 1) * 0.5f);
	// Clip the radius to the diagonal of the image to avoid autobuffer too large exception
	radius = min(radius, (int)sqrt((double)img->width*img->width + img->height*img->height));
	cos_t /= hist_width;
	sin_t /= hist_width;

	int i, j, k, len = (radius * 2 + 1)*(radius * 2 + 1), histlen = (d + 2)*(d + 2)*(n + 2);
	int rows = img->height, cols = img->width;
	int step = img->width*img->channels;

	float* buf = (float*)malloc(sizeof(float)*(len * 6 + histlen));
	float *X = buf, *Y = X + len, *Mag = Y, *Ori = Mag + len, *W = Ori + len;
	float *RBin = W + len, *CBin = RBin + len, *hist = CBin + len;

	for (i = 0; i < d + 2; i++)
	{
		for (j = 0; j < d + 2; j++)
		for (k = 0; k < n + 2; k++)
			hist[(i*(d + 2) + j)*(n + 2) + k] = 0.;
	}

	for (i = -radius, k = 0; i <= radius; i++)
	for (j = -radius; j <= radius; j++)
	{
		// Calculate sample's histogram array coords rotated relative to ori.
		// Subtract 0.5 so samples that fall e.g. in the center of row 1 (i.e.
		// r_rot = 1.5) have full weight placed in row 1 after interpolation.
		float c_rot = j * cos_t - i * sin_t;
		float r_rot = j * sin_t + i * cos_t;
		float rbin = r_rot + d / 2 - 0.5f;
		float cbin = c_rot + d / 2 - 0.5f;
		int r = pt.y + i, c = pt.x + j;

		if (rbin > -1 && rbin < d && cbin > -1 && cbin < d &&
			r > 0 && r < rows - 1 && c > 0 && c < cols - 1)
		{
			float dx = (float)((float)img->imageData[r*step + c + 1] - (float)img->imageData[r*step + c - 1]);
			float dy = (float)((float)img->imageData[(r - 1)*step + c] - (float)img->imageData[(r + 1)*step + c]);
			X[k] = dx; Y[k] = dy; RBin[k] = rbin; CBin[k] = cbin;
			W[k] = (c_rot * c_rot + r_rot * r_rot)*exp_scale;
			k++;
		}
	}

	len = k;
	for (i = 0; i < len; i++)
	{
		Ori[i] = atan2f(Y[i], X[i]) * 180 / M_PI;
		Mag[i] = sqrtf(X[i] * X[i] + Y[i] * Y[i]);
		W[i] = expf(W[i]);
	}


	for (k = 0; k < len; k++)
	{
		float rbin = RBin[k], cbin = CBin[k];
		float obin = (Ori[k] - ori)*bins_per_rad;
		float mag = Mag[k] * W[k];

		int r0 = floor(rbin);
		int c0 = floor(cbin);
		int o0 = floor(obin);
		rbin -= r0;
		cbin -= c0;
		obin -= o0;

		if (o0 < 0)
			o0 += n;
		if (o0 >= n)
			o0 -= n;

		// histogram update using tri-linear interpolation
		float v_r1 = mag*rbin, v_r0 = mag - v_r1;
		float v_rc11 = v_r1*cbin, v_rc10 = v_r1 - v_rc11;
		float v_rc01 = v_r0*cbin, v_rc00 = v_r0 - v_rc01;
		float v_rco111 = v_rc11*obin, v_rco110 = v_rc11 - v_rco111;
		float v_rco101 = v_rc10*obin, v_rco100 = v_rc10 - v_rco101;
		float v_rco011 = v_rc01*obin, v_rco010 = v_rc01 - v_rco011;
		float v_rco001 = v_rc00*obin, v_rco000 = v_rc00 - v_rco001;

		int idx = ((r0 + 1)*(d + 2) + c0 + 1)*(n + 2) + o0;
		hist[idx] += v_rco000;
		hist[idx + 1] += v_rco001;
		hist[idx + (n + 2)] += v_rco010;
		hist[idx + (n + 3)] += v_rco011;
		hist[idx + (d + 2)*(n + 2)] += v_rco100;
		hist[idx + (d + 2)*(n + 2) + 1] += v_rco101;
		hist[idx + (d + 3)*(n + 2)] += v_rco110;
		hist[idx + (d + 3)*(n + 2) + 1] += v_rco111;
	}

	// finalize histogram, since the orientation histograms are circular
	for (i = 0; i < d; i++)
	for (j = 0; j < d; j++)
	{
		int idx = ((i + 1)*(d + 2) + (j + 1))*(n + 2);
		hist[idx] += hist[idx + n];
		hist[idx + 1] += hist[idx + n + 1];
		for (k = 0; k < n; k++)
			dst[(i*d + j)*n + k] = hist[idx + k];
	}

	free(buf);
	// copy histogram to the descriptor,
	// apply hysteresis thresholding
	// and scale the result, so that it can be easily converted
	// to byte array
	float nrm2 = 0;
	len = d*d*n;
	for (k = 0; k < len; k++)
		nrm2 += dst[k] * dst[k];
	float thr = sqrtf(nrm2)*SIFT_DESCR_MAG_THR;
	for (i = 0, nrm2 = 0; i < k; i++)
	{
		float val = min(dst[i], thr);
		dst[i] = val;
		nrm2 += val*val;
	}
	nrm2 = SIFT_INT_DESCR_FCTR / max(sqrtf(nrm2), FLT_EPSILON);

#if 1
	for (k = 0; k < len; k++)
	{
		dst[k] = (uchar)(dst[k] * nrm2);
	}
#else
	float nrm1 = 0;
	for (k = 0; k < len; k++)
	{
		dst[k] *= nrm2;
		nrm1 += dst[k];
	}
	nrm1 = 1.f / std::max(nrm1, FLT_EPSILON);
	for (k = 0; k < len; k++)
	{
		dst[k] = std::sqrt(dst[k] * nrm1);//saturate_cast<uchar>(std::sqrt(dst[k] * nrm1)*SIFT_INT_DESCR_FCTR);
	}
#endif
}

Descriptor* calcDescriptors(ClImage** gpyr, KeyPoint* keypoints, int kpnum, 
	int nOctaveLayers, int firstOctave)
{
	int d = SIFT_DESCR_WIDTH, n = SIFT_DESCR_HIST_BINS;
	Descriptor* descriptor = (Descriptor*)malloc(sizeof(Descriptor));
	descriptor->channels = 1;
	descriptor->height = kpnum;
	descriptor->width = SIFT_DESCR_WIDTH*SIFT_DESCR_WIDTH*SIFT_DESCR_HIST_BINS;
	descriptor->imageData = (float*)malloc(sizeof(float)*descriptor->height*descriptor->width*descriptor->channels);

	int i;
	for (i = 0; i < kpnum; i++)
	{
		KeyPoint kpt = keypoints[i];
		int octave, layer;
		float scale;
		unpackOctave(kpt, &octave, &layer, &scale);
		float size = kpt.size*scale;
		Point2f ptf = { kpt.pt.x*scale, kpt.pt.y*scale };
		ClImage* img = gpyr[(octave - firstOctave)*(nOctaveLayers + 3) + layer];
		float angle = 360.f - kpt.angle;
		if (fabs(angle - 360.f) < FLT_EPSILON)
			angle = 0.f;
		calcSIFTDescriptor(img, ptf, angle, size*0.5f, d, n, descriptor->imageData + i*descriptor->width, i);

	}
	return descriptor;
}


Descriptor* detect(ClImage* _image, KeyPoint* keypoints, bool useProvidedKeypoints)
{
	int firstOctave = -1, actualNOctaves = 0, actualNLayers = 0;
	int i;
	ClImage* base = createInitialImage(_image, firstOctave < 0, (float)sigma);
	int nOctaves = actualNOctaves > 0 ? actualNOctaves : round(log((double)min(base->width, base->height)) / log(2.) - 2) - firstOctave;

	ClImage** gpyr = buildGaussianPyramid(base, nOctaves);
	printf("Build gaussian pyramid... Done\n");
	FtImage** dogpyr = buildDoGPyramid(gpyr, nOctaves);
	printf("Build DoG pyramid... Done\n");

	int kpnum = 0;

	findScaleSpaceExtrema(gpyr, dogpyr, keypoints, &kpnum, nOctaves);
	removeDuplicated(keypoints, &kpnum);
	printf("Detect feature points... Done\n");
	printf("Get %d feature points.\n", kpnum);

	if (firstOctave < 0)
	{
		for (i = 0; i < kpnum; i++)
		{
			float scale = 1.f / (float)(1 << -firstOctave);
			keypoints[i].octave = (keypoints[i].octave&~255) | ((keypoints[i].octave + firstOctave) & 255);
			keypoints[i].pt.x *= scale;
			keypoints[i].pt.y *= scale;
			keypoints[i].size *= scale;
		}
	}

	Descriptor* descriptor = calcDescriptors(gpyr, keypoints, kpnum, nOctaveLayers, firstOctave);
	printf("Calculate descriptors... Done\n");

	return descriptor;
	

}


int match(ClImage* img1, ClImage* img2, Descriptor* descp1, Descriptor* descp2, KeyPoint* keypoints1, KeyPoint* keypoints2, int** matchp)
{
	int pnum1 = descp1->height, pnum2 = descp2->height;
	int dim = descp1->width;

	int i, j, p = -1, q = -1, step = descp1->width;
	float thresh = 0.4;
	int len = 0;

	for (i = 0; i < pnum1; i++)
	{
		float MIN = INT_MAX; // MAX = -INT_MAX;
		float sec_MIN = INT_MAX;
		for (j = 0; j < pnum2; j++)
		{

			float dis = L2dis(descp1->imageData + i*step, descp2->imageData + j*step, dim);
			if (MIN > dis) // 比最近的还要小
			{

				sec_MIN = MIN; // 将最小的赋给第二小的
				q = p;
				MIN = dis;
				p = j;
			}
			else if (dis < sec_MIN)
			{
				sec_MIN = dis;
				q = j;
			}
		}
		if (MIN / sec_MIN < thresh)
		{
			*((int*)matchp + 2 * len) = i;
			*((int*)matchp + 2 * len + 1) = p;
			len++;
		}

	}
	return len;
}

