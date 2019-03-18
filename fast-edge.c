#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "io.h"
#include "fast-edge.h"

#ifdef WIDTH
int g[WIDTH  * HEIGHT], dir[WIDTH  * HEIGHT] = {0};
unsigned char img_scratch_data[WIDTH  * HEIGHT] = {0};
#endif
ClImage* canny_edge_detect(ClImage * img_in) {
	ClImage* img_out = (ClImage*)malloc(sizeof(ClImage));
	img_out->height = img_in->height;
	img_out->width = img_in->width;
	img_out->channels = 1;
	img_out->imageData = malloc(img_in->height*img_in->width * sizeof(unsigned char));
	ClImage img_scratch;
	int high, low;
	#ifndef WIDTH
	int * g = (int *)malloc(img_in->width * img_in->height * sizeof(int));
	int * dir = malloc(img_in->width * img_in->height * sizeof(int));
	unsigned char * img_scratch_data = malloc(img_in->width * img_in->height * sizeof(unsigned char));
	#endif
	img_scratch.width = img_in->width;
	img_scratch.height = img_in->height;
	img_scratch.imageData = img_scratch_data;
	calc_gradient_sobel(img_in, g, dir);
	printf("*** performing non-maximum suppression ***\n");
	non_max_suppression(&img_scratch, g, dir);
	estimate_threshold(&img_scratch, &high, &low);
	hysteresis(high, low, &img_scratch, img_out);
	#ifndef WIDTH
	free(g);
	free(dir);
	free(img_scratch_data);
	#endif
	return img_out;
}
void gaussian_noise_reduce(ClImage * img_in, ClImage * img_out)
{
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	int w, h, x, y, max_x, max_y;
	w = img_in->width;
	h = img_in->height;
	img_out->width = w;
	img_out->height = h;
	max_x = w - 2;
	max_y = w * (h - 2);
	for (y = w * 2; y < max_y; y += w) {
		for (x = 2; x < max_x; x++) {
			img_out->imageData[x + y] = (2 * img_in->imageData[x + y - 2 - w - w] + 
			4 * img_in->imageData[x + y - 1 - w - w] + 
			5 * img_in->imageData[x + y - w - w] + 
			4 * img_in->imageData[x + y + 1 - w - w] + 
			2 * img_in->imageData[x + y + 2 - w - w] + 
			4 * img_in->imageData[x + y - 2 - w] + 
			9 * img_in->imageData[x + y - 1 - w] + 
			12 * img_in->imageData[x + y - w] + 
			9 * img_in->imageData[x + y + 1 - w] + 
			4 * img_in->imageData[x + y + 2 - w] + 
			5 * img_in->imageData[x + y - 2] + 
			12 * img_in->imageData[x + y - 1] + 
			15 * img_in->imageData[x + y] + 
			12 * img_in->imageData[x + y + 1] + 
			5 * img_in->imageData[x + y + 2] + 
			4 * img_in->imageData[x + y - 2 + w] + 
			9 * img_in->imageData[x + y - 1 + w] + 
			12 * img_in->imageData[x + y + w] + 
			9 * img_in->imageData[x + y + 1 + w] + 
			4 * img_in->imageData[x + y + 2 + w] + 
			2 * img_in->imageData[x + y - 2 + w + w] + 
			4 * img_in->imageData[x + y - 1 + w + w] + 
			5 * img_in->imageData[x + y + w + w] + 
			4 * img_in->imageData[x + y + 1 + w + w] + 
			2 * img_in->imageData[x + y + 2 + w + w]) / 159;
		}
	}
	#ifdef CLOCK
	printf("Gaussian noise reduction - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}

void calc_gradient_sobel(ClImage * img_in, int g[], int dir[]) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	int w, h, x, y, max_x, max_y, g_x, g_y;
	float g_div;
	w = img_in->width;
	h = img_in->height;
	max_x = w - 3;
	max_y = w * (h - 3);
	for (y = w * 3; y < max_y; y += w) {
		for (x = 3; x < max_x; x++) {
			g_x = (2 * img_in->imageData[x + y + 1] 
				+ img_in->imageData[x + y - w + 1]
				+ img_in->imageData[x + y + w + 1]
				- 2 * img_in->imageData[x + y - 1] 
				- img_in->imageData[x + y - w - 1]
				- img_in->imageData[x + y + w - 1]);
			g_y = 2 * img_in->imageData[x + y - w] 
				+ img_in->imageData[x + y - w + 1]
				+ img_in->imageData[x + y - w - 1]
				- 2 * img_in->imageData[x + y + w] 
				- img_in->imageData[x + y + w + 1]
				- img_in->imageData[x + y + w - 1];
			#ifndef ABS_APPROX
			g[x + y] = sqrt(g_x * g_x + g_y * g_y);
			#endif
			#ifdef ABS_APPROX
			g[x + y] = abs(g_x[x + y]) + abs(g_y[x + y]);
			#endif
			if (g_x == 0) {
				dir[x + y] = 2;
			} else {
				g_div = g_y / (float) g_x;
				/* the following commented-out code is slightly faster than the code that follows, but is a slightly worse approximation for determining the edge direction angle
				if (g_div < 0) {
					if (g_div < -1) {
						dir[n] = 0;
					} else {
						dir[n] = 1;
					}
				} else {
					if (g_div > 1) {
						dir[n] = 0;
					} else {
						dir[n] = 3;
					}
				}
				*/
				if (g_div < 0) {
					if (g_div < -2.41421356237) {
						dir[x + y] = 0;
					} else {
						if (g_div < -0.414213562373) {
							dir[x + y] = 1;
						} else {
							dir[x + y] = 2;
						}
					}
				} else {
					if (g_div > 2.41421356237) {
						dir[x + y] = 0;
					} else {
						if (g_div > 0.414213562373) {
							dir[x + y] = 3;
						} else {
							dir[x + y] = 2;
						}
					}
				}
			}
		}
		
	}	
	#ifdef CLOCK
	printf("Calculate gradient Sobel - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}

/*
	CALC_GRADIENT_SCHARR
	calculates the result of the Scharr version of the Sobel operator - http://en.wikipedia.org/wiki/Sobel_operator - and estimates edge direction angle
	may have better rotational symmetry
*/
void calc_gradient_scharr(ClImage * img_in, int g_x[], int g_y[], int g[], int dir[]) {//float theta[]) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	int w, h, x, y, max_x, max_y, n;
	float g_div;
	w = img_in->width;
	h = img_in->height;
	max_x = w - 1;
	max_y = w * (h - 1);
	n = 0;
	for (y = w; y < max_y; y += w) {
		for (x = 1; x < max_x; x++) {
			g_x[n] = (10 * img_in->imageData[x + y + 1] 
				+ 3 * img_in->imageData[x + y - w + 1]
				+ 3 * img_in->imageData[x + y + w + 1]
				- 10 * img_in->imageData[x + y - 1] 
				- 3 * img_in->imageData[x + y - w - 1]
				- 3 * img_in->imageData[x + y + w - 1]);
			g_y[n] = 10 * img_in->imageData[x + y - w] 
				+ 3 * img_in->imageData[x + y - w + 1]
				+ 3 * img_in->imageData[x + y - w - 1]
				- 10 * img_in->imageData[x + y + w] 
				- 3 * img_in->imageData[x + y + w + 1]
				- 3 * img_in->imageData[x + y + w - 1];
			#ifndef ABS_APPROX
			g[n] = sqrt(g_x[n] * g_x[n] + g_y[n] * g_y[n]);
			#endif
			#ifdef ABS_APPROX
			g[n] = abs(g_x[n]) + abs(g_y[n]);
			#endif
			if (g_x[n] == 0) {
				dir[n] = 2;
			} else {
				g_div = g_y[n] / (float) g_x[n];
				if (g_div < 0) {
					if (g_div < -2.41421356237) {
						dir[n] = 0;
					} else {
						if (g_div < -0.414213562373) {
							dir[n] = 1;
						} else {
							dir[n] = 2;
						}
					}
				} else {
					if (g_div > 2.41421356237) {
						dir[n] = 0;
					} else {
						if (g_div > 0.414213562373) {
							dir[n] = 3;
						} else {
							dir[n] = 2;
						}
					}
				}
			}
			n++;
		}
	}	
	#ifdef CLOCK
	printf("Calculate gradient Scharr - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}
/*
	NON_MAX_SUPPRESSION
	using the estimates of the Gx and Gy image gradients and the edge direction angle determines whether the magnitude of the gradient assumes a local  maximum in the gradient direction
	if the rounded edge direction angle is 0 degrees, checks the north and south directions
	if the rounded edge direction angle is 45 degrees, checks the northwest and southeast directions
	if the rounded edge direction angle is 90 degrees, checks the east and west directions
	if the rounded edge direction angle is 135 degrees, checks the northeast and southwest directions
*/
void non_max_suppression(ClImage * img, int g[], int dir[]) {//float theta[]) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	int w, h, x, y, max_x, max_y;
	w = img->width;
	h = img->height;
	max_x = w;
	max_y = w * h;
	for (y = 0; y < max_y; y += w) {
		for (x = 0; x < max_x; x++) {
			switch (dir[x + y]) {
				case 0:
					if (g[x + y] > g[x + y - w] && g[x + y] > g[x + y + w]) {
						if (g[x + y] > 255) {
						img->imageData[x + y] = 0xFF;
						} else {
							img->imageData[x + y] = g[x + y];
						}
					} else {
						img->imageData[x + y] = 0x00;
					}
					break;
				case 1:
					if (g[x + y] > g[x + y - w - 1] && g[x + y] > g[x + y + w + 1]) {
						if (g[x + y] > 255) {
						img->imageData[x + y] = 0xFF;
						} else {
							img->imageData[x + y] = g[x + y];
						}
					} else {
						img->imageData[x + y] = 0x00;
					}
					break;
				case 2:
					if (g[x + y] > g[x + y - 1] && g[x + y] > g[x + y + 1]) {
						if (g[x + y] > 255) {
						img->imageData[x + y] = 0xFF;
						} else {
							img->imageData[x + y] = g[x + y];
						}
					} else { 
						img->imageData[x + y] = 0x00;
					}
					break;
				case 3:
					if (g[x + y] > g[x + y - w + 1] && g[x + y] > g[x + y + w - 1]) {
						if (g[x + y] > 255) {
						img->imageData[x + y] = 0xFF;
						} else {
							img->imageData[x + y] = g[x + y];
						}
					} else {
						img->imageData[x + y] = 0x00;
					}
					break;
				default:
					printf("ERROR - direction outside range 0 to 3");
					break;
			}
		}
	}
	#ifdef CLOCK
	printf("Non-maximum suppression - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}
/*
	ESTIMATE_THRESHOLD
	estimates hysteresis threshold, assuming that the top X% (as defined by the HIGH_THRESHOLD_PERCENTAGE) of edge pixels with the greatest intesity are true edges
	and that the low threshold is equal to the quantity of the high threshold plus the total number of 0s at the low end of the histogram divided by 2
*/
void estimate_threshold(ClImage * img, int * high, int * low) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	int i, max, pixels, high_cutoff;
	int histogram[256];
	max = img->width * img->height;
	for (i = 0; i < 256; i++) {
		histogram[i] = 0;
	}
	for (i = 0; i < max; i++) {
		histogram[img->imageData[i]]++;
	}
	pixels = (max - histogram[0]) * HIGH_THRESHOLD_PERCENTAGE;
	high_cutoff = 0;
	i = 255;
	while (high_cutoff < pixels) {
		high_cutoff += histogram[i];
		i--;
	}
	*high = i;
	i = 1;
	while (histogram[i] == 0) {
		i++;
	}
	*low = (*high + i) * LOW_THRESHOLD_PERCENTAGE;
	#ifdef PRINT_HISTOGRAM
	for (i = 0; i < 256; i++) {
		printf("i %d count %d\n", i, histogram[i]);
	}
	#endif
	
	#ifdef CLOCK
	printf("Estimate threshold - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}

void hysteresis (int high, int low, ClImage * img_in, ClImage * img_out)
{
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	int x, y, n, max;
	max = img_in->width * img_in->height;
	for (n = 0; n < max; n++) {
		img_out->imageData[n] = 0x00;
	}
	for (y=0; y < img_out->height; y++) {
	  for (x=0; x < img_out->width; x++) {
			if (img_in->imageData[y * img_out->width + x] >= high) {
				trace (x, y, low, img_in, img_out);
			}
		}
	}
	#ifdef CLOCK
	printf("Hysteresis - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}

int trace(int x, int y, int low, ClImage * img_in, ClImage * img_out)
{
	int y_off, x_off;//, flag;
	if (img_out->imageData[y * img_out->width + x] == 0)
	{
		img_out->imageData[y * img_out->width + x] = 255;
		for (y_off = -1; y_off <=1; y_off++)
		{
		    for(x_off = -1; x_off <= 1; x_off++)
		    {
				if (!(y == 0 && x_off == 0) && range(img_in, x + x_off, y + y_off) && img_in->imageData[(y + y_off) * img_out->width + x + x_off] >= low) {
					if (trace(x + x_off, y + y_off, low, img_in, img_out))
					{
					    return(1);
					}
				}
		    }
		}
		return(1);
	}
	return(0);
}

int range(ClImage * img, int x, int y)
{
	if ((x < 0) || (x >= img->width)) {
		return(0);
	}
	if ((y < 0) || (y >= img->height)) {
		return(0);
	}
	return(1);
}

void dilate_1d_h(ClImage * img, ClImage * img_out) {
	int x, y, offset, y_max;
	y_max = img->height * (img->width - 2);
	for (y = 2 * img->width; y < y_max; y += img->width) {
		for (x = 2; x < img->width - 2; x++) {
			offset = x + y;
			img_out->imageData[offset] = max(max(max(max(img->imageData[offset-2], img->imageData[offset-1]), img->imageData[offset]), img->imageData[offset+1]), img->imageData[offset+2]);	
		}
	}
}

void dilate_1d_v(ClImage * img, ClImage * img_out) {
	int x, y, offset, y_max;
	y_max = img->height * (img->width - 2);
	for (y = 2 * img->width; y < y_max; y += img->width) {
		for (x = 2; x < img->width - 2; x++) {
			offset = x + y;
			img_out->imageData[offset] = max(max(max(max(img->imageData[offset-2 * img->width], img->imageData[offset-img->width]), img->imageData[offset]), img->imageData[offset+img->width]), img->imageData[offset+2*img->width]);	
		}
	}
}

void erode_1d_h(ClImage * img, ClImage * img_out) {
	int x, y, offset, y_max;
	y_max = img->height * (img->width - 2);
	for (y = 2 * img->width; y < y_max; y += img->width) {
		for (x = 2; x < img->width - 2; x++) {
			offset = x + y;
			img_out->imageData[offset] = min(min(min(min(img->imageData[offset-2], img->imageData[offset-1]), img->imageData[offset]), img->imageData[offset+1]), img->imageData[offset+2]);	
		}
	}
}

void erode_1d_v(ClImage * img, ClImage * img_out) {
	int x, y, offset, y_max;
	y_max = img->height * (img->width - 2);
	for (y = 2 * img->width; y < y_max; y += img->width) {
		for (x = 2; x < img->width - 2; x++) {
			offset = x + y;
			img_out->imageData[offset] = min(min(min(min(img->imageData[offset-2 * img->width], img->imageData[offset-img->width]), img->imageData[offset]), img->imageData[offset+img->width]), img->imageData[offset+2*img->width]);	
		}
	}
}

void erode(ClImage * img_in, ClImage * img_scratch, ClImage * img_out) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	erode_1d_h(img_in, img_scratch);
	erode_1d_v(img_scratch, img_out);
	#ifdef CLOCK
	printf("Erosion - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}

void dilate(ClImage * img_in, ClImage * img_scratch, ClImage * img_out) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	dilate_1d_h(img_in, img_scratch);
	dilate_1d_v(img_scratch, img_out);
	#ifdef CLOCK
	printf("Dilation - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}

void morph_open(ClImage * img_in, ClImage * img_scratch, ClImage * img_scratch2, ClImage * img_out) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	erode(img_in, img_scratch, img_scratch2);
	dilate(img_scratch2, img_scratch, img_out);
	#ifdef CLOCK
	printf("Morphological opening - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}

void morph_close(ClImage * img_in, ClImage * img_scratch, ClImage * img_scratch2, ClImage * img_out) {
	#ifdef CLOCK
	clock_t start = clock();
	#endif
	dilate(img_in, img_scratch, img_scratch2);
	erode(img_scratch2, img_scratch, img_out);
	#ifdef CLOCK
	printf("Morphological closing - time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	#endif
}
