/*
	FAST-EDGE
	Copyright (c) 2009 Benjamin C. Haynor

	Permission is hereby granted, free of charge, to any person
	obtaining a copy of this software and associated documentation
	files (the "Software"), to deal in the Software without
	restriction, including without limitation the rights to use,
	copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following
	conditions:

	The above copyright notice and this permission notice shall be
	included in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
	OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FASTEDGE
#define _FASTEDGE
#define LOW_THRESHOLD_PERCENTAGE 0.8 // percentage of the high threshold value that the low threshold shall be set at
#define PI 3.14159265
#define HIGH_THRESHOLD_PERCENTAGE 0.20 // percentage of pixels that meet the high threshold - for example 0.15 will ensure that at least 15% of edge pixels are considered to meet the high threshold

#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) < (Y) ? (Y) : (X))

//#define WIDTH 640			// uncomment to define width for situations where width is always known
//#define HEIGHT 480		// uncomment to define heigh for situations where height is always known

//#define CLOCK			// uncomment to show running times of image processing functions (in seconds)
//#define ABS_APPROX		// uncomment to use the absolute value approximation of sqrt(Gx ^ 2 + Gy ^2)
//#define PRINT_HISTOGRAM	// uncomment to print the histogram used to estimate the threshold

ClImage* canny_edge_detect(ClImage * img_int);
void gaussian_noise_reduce(ClImage * img_in, ClImage * img_out);
void calc_gradient_sobel(ClImage * img_in, int g[], int dir[]);
void calc_gradient_scharr(ClImage * img_in, int g_x[], int g_y[], int g[], int dir[]);
void non_max_suppression(ClImage * img, int g[], int dir[]);
void estimate_threshold(ClImage * img, int * high, int * low);
void hysteresis (int high, int low, ClImage * img_in, ClImage * img_out);
int trace (int x, int y, int low, ClImage * img_in, ClImage * img_out);
int range (ClImage * img, int x, int y);
void dilate_1d_h(ClImage * img, ClImage * img_out);
void dilate_1d_v(ClImage * img, ClImage * img_out);
void erode_1d_h(ClImage * img, ClImage * img_out);
void erode_1d_v(ClImage * img, ClImage * img_out);
void erode(ClImage * img_in, ClImage * img_scratch, ClImage * img_out);
void dilate(ClImage * img_in, ClImage * img_scratch, ClImage * img_out);
void morph_open(ClImage * img_in, ClImage * img_scratch, ClImage * img_scratch2, ClImage * img_out);
void morph_close(ClImage * img_in, ClImage * img_scratch, ClImage * img_scratch2, ClImage * img_out);
#endif
