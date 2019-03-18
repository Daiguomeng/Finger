#include "ds.h"
#include "io.h"

ClImage* createInitialImage(ClImage* img, bool doubleImageSize, float sigma);
ClImage** buildGaussianPyramid(ClImage* base, int nOctaves);
//ClImage** buildDoGPyramid(ClImage** gpyr, int nOctaves);
//void calcSIFTDescriptor(ClImage* img, Point2f ptf, float ori, float scl, int d, int n, float* dst);
int solve(float** H, float* b, float* X_new, int n);

Descriptor* detect(ClImage* _image, KeyPoint* keypoints, bool useProvidedKeypoints);
int match(ClImage* img1, ClImage* img2, Descriptor* descp1, Descriptor* descp2, KeyPoint* keypoints1, KeyPoint* keypoints2, int** matchp);
