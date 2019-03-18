#ifndef _DS_H_
#define _DS_H_

typedef struct
{
	int x;
	int y;
}Point2d;

typedef struct
{
	float x;
	float y;
}Point2f;

typedef struct
{
	Point2f pt;
	float size;
	float angle;
	float response;
	int octave;
	int class_id;

}KeyPoint;


typedef struct
{
	int height;
	int width;
}Size;

typedef unsigned char uchar;



#endif