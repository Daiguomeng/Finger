#include "io.h"


ClImage* clLoadImage(char* path)
{
	ClImage* bmpImg;
	FILE* pFile;
	unsigned short fileType;
	ClBitMapFileHeader bmpFileHeader;
	ClBitMapInfoHeader bmpInfoHeader;
	int channels = 1;
	int width = 0;
	int height = 0;
	int step = 0;
	int offset = 0;
	unsigned char pixVal;
	ClRgbQuad* quad;
	int i, j, k;

	bmpImg = (ClImage*)malloc(sizeof(ClImage));
	pFile = fopen(path, "rb");
	if (!pFile)
	{
		free(bmpImg);
		return NULL;
	}

	fread(&fileType, sizeof(unsigned short), 1, pFile);
	if (fileType == 0x4D42)
	{
		//printf("bmp file! \n");

		fread(&bmpFileHeader, sizeof(ClBitMapFileHeader), 1, pFile);
		/*printf("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");
		printf("bmp�ļ�ͷ��Ϣ��\n");
		printf("�ļ���С��%d \n", bmpFileHeader.bfSize);
		printf("�����֣�%d \n", bmpFileHeader.bfReserved1);
		printf("�����֣�%d \n", bmpFileHeader.bfReserved2);
		printf("λͼ����ƫ���ֽ�����%d \n", bmpFileHeader.bfOffBits);*/

		fread(&bmpInfoHeader, sizeof(ClBitMapInfoHeader), 1, pFile);
		/*printf("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");
		printf("bmp�ļ���Ϣͷ\n");
		printf("�ṹ�峤�ȣ�%d \n", bmpInfoHeader.biSize);
		printf("λͼ��ȣ�%d \n", bmpInfoHeader.biWidth);
		printf("λͼ�߶ȣ�%d \n", bmpInfoHeader.biHeight);
		printf("λͼƽ������%d \n", bmpInfoHeader.biPlanes);
		printf("��ɫλ����%d \n", bmpInfoHeader.biBitCount);
		printf("ѹ����ʽ��%d \n", bmpInfoHeader.biCompression);
		printf("ʵ��λͼ����ռ�õ��ֽ�����%d \n", bmpInfoHeader.biSizeImage);
		printf("X����ֱ��ʣ�%d \n", bmpInfoHeader.biXPelsPerMeter);
		printf("Y����ֱ��ʣ�%d \n", bmpInfoHeader.biYPelsPerMeter);
		printf("ʹ�õ���ɫ����%d \n", bmpInfoHeader.biClrUsed);
		printf("��Ҫ��ɫ����%d \n", bmpInfoHeader.biClrImportant);
		printf("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");*/

		if (bmpInfoHeader.biBitCount == 8)
		{
			//printf("���ļ��е�ɫ�壬����λͼΪ�����ɫ\n\n");
			channels = 1;
			width = bmpInfoHeader.biWidth;
			height = bmpInfoHeader.biHeight;
			offset = (channels*width) % 4;
			if (offset != 0)
			{
				offset = 4 - offset;
			}
			//bmpImg->mat = kzCreateMat(height, width, 1, 0);
			bmpImg->width = width;
			bmpImg->height = height;
			bmpImg->channels = 1;
			bmpImg->imageData = (unsigned char*)malloc(sizeof(unsigned char)*width*height);
			step = channels*width;

			quad = (ClRgbQuad*)malloc(sizeof(ClRgbQuad)* 256);
			fread(quad, sizeof(ClRgbQuad), 256, pFile);
			free(quad);

			for (i = 0; i<height; i++)
			{
				for (j = 0; j<width; j++)
				{
					fread(&pixVal, sizeof(unsigned char), 1, pFile);
					bmpImg->imageData[(height - 1 - i)*step + j] = pixVal;
				}
				if (offset != 0)
				{
					for (j = 0; j<offset; j++)
					{
						fread(&pixVal, sizeof(unsigned char), 1, pFile);
					}
				}
			}
		}
		else if (bmpInfoHeader.biBitCount == 24)
		{
			//printf("��λͼΪλ���ɫ\n\n");
			channels = 1;
			width = bmpInfoHeader.biWidth;
			height = bmpInfoHeader.biHeight;

			bmpImg->width = width;
			bmpImg->height = height;
			bmpImg->channels = 1;
			bmpImg->imageData = (unsigned char*)malloc(sizeof(unsigned char)*width * height);
			step = channels*width;

			offset = (channels*width) % 4;
			if (offset != 0)
			{
				offset = 4 - offset;
			}

			for (i = 0; i<height; i++)
			{
				for (j = 0; j<width; j++)
				{
					for (k = 0; k<3; k++)
					{
						fread(&pixVal, sizeof(unsigned char), 1, pFile);
						//bmpImg->imageData[(height - 1 - i)*step + j * 3 + k] = pixVal;
					}
					bmpImg->imageData[(height - 1 - i)*step + j] = pixVal;
					//kzSetMat(bmpImg->mat, height-1-i, j, kzScalar(pixVal[0], pixVal[1], pixVal[2]));
				}
				if (offset != 0)
				{
					for (j = 0; j<offset; j++)
					{
						fread(&pixVal, sizeof(unsigned char), 1, pFile);
					}
				}
			}
		}
	}

	return bmpImg;
}

bool clSaveImage(char* path, ClImage* bmpImg)
{
	FILE *pFile;
	unsigned short fileType;
	ClBitMapFileHeader bmpFileHeader;
	ClBitMapInfoHeader bmpInfoHeader;
	int step;
	int offset;
	unsigned char pixVal = '\0';
	int i, j;
	ClRgbQuad* quad;

	pFile = fopen(path, "wb");
	if (!pFile)
	{
		return false;
	}

	fileType = 0x4D42;
	fwrite(&fileType, sizeof(unsigned short), 1, pFile);

	if (bmpImg->channels == 3)//24λ��ͨ������ͼ
	{
		step = bmpImg->channels*bmpImg->width;
		offset = step % 4;
		if (offset != 4)
		{
			step += 4 - offset;
		}

		bmpFileHeader.bfSize = bmpImg->height*step + 54;
		bmpFileHeader.bfReserved1 = 0;
		bmpFileHeader.bfReserved2 = 0;
		bmpFileHeader.bfOffBits = 54;
		fwrite(&bmpFileHeader, sizeof(ClBitMapFileHeader), 1, pFile);

		bmpInfoHeader.biSize = 40;
		bmpInfoHeader.biWidth = bmpImg->width;
		bmpInfoHeader.biHeight = bmpImg->height;
		bmpInfoHeader.biPlanes = 1;
		bmpInfoHeader.biBitCount = 24;
		bmpInfoHeader.biCompression = 0;
		bmpInfoHeader.biSizeImage = bmpImg->height*step;
		bmpInfoHeader.biXPelsPerMeter = 0;
		bmpInfoHeader.biYPelsPerMeter = 0;
		bmpInfoHeader.biClrUsed = 0;
		bmpInfoHeader.biClrImportant = 0;
		fwrite(&bmpInfoHeader, sizeof(ClBitMapInfoHeader), 1, pFile);

		for (i = bmpImg->height - 1; i>-1; i--)
		{
			for (j = 0; j<bmpImg->width; j++)
			{
				pixVal = bmpImg->imageData[i*bmpImg->width * 3 + j * 3];
				fwrite(&pixVal, sizeof(unsigned char), 1, pFile);
				pixVal = bmpImg->imageData[i*bmpImg->width * 3 + j * 3 + 1];
				fwrite(&pixVal, sizeof(unsigned char), 1, pFile);
				pixVal = bmpImg->imageData[i*bmpImg->width * 3 + j * 3 + 2];
				fwrite(&pixVal, sizeof(unsigned char), 1, pFile);
			}
			if (offset != 0)
			{
				for (j = 0; j<offset; j++)
				{
					pixVal = 0;
					fwrite(&pixVal, sizeof(unsigned char), 1, pFile);
				}
			}
		}
	}
	else if (bmpImg->channels == 1)//8λ����ͨ�����Ҷ�ͼ
	{
		step = bmpImg->width;
		offset = step % 4;
		if (offset != 4)
		{
			step += 4 - offset;
		}

		bmpFileHeader.bfSize = 54 + 256 * 4 + bmpImg->width;
		bmpFileHeader.bfReserved1 = 0;
		bmpFileHeader.bfReserved2 = 0;
		bmpFileHeader.bfOffBits = 54 + 256 * 4;
		fwrite(&bmpFileHeader, sizeof(ClBitMapFileHeader), 1, pFile);

		bmpInfoHeader.biSize = 40;
		bmpInfoHeader.biWidth = bmpImg->width;
		bmpInfoHeader.biHeight = bmpImg->height;
		bmpInfoHeader.biPlanes = 1;
		bmpInfoHeader.biBitCount = 8;
		bmpInfoHeader.biCompression = 0;
		bmpInfoHeader.biSizeImage = bmpImg->height*step;
		bmpInfoHeader.biXPelsPerMeter = 0;
		bmpInfoHeader.biYPelsPerMeter = 0;
		bmpInfoHeader.biClrUsed = 256;
		bmpInfoHeader.biClrImportant = 256;
		fwrite(&bmpInfoHeader, sizeof(ClBitMapInfoHeader), 1, pFile);

		quad = (ClRgbQuad*)malloc(sizeof(ClRgbQuad)* 256);
		for (i = 0; i<256; i++)
		{
			quad[i].rgbBlue = i;
			quad[i].rgbGreen = i;
			quad[i].rgbRed = i;
			quad[i].rgbReserved = 0;
		}
		fwrite(quad, sizeof(ClRgbQuad), 256, pFile);
		free(quad);

		for (i = bmpImg->height - 1; i>-1; i--)
		{
			for (j = 0; j<bmpImg->width; j++)
			{
				pixVal = bmpImg->imageData[i*bmpImg->width + j];
				fwrite(&pixVal, sizeof(unsigned char), 1, pFile);
			}
			if (offset != 0)
			{
				for (j = 0; j<offset; j++)
				{
					pixVal = 0;
					fwrite(&pixVal, sizeof(unsigned char), 1, pFile);
				}
			}
		}
	}
	fclose(pFile);

	return true;
}

