#pragma once
#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
struct HogFeatureParameter
{
	cv::Size 	winSize;				//樣本大小
	cv::Size 	blockSize;				//Block大小
	cv::Size 	blockStride;			//掃描移動格數
	cv::Size 	cellSize;				//單元大小
	int 	nbins;					//直方圖bin的數量
	int 	derivAperture;		//聽說跟Sobel的Kernal有關
	double 	winSigma;				//高斯函數的方差
	int 	histogramNormType;		//歸一化的類型 
	double 	L2HysThreshold;		//Block内直方圖規一化類型L2-Hys的歸一化收缩率最大值為0.2
	bool 	gammaCorrection;		//是否Gamma校正
	int 	nlevels;				//樣本的最大数量
	bool 	signedGradient;		//梯度是否有號
};