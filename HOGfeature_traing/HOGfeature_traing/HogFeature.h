#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>


#include "PublicDataStructure.h"

using namespace std;
using namespace cv;

class HogFeature
{
	public:
		HogFeature();
		HogFeature(HogFeatureParameter *hogFeatureParameter);
		~HogFeature();
		Mat calculateResidualedEdgeImage(Mat &inputMat);
		void calculateOrientation(Mat &x_Residuals, Mat &y_Residuals, vector<float> &hogFeatureVector, int &blockX, int &blockY, int &vectorIndex);
		void calculateHogFeature(Mat &inputMat, vector<float> &returnHogFeature);
		void getWindowBlocks(Mat &frame, vector<float> &hogFeatures, int &index_X, int &index_Y, vector<float> &returnFeature, float &trainScale);
	private:
		HogFeatureParameter  *_hogFeatureParameter;
		float _radiusToDegree;
		int cellsInBlock = 2;

		float squareRoot(float x);
		float myATan2(int y, int x);

		void calculateResiduale_X(Mat &input, Mat &outputMat);
		void calculateResiduale_Y(Mat &input, Mat &outputMat);
		int pixelSubtract(uchar firstPixel, uchar secPixel);
		void normalizeHogFeature(Mat &inputMat, vector<float> &hogFeatureVector, vector<float> &normalizedFeature);
};

