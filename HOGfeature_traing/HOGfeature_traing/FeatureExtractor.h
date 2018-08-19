#ifndef FEATURE_EXTRACTOR_H
#define FEATURE_EXTRACTOR_H

#include <string>
#include <fstream>
#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include "PublicDataStructure.h"
#include "HogFeature.h"
#include "PrimalSVM.h"

using namespace std;
using namespace cv;

class FeatureExtractor
{
private:
	HogFeatureParameter* _hogFeature;
	string dir;	
	vector<vector<float>> descriptorValueList;
	vector<Mat> PosTrainingData;
	vector<Mat> NegTrainingData;
	HOGDescriptor *_descriptor;
	Ptr<ml::SVM> svm;

	Mat _trainingMat;
	
public:	
	FeatureExtractor();
	FeatureExtractor(string modelFile);
	~FeatureExtractor();
	Mat ShowHOGFeature(string fileName);
	Mat get_hogdescriptor_visu(const Mat & color_origImg, vector<float>& descriptorValues, const Size & size);
	
	void run();

	void loadTrainingData();
	void runTrainProcessOpenCV();
	void writeTrainingParameter(string fileName);
	
	void runTrainProcess();
	void extractSample();
	const HogFeatureParameter* getHogFeature();
	Mat ExtractorPositiveSample();
	Mat ExtractorNegativeSample();
	void generateHogTrainingFeature(vector<Mat>& inputsImages, float scale);
	Mat fillTrainingLabel();

	void KFoldCrossValidation(int k);
	void flipImage();

	void verifyNewHOG();

	void verifyHelmet();
	Rect checkROI(Mat frame, string type);	
	struct Head 
	{
		Point center;
		double radius;
	};
	bool detectedHeadHoughCircles(Mat grayFrame, Head * head);
};

#endif