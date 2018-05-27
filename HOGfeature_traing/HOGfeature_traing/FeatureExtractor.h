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
	const HogFeatureParameter* getHogFeature();
	Mat ExtractorPositiveSample();
	Mat ExtractorNegativeSample();
	Mat ShowHOGFeature(string fileName);
	Mat get_hogdescriptor_visu(const Mat & color_origImg, vector<float>& descriptorValues, const Size & size);	
	void loadTrainingData(string posPath, string negPath);
	void extractSample();
	void run();
	void runTrainProcess();
	void runTrainProcessOpenCV();
	void writeTrainingParameter(string fileName);
	void test();
	void generateHogTrainingFeature(vector<Mat>& inputsImages, float scale);
	Mat fillTrainingLabel();
};

#endif