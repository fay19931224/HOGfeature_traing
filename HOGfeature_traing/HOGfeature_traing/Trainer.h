#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>

#include "FeatureExtractor.h"


using namespace std;
using namespace cv;

class Trainer
{
	public:
		Trainer(FeatureExtractor* featureExtractor);		
		~Trainer();
		void runTrainProcess();

		void loadTrainData(string posSamplesDir, string negSamplesDir);
		void resizeTrainingImages(int targetWidth, int targetHeight);
		void generateHogTrainingFeature(vector<Mat> &inputsImages, float scale);
		void fillTrainingLabel();
		void TrainLibSVM(string outputFileName);
		
		void saveHOGDescriptor(string outputFileName);
		

		void getFilesInDirectory(const string& dirName, vector<cv::Mat>& files, const vector<string>& validExtensions);
		FeatureExtractor * _featureExtractor;		
	private:		
		
		vector<Mat> _positiveTrainingImages;
		vector<Mat> _negativeTrainingImages;
		vector<Mat> _loadedImages;
		vector<Mat> _allTrainningImages;
		Mat _trainingMat;
		Mat _labelMat;
};