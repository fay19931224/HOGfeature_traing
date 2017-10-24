#include <iostream>

#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>

#include "FeatureExtractor.h"


int main()
{
	FeatureExtractor extractor;
	Mat positive = extractor.ExtractorPositiveSample();
	printf("positive extract finish!\n");
	Mat negative = extractor.ExtractorNegativeSample();
	printf("negative extract finish!\n");

	Mat trainingDataMat(positive.rows + negative.rows, positive.cols, CV_32FC1);
	//將正樣本資料放入mixData
	memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
	//將負樣本資料放入mixData
	memcpy(&trainingDataMat.data[sizeof(float)*positive.rows*positive.cols], negative.data, sizeof(float)*negative.rows*negative.cols);
	Mat dataProperty(positive.rows + negative.rows, 1, CV_32SC1, Scalar(-1.0));
	dataProperty.rowRange(0, positive.rows) = Scalar(1.0);

	std::cout << "memcpy finish!"<<std::endl;
	

	positive.release();
	negative.release();

	Ptr<ml::SVM> svm = ml::SVM::create();
	svm->setType(cv::ml::SVM::Types::C_SVC);
	svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
	svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, 100000, 1e-6));
	cv::Ptr<cv::ml::TrainData> td = cv::ml::TrainData::create(trainingDataMat, cv::ml::SampleTypes::ROW_SAMPLE, dataProperty);
	svm->train(td);
	svm->save("1024.xml");
	/*CvSVM svm;
	CvSVMParams params;
	params.svm_type = SVM::C_SVC;
	params.kernel_type = SVM::LINEAR;
	params.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 100000, 1e-6);	
	svm.train(mixData, dataProperty, Mat(), Mat(), params);//訓練數據、分類結果、參數設定
	svm.save("svmFeature1002.xml");*/
	std::cout << "svm training is done!" << std::endl;	
	system("pause");
	return 0;
}

