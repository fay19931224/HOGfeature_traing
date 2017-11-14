#include <iostream>

#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>

#include "FeatureExtractor.h"


int main()
{
	FeatureExtractor* extractor=new FeatureExtractor();
	Mat positive = extractor->ExtractorPositiveSample();	
	std::cout << "positive extract finish!" << std::endl;
	Mat negative = extractor->ExtractorNegativeSample();
	std::cout << "negative extract finish!" << std::endl;
	

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
	//svm->setKernel(cv::ml::SVM::KernelTypes::RBF);
	svm->setTermCriteria(cv::TermCriteria(CV_TERMCRIT_ITER, 100000, 1e-6));
	svm->train(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);
	svm->save("側面全身C_SVC_LINEAR.xml");	
	std::cout << "svm training is done!" << std::endl;	
	system("pause");
	return 0;
}
/*int main()
{
	FeatureExtractor extractor;
	Mat positive = extractor.ExtractorPositiveSample();
	std::cout << "positive extract finish!" << std::endl;

	Mat trainingDataMat(positive.rows , positive.cols, CV_32FC1);
	//將正樣本資料放入mixData
	memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);	
	Mat dataProperty(positive.rows, 1, CV_32SC1, Scalar(-1.0));
	dataProperty.rowRange(0, positive.rows) = Scalar(1.0);

	std::cout << "memcpy finish!" << std::endl;

	positive.release();	

	Ptr<ml::SVM> svm = ml::SVM::create();
	svm->setType(cv::ml::SVM::Types::ONE_CLASS);
	svm->setKernel(cv::ml::SVM::KernelTypes::RBF);
	
	svm->setNu(0.1);
	
	svm->setTermCriteria(cvTermCriteria(CV_TERMCRIT_ITER,10000, FLT_EPSILON));

	svm->train(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);
	svm->save("側面全身ONE_CLASS_RBF.xml");
	std::cout << "svm training is done!" << std::endl;
	system("pause");
	return 0;
}*/
