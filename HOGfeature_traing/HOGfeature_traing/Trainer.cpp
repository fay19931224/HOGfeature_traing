#include "Trainer.h"

Trainer::Trainer(FeatureExtractor* featureExtractor)
{
	_featureExtractor = featureExtractor;
}


Trainer::~Trainer()
{
}

//void Trainer::loadTrainData(string posSamplesDir, string negSamplesDir)
//{
//	vector<string> validExtensions;
//	validExtensions.push_back("jpg");
//	validExtensions.push_back("png");
//
//	vector<Mat> tempNegVector;
//	_utility.getFilesInDirectory(posSamplesDir, _positiveTrainingImages, validExtensions);
//	_utility.getFilesInDirectory(negSamplesDir, tempNegVector, validExtensions);
//	
//	//for (int negIndex = 0; negIndex < tempNegVector.size(); negIndex++)
//	//{
//	//	for (int heightIndex = 0; heightIndex <= tempNegVector[negIndex].rows - _programAttributes.windowHeight; heightIndex += _programAttributes.pixelsInCell)
//	//	{
//	//		for (int widthIndex = 0; widthIndex <= tempNegVector[negIndex].cols - _programAttributes.windowWidth; widthIndex += _programAttributes.pixelsInCell)
//	//		{
//	//			Mat frame_roi = tempNegVector[negIndex](Range(heightIndex, heightIndex + _programAttributes.windowHeight), cv::Range(widthIndex, widthIndex + _programAttributes.windowWidth));
//
//	//			_negativeTrainingImages.push_back(frame_roi);
//	//		}
//	//	}
//	//}
//
//
//	//如果負樣本比windowsize大，擷取部分位置
//	/*
//	|一一一一一一一一一|
//	||一一| 		   |
//	||	  |			   |
//	||一一|			   |
//	|				   |
//	|一一一一一一一一一|
//	*/
//	for (int negIndex = 0; negIndex < tempNegVector.size(); negIndex++)
//	{
//		if (tempNegVector[negIndex].rows > _programAttributes.windowHeight && tempNegVector[negIndex].cols > _programAttributes.windowWidth)
//		{
//			for (int heightIndex = 0; heightIndex <= tempNegVector[negIndex].rows - _programAttributes.windowHeight; heightIndex += (_programAttributes.windowHeight / 2))
//			{
//				for (int widthIndex = 0; widthIndex <= tempNegVector[negIndex].cols - _programAttributes.windowWidth; widthIndex += (_programAttributes.windowWidth / 2))
//				{
//					Mat frame_roi = tempNegVector[negIndex](Range(heightIndex, heightIndex + _programAttributes.windowHeight), Range(widthIndex, widthIndex + _programAttributes.windowWidth));
//
//					_negativeTrainingImages.push_back(frame_roi);
//				}
//			}
//		}
//		else
//		{
//			_negativeTrainingImages.push_back(tempNegVector[negIndex]);
//		}
//	}
//	
//
//	for (int index = 0; index < _positiveTrainingImages.size(); index++)
//	{
//		_loadedImages.push_back(_positiveTrainingImages[index]);
//	}
//	for (int index = 0; index < _negativeTrainingImages.size(); index++)
//	{
//		_loadedImages.push_back(_negativeTrainingImages[index]);
//	}
//}

//void Trainer::generateHogTrainingFeature(vector<Mat>& inputsImages, float scale)
//{
//	int X_Blocks_Num = (_programAttributes.windowWidth * scale) / _programAttributes.pixelsInCell; 
//	int Y_Blocks_Num = (_programAttributes.windowHeight * scale) / _programAttributes.pixelsInCell;
//	int blockStep = _programAttributes.cellsInBlock * _programAttributes.cellsInBlock * _programAttributes.orientations;
//
//	//(X_Blocks_Num - 1) * (Y_Blocks_Num - 1) * blockStep 的值為維度
//	_trainingMat = Mat(inputsImages.size(), (X_Blocks_Num - 1) * (Y_Blocks_Num - 1) * blockStep, CV_32FC1);
//
//	HogFeature mHog(_programAttributes);
//
//	cout << "Calculating Hog Feature..." << endl;
//
//	for (int index = 0; index < inputsImages.size(); index++)
//	{
//		int innerIndex = 0;
//
//		vector<float> descriptors;
//		mHog.calculateHogFeature(inputsImages[index], descriptors);
//
//		cout << "\rNow Current : " << index << flush;
//
//		for (int j = 0; j < descriptors.size(); j++)
//		{
//			_trainingMat.at<float>(index, innerIndex++) = descriptors.at(j);
//		}
//	}
//	cout << endl;
//}
//
//void Trainer::fillTrainingLabel()
//{
//	int num_total_files = _positiveTrainingImages.size() + _negativeTrainingImages.size();
//	float* label = new float[num_total_files];
//	for (int i = 0; i < _positiveTrainingImages.size(); i++)
//	{
//		label[i] = 1.0;
//	}
//
//	for (int i = _positiveTrainingImages.size(); i < num_total_files; i++)
//	{
//		label[i] = -1.0;
//	}
//	_labelMat = Mat(num_total_files, 1, CV_32FC1, label);
//}
//
//void Trainer::TrainLibSVM(string outputFileName)
//{
//	svm_train_class libSVM;
//
//	cout << "Training SVM..." << endl;
//	libSVM.customTrain(outputFileName, _trainingMat, _labelMat);
//
//	svm_predict_class svm;
//	svm.loadModel(outputFileName);
//	SVMSimplification(svm.model, outputFileName);
//}
//
//void Trainer::SVMSimplification(svm_model* model, string outputFileName)
//{
//	vector<float> decisionFunction(model->svLength);
//	for (int innerIndex = 0; innerIndex < decisionFunction.size(); innerIndex++)
//	{
//		decisionFunction[innerIndex] = 0;
//	}
//
//	for (int svIndex = 0; svIndex < model->l; svIndex++)
//	{
//		for (int innerIndex = 0; innerIndex < model->svLength; innerIndex++)
//		{
//			//cout << (model->SV[svIndex] + innerIndex)->value << endl;
//			decisionFunction[innerIndex] += ((model->SV[svIndex] + innerIndex)->value * model->sv_coef[0][svIndex]);
//		}
//	}
//
//	ofstream myfile(outputFileName);
//	if (myfile.is_open())
//	{
//		myfile << *(model->rho) << "\n";
//
//		for (int index = 0; index < decisionFunction.size(); index++)
//		{
//			if (index != decisionFunction.size() - 1)
//			{
//				myfile << decisionFunction[index] << ", ";
//			}
//			else
//			{
//				myfile << decisionFunction[index];
//			}
//		}
//		//myfile << "This is a line.\n";
//		//myfile << "This is another line.\n";
//		myfile.close();
//	}
//}
//
//
//void Trainer::resizeTrainingImages(int targetWidth, int targetHeight)
//{
//	_allTrainningImages.clear();
//
//	for (int index = 0; index < _loadedImages.size(); index++)
//	{
//		Mat resizedMat;
//		resize(_loadedImages[index], resizedMat, cv::Size(targetWidth, targetHeight));
//		_allTrainningImages.push_back(resizedMat);
//	}
//}

void Trainer::runTrainProcess()
{	
	_featureExtractor->run();
}

//void Trainer::saveHOGDescriptor(string outputFileName)
//{
//	ofstream outFile;
//	outFile.open(outputFileName, std::ofstream::out | std::ofstream::trunc);
//
//	Tools utility;
//
//	for (int index = 0; index < _trainingMat.rows; index++)
//	{
//		string outputString;
//		string s;
//		stringstream ss(s);
//		int label = _labelMat.at<float>(index, 0);
//		if (label > 0)
//		{
//			//outputString += "+";
//			//outputString += utility.itostr(label);
//			//outputString += " ";
//			ss << "+" << label << " ";
//		}
//		else if (label == -1)
//		{
//			//outputString += label;
//			//outputString += " ";
//			ss << label << " ";
//		}
//
//		for (int colIndex = 0; colIndex < _trainingMat.cols; colIndex++)
//		{
//			string temp;
//			int innerIndex = colIndex + 1;
//			float value = _trainingMat.at<float>(index, colIndex);
//
//			//outputString += utility.itostr(innerIndex);
//			//outputString += ":";
//			//outputString += utility.itostr((int)(value * 1000000000));
//			//outputString += " ";
//			ss << innerIndex << ":" << value << " ";
//		}
//		//outputString += "\n";
//		ss << "\n";
//		outputString += ss.str();
//		outFile << outputString;
//
//		cout << "\rNow Current : " << index << flush;
//	}
//
//	outFile.close();
//}
