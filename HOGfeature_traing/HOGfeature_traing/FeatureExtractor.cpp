#include "FeatureExtractor.h"

FeatureExtractor::FeatureExtractor()
{	
	_hogFeature = new HogFeatureParameter();

	/*dir = "ROI/";
	_winSize = Size(72, 104);*/

	
	dir = "�������I��/";
	_hogFeature->winSize = Size(48, 104);
	
	
	/*dir = "������������/";
	_hogFeature->winSize = Size(72, 88);*/

	_hogFeature->cellSize = Size(8, 8);
	_hogFeature->blockStride = Size(8, 8);
	_hogFeature->blockSize= Size(_hogFeature->cellSize.width * 2, _hogFeature->cellSize.height * 2);
	_hogFeature->nbins = 9;
	_hogFeature->derivAperture = 1;
	_hogFeature->winSigma = -1;
	_hogFeature->histogramNormType = HOGDescriptor::L2Hys;
	_hogFeature->L2HysThreshold = 0.2;
	_hogFeature->gammaCorrection = true;
	_hogFeature->nlevels = HOGDescriptor::DEFAULT_NLEVELS;
	_hogFeature->signedGradient = false;
		
	_descriptor = new HOGDescriptor(_hogFeature->winSize, _hogFeature->blockSize, _hogFeature->blockStride, _hogFeature->cellSize, _hogFeature->nbins, _hogFeature->derivAperture, _hogFeature->winSigma, _hogFeature->histogramNormType, _hogFeature->L2HysThreshold, _hogFeature->gammaCorrection, _hogFeature->nlevels, _hogFeature->signedGradient);		
}

FeatureExtractor::FeatureExtractor(string modelFileName)
{		
	_hogFeature = new HogFeatureParameter();
	ifstream modelFile(modelFileName, ifstream::in);	
	if (!modelFile.is_open())
	{
		cout << "fail to open file: " << modelFileName<< endl;
	}
	string temp;
	string token;	
	stringstream ss;
	
	if (getline(modelFile, temp)) 
	{			
		token = temp.substr(temp.find(":")+1,temp.size());		
		dir = token;
		cout << "dir:" << dir << endl;
	}
	if (getline(modelFile, temp))
	{		
		token = temp.substr(temp.find(":") + 1, temp.size());				
		int width, height;		
		ss << token.substr(0, token.find(","));
		ss >> width;
		ss.clear();
		ss << token.substr(token.find(",")+1, token.size());
		ss >> height;
		_hogFeature->winSize = Size(width, height);
		cout <<"winSize:"<< width <<","<<height<< endl;		
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		int width, height;
		ss.clear();
		ss << token.substr(0, token.find(","));
		ss >> width;
		ss.clear();
		ss << token.substr(token.find(",") + 1, token.size());
		ss >> height;
		_hogFeature->blockSize = Size(width, height);
		cout << "blockSize:" << width << "," << height << endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		int width, height;
		ss.clear();		
		ss << token.substr(0, token.find(","));
		ss >> width;
		ss.clear();
		ss << token.substr(token.find(",") + 1, token.size());
		ss >> height;
		_hogFeature->blockStride = Size(width, height);
		cout << "blockStride :" << width << "," << height << endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		int width, height;
		ss.clear();
		ss << token.substr(0, token.find(","));
		ss >> width;
		ss.clear();
		ss << token.substr(token.find(",") + 1, token.size());
		ss >> height;
		_hogFeature->cellSize = Size(width, height);
		cout << "cellSize :" << width << "," << height << endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());		
		int tempDigit;
		ss.clear();
		ss << token;
		ss >> tempDigit;
		_hogFeature->nbins = tempDigit;
		cout << "nbins :" << _hogFeature->nbins << endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		int tempDigit;
		ss.clear();
		ss << token;
		ss >> tempDigit;
		_hogFeature->derivAperture = tempDigit;
		cout << "derivAperture  :" << _hogFeature->derivAperture << endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		int tempDigit;
		ss.clear();
		ss << token;
		ss >> tempDigit;
		_hogFeature->winSigma = tempDigit;
		cout << "winSigma:" << _hogFeature->winSigma<< endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		int tempDigit;
		ss.clear();
		ss << token;
		ss >> tempDigit;
		_hogFeature->histogramNormType = tempDigit;
		cout << "histogramNormType :" << _hogFeature->histogramNormType << endl;
	}
	if (getline(modelFile, temp))
	{
		float tempDigit;
		token = temp.substr(temp.find(":") + 1, temp.size());
		ss.clear();
		ss << token;
		ss >> tempDigit;
		_hogFeature->L2HysThreshold = tempDigit;
		cout << "L2HysThreshold:" << _hogFeature->L2HysThreshold << endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		if (!strcmp(token.c_str(), "true")) 
		{
			_hogFeature->gammaCorrection = true;
			cout << "gammaCorrection :true"<< endl;
		}
		else if (!strcmp(token.c_str(), "false")) 
		{
			_hogFeature->gammaCorrection = false;
			cout << "gammaCorrection :false" << endl;
		}
		else
		{
			throw exception("gammaCorrection enter error");
		}		
	}
	if (getline(modelFile, temp))
	{
		int tempDigit;
		token = temp.substr(temp.find(":") + 1, temp.size());
		ss.clear();
		ss << token;
		ss >> tempDigit;		
		_hogFeature->nlevels = tempDigit;
		cout << "nlevels:" << _hogFeature->nlevels << endl;
	}
	if (getline(modelFile, temp))
	{
		token = temp.substr(temp.find(":") + 1, temp.size());
		if (!strcmp(token.c_str(), "true"))
		{
			_hogFeature->signedGradient = true;
			cout << "signedGradient :true" << endl;
		}
		else if (!strcmp(token.c_str(), "false"))
		{
			_hogFeature->signedGradient = false;
			cout << "signedGradient :false" << endl;
		}
		else
		{
			throw exception("signedGradient enter error");
		}
	}
	modelFile.close();
	cout << endl;
	_descriptor = new HOGDescriptor(_hogFeature->winSize, _hogFeature->blockSize, _hogFeature->blockStride, _hogFeature->cellSize, _hogFeature->nbins, _hogFeature->derivAperture, _hogFeature->winSigma, _hogFeature->histogramNormType, _hogFeature->L2HysThreshold, _hogFeature->gammaCorrection, _hogFeature->nlevels, _hogFeature->signedGradient);		
}

FeatureExtractor::~FeatureExtractor()
{
}

const HogFeatureParameter * FeatureExtractor::getHogFeature()
{
	return _hogFeature;
}

Mat FeatureExtractor::ExtractorPositiveSample()
{	                	
	string path = "pos/" + dir;
	ifstream ifs(path + "000positive.txt", ios::in);
	//vector<vector<Point>> locationList;
	vector<vector<float>> descriptorValueList;

	//fstream fp;
	//fp.open("���I��pos.txt", ios::out);//�}���ɮ�
	//if (!fp)
	//{//�p�G�}���ɮץ��ѡAfp��0�F���\�Afp���D0
	//	cout << "Fail to open file: " << endl;
	//}
	

	if (!ifs.is_open())
	{
		cout << "fail" << endl;
		system("PAUSE");
		return Mat();
	}
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt")
		{
			Mat img = imread(path + fileName, 0);					
			if (img.empty())
			{
				cout << fileName << "error" << endl;
				system("PAUSE");
				continue;
			}
			else {
				//cout << fileName << endl;
			}
			//resize(img, img, WINDOW_SIZE);
			//imshow("img", img);
			//waitKey(0);						
			vector<float> descriptorValue;
			
			_descriptor->compute(img, descriptorValue);
			descriptorValueList.push_back(descriptorValue);
			flip(img, img,1);

			_descriptor->compute(img, descriptorValue);
			descriptorValueList.push_back(descriptorValue);
			
		}
	}	
	//fp.close();//�����ɮ�	

	int row = descriptorValueList.size();
	int col = descriptorValueList[0].size();
	Mat sample(row, col, CV_32F);
	for (int i = 0; i < row; i++)
	{
		memcpy(&(sample.data[col*i * sizeof(float)]), descriptorValueList[i].data(), col * sizeof(float));
	}
	return sample;
}

Mat FeatureExtractor::ExtractorNegativeSample()
{
	string path = "neg\\"+dir;
	
	ifstream ifs(path + "\\000negative.txt", ios::in);
	//vector<vector<Point>> locationList;
	vector<vector<float>> descriptorValueList;

	//fstream fp;
	//fp.open("���I��neg.txt", ios::out);//�}���ɮ�
	//if (!fp)
	//{//�p�G�}���ɮץ��ѡAfp��0�F���\�Afp���D0
	//	cout << "Fail to open file: " << endl;
	//}


	if (!ifs.is_open())
	{
		cout << "fail" << endl;
		return Mat();
	}
	
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt")
		{			
			Mat img = imread(path + fileName, 0);
			
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				system("PAUSE");
				continue;
			}
			else {
				//cout << fileName << endl;
			}
			//resize(img, img, WINDOW_SIZE);
			//imshow("img", img);
			//waitKey(0);			
		
			vector<float> descriptorValue;
			//_descriptor->compute(img, descriptorValue, _winSize, _cellSize);
			_descriptor->compute(img, descriptorValue);
			descriptorValueList.push_back(descriptorValue);
			//cout << descriptorValue.size() << endl;
			/*for (int i = 0; i<descriptorValue.size(); i++)
			{
				fp << descriptorValue[i] << ",";
			}
			fp << endl;*/
		}
	}
	//fp.close();//�����ɮ�	
	int row = descriptorValueList.size();
	int col = descriptorValueList[0].size();
	Mat sample(row, col, CV_32F);
	for (int i = 0; i < row; i++)
	{
		memcpy(&(sample.data[col*i * sizeof(float)]), descriptorValueList[i].data(), col * sizeof(float));
	}
	return sample;
}

Mat FeatureExtractor::ShowHOGFeature(string fileName)
{
	Mat pic = imread(fileName, 0);	
	imshow("Orignal", pic);		
	
	vector<float> features;
	_descriptor->setVersion(false);
	_descriptor->compute(pic, features);
	fstream file2;
	file2.open("before.txt", ios::out);
	for (int i = 0; i<features.size(); i++)
	{
		file2 << features[i] << endl;
	}
	file2.close();
	Mat visualization = get_hogdescriptor_visu(imread(fileName, 1), features, _hogFeature->winSize);
	cv::imshow("Visualization before", visualization);
	cv::moveWindow("Visualization before", 500, 500);

	vector<float> features_new;
	_descriptor = new HOGDescriptor(_hogFeature->winSize, _hogFeature->blockSize, 
		_hogFeature->blockStride, _hogFeature->cellSize, _hogFeature->nbins,
		_hogFeature->derivAperture, _hogFeature->winSigma, 0,
		_hogFeature->L2HysThreshold, _hogFeature->gammaCorrection, _hogFeature->nlevels,
		_hogFeature->signedGradient);
	_descriptor->setVersion(true);
	_descriptor->compute(pic, features_new);
	fstream file;
	file.open("after.txt", ios::out);
	for (int i = 0; i<features_new.size(); i++)
	{
		file << features_new[i] << endl;
	}
	file.close();	
	Mat visualization_new = get_hogdescriptor_visu(imread(fileName, 1), features_new, _hogFeature->winSize);
	
	cv::imshow("Visualization after", visualization_new);
	cv::moveWindow("Visualization after", 750, 500);
	cv::waitKey();
	return pic;
}

Mat FeatureExtractor::get_hogdescriptor_visu(const Mat& color_origImg, vector<float>& descriptorValues, const Size & size)
{
	const int DIMX = size.width;
	const int DIMY = size.height;
	float zoomFac = 3;
	Mat visu;
	resize(color_origImg, visu, Size((int)(color_origImg.cols*zoomFac), (int)(color_origImg.rows*zoomFac)));

	Size cvcellSize = _hogFeature->cellSize;
	int cellSize = cvcellSize.height;
	int gradientBinSize = _hogFeature->nbins;
	float radRangeForOneBin = (float)(CV_PI / (float)gradientBinSize); // dividing 180 into 9 bins, how large (in rad) is one bin?

																	   // prepare data structure: 9 orientation / gradient strenghts for each cell
	int cells_in_x_dir = DIMX / cellSize;
	int cells_in_y_dir = DIMY / cellSize;
	float*** gradientStrengths = new float**[cells_in_y_dir];
	int** cellUpdateCounter = new int*[cells_in_y_dir];
	for (int y = 0; y<cells_in_y_dir; y++)
	{
		gradientStrengths[y] = new float*[cells_in_x_dir];
		cellUpdateCounter[y] = new int[cells_in_x_dir];
		for (int x = 0; x<cells_in_x_dir; x++)
		{
			gradientStrengths[y][x] = new float[gradientBinSize];
			cellUpdateCounter[y][x] = 0;

			for (int bin = 0; bin<gradientBinSize; bin++)
				gradientStrengths[y][x][bin] = 0.0;
		}
	}

	// nr of blocks = nr of cells - 1
	// since there is a new block on each cell (overlapping blocks!) but the last one
	int blocks_in_x_dir = cells_in_x_dir - 1;
	int blocks_in_y_dir = cells_in_y_dir - 1;

	// compute gradient strengths per cell
	int descriptorDataIdx = 0;
	int cellx = 0;
	int celly = 0;

	for (int blockx = 0; blockx<blocks_in_x_dir; blockx++)
	{
		for (int blocky = 0; blocky<blocks_in_y_dir; blocky++)
		{
			// 4 cells per block ...
			for (int cellNr = 0; cellNr<4; cellNr++)
			{
				// compute corresponding cell nr
				cellx = blockx;
				celly = blocky;
				if (cellNr == 1) celly++;
				if (cellNr == 2) cellx++;
				if (cellNr == 3)
				{
					cellx++;
					celly++;
				}

				for (int bin = 0; bin<gradientBinSize; bin++)
				{
					float gradientStrength = descriptorValues[descriptorDataIdx];
					descriptorDataIdx++;

					gradientStrengths[celly][cellx][bin] += gradientStrength;

				} // for (all bins)


				  // note: overlapping blocks lead to multiple updates of this sum!
				  // we therefore keep track how often a cell was updated,
				  // to compute average gradient strengths
				cellUpdateCounter[celly][cellx]++;

			} // for (all cells)


		} // for (all block x pos)
	} // for (all block y pos)


	  // compute average gradient strengths
	for (celly = 0; celly<cells_in_y_dir; celly++)
	{
		for (cellx = 0; cellx<cells_in_x_dir; cellx++)
		{

			float NrUpdatesForThisCell = (float)cellUpdateCounter[celly][cellx];

			// compute average gradient strenghts for each gradient bin direction
			for (int bin = 0; bin<gradientBinSize; bin++)
			{
				gradientStrengths[celly][cellx][bin] /= NrUpdatesForThisCell;
			}
		}
	}

	// draw cells
	for (celly = 0; celly<cells_in_y_dir; celly++)
	{
		for (cellx = 0; cellx<cells_in_x_dir; cellx++)
		{
			int drawX = cellx * cellSize;
			int drawY = celly * cellSize;

			int mx = drawX + cellSize / 2;
			int my = drawY + cellSize / 2;

			rectangle(visu, Point((int)(drawX*zoomFac), (int)(drawY*zoomFac)), Point((int)((drawX + cellSize)*zoomFac), (int)((drawY + cellSize)*zoomFac)), Scalar(0, 0, 255), 1);

			// draw in each cell all 9 gradient strengths
			for (int bin = 0; bin<gradientBinSize; bin++)
			{
				float currentGradStrength = gradientStrengths[celly][cellx][bin];

				// no line to draw?
				if (currentGradStrength == 0)
					continue;

				float currRad = bin * radRangeForOneBin + radRangeForOneBin / 2;

				float dirVecX = cos(currRad);
				float dirVecY = sin(currRad);
				float maxVecLen = (float)(cellSize / 2.f);
				float scale = 2.5; // just a visualization scale, to see the lines better

								   // compute line coordinates
				float x1 = mx - dirVecX * currentGradStrength * maxVecLen * scale;
				float y1 = my - dirVecY * currentGradStrength * maxVecLen * scale;
				float x2 = mx + dirVecX * currentGradStrength * maxVecLen * scale;
				float y2 = my + dirVecY * currentGradStrength * maxVecLen * scale;

				// draw gradient visualization
				line(visu, Point((int)(x1*zoomFac), (int)(y1*zoomFac)), Point((int)(x2*zoomFac), (int)(y2*zoomFac)), Scalar(255, 0, 0), 1);

			} // for (all bins)

		} // for (cellx)
	} // for (celly)


	  // don't forget to free memory allocated by helper data structures!
	for (int y = 0; y<cells_in_y_dir; y++)
	{
		for (int x = 0; x<cells_in_x_dir; x++)
		{
			delete[] gradientStrengths[y][x];
		}
		delete[] gradientStrengths[y];
		delete[] cellUpdateCounter[y];
	}
	delete[] gradientStrengths;
	delete[] cellUpdateCounter;

	return visu;

} // get_hogdescriptor_visu

void FeatureExtractor::loadTrainingData(string posPath,string negPath)
{		
	//load pos data
	ifstream ifs(posPath + "000positive.txt", ios::in);
	if (!ifs.is_open())
	{
		throw exception("open POS data fail");		
	}
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt")
		{			
			Mat img = imread(posPath + fileName, IMREAD_GRAYSCALE);
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				throw exception("load POS data fail");				
			}					
			PosTrainingData.push_back(img);			
			flip(img, img, 1);			
			PosTrainingData.push_back(img);
		}
	}
	ifs.close();
	cout << "load POS data success" << endl;

	//load neg data
	ifs.open(negPath + "000negative.txt", ios::in);
	if (!ifs.is_open())
	{
		throw exception("open NEG data fail");
	}
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt")
		{
			Mat img = imread(negPath + fileName, IMREAD_GRAYSCALE);
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				throw exception("load NEG data fail");
			}
			NegTrainingData.push_back(img);			
		}
	}
	ifs.close();
	cout << "load NEG data success" << endl;
}

void FeatureExtractor::extractSample()
{	
	Mat positive = ExtractorPositiveSample();
	std::cout << "positive extract finish!" << std::endl;
	Mat negative = ExtractorNegativeSample();
	std::cout << "negative extract finish!" << std::endl;

	Mat trainingDataMat(positive.rows + negative.rows, positive.cols, CV_32FC1);
	//�N���˥���Ʃ�JmixData
	memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
	//�N�t�˥���Ʃ�JmixData
	memcpy(&trainingDataMat.data[sizeof(float)*positive.rows*positive.cols], negative.data, sizeof(float)*negative.rows*negative.cols);
	Mat dataProperty(positive.rows + negative.rows, 1, CV_32SC1, Scalar(-1.0));
	dataProperty.rowRange(0, positive.rows) = Scalar(1.0);

	std::cout << "memcpy finish!" << std::endl;

	positive.release();
	negative.release();


	/*string  filename = "���I��0123C_SVC_LINEAR";
	Ptr<ml::SVM> svm = ml::SVM::create();
	svm->setType(cv::ml::SVM::Types::C_SVC);
	svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
	svm->setTermCriteria(cv::TermCriteria(CV_TERMCRIT_ITER, 100000, 1e-6));
	svm->train(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);	*/


	string  filename = "����0308_auto_LINEAR_gamm_new2";
	svm = ml::SVM::create();
	svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
	Ptr<ml::TrainData> data = ml::TrainData::create(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);
	svm->trainAuto(data);


	svm->save(filename + ".xml");

	writeTrainingParameter(filename);

	std::cout << "svm training is done!" << std::endl;
}

void FeatureExtractor::run()
{	
	string posPath = "pos/" + dir;
	string negPath = "neg/" + dir;
	loadTrainingData(posPath, negPath);
	double t = (double)cv::getTickCount();
	runTrainProcessOpenCV();
	//runTrainProcess();
	t = (double)cv::getTickCount() - t;	
	cout << "Spending time:" <<t<< endl;
}

//void FeatureExtractor::runTrainProcess()
//{	
//	vector<Mat> tempNegVector = NegTrainingData;
//	NegTrainingData.clear();
//	for (int negIndex = 0; negIndex < tempNegVector.size(); negIndex++)
//	{
//		if (tempNegVector[negIndex].rows > _hogFeature->winSize.height && tempNegVector[negIndex].cols > _hogFeature->winSize.width)
//		{
//			for (int heightIndex = 0; heightIndex <= tempNegVector[negIndex].rows - _hogFeature->winSize.height; heightIndex += (_hogFeature->winSize.height / 2))
//			{
//				for (int widthIndex = 0; widthIndex <= tempNegVector[negIndex].cols - _hogFeature->winSize.width; widthIndex += (_hogFeature->winSize.width / 2))
//				{
//					Mat frame_roi = tempNegVector[negIndex](Range(heightIndex, heightIndex + _hogFeature->winSize.height), Range(widthIndex, widthIndex + _hogFeature->winSize.width));
//					NegTrainingData.push_back(frame_roi);
//				}
//			}
//		}
//		else
//		{
//			NegTrainingData.push_back(tempNegVector[negIndex]);
//		}
//	}
//
//	for (int index = 0; index < PosTrainingData.size(); index++)
//	{
//		_loadedImages.push_back(PosTrainingData[index]);
//	}
//	for (int index = 0; index < NegTrainingData.size(); index++)
//	{
//		_loadedImages.push_back(NegTrainingData[index]);
//	}
//
//	Mat labelMat=fillTrainingLabel();
//
//	float trainScaleStart = 1;
//	float trainScaleTarget = 0.7;
//	float trainScaleStep = 0.1;
//
//	
//
//	for (float trainScale = trainScaleStart; trainScale>= trainScaleTarget; trainScale -= trainScaleStep)
//	{
//		cout << "Now Training " << trainScale << endl;
//		vector<Mat> scaledloadedImages;
//		int scaledheight = _loadedImages[0].cols*trainScale;
//		int sacledwidth = _loadedImages[0].rows*trainScale;
//		for (int index = 0; index < _loadedImages.size(); index++) 
//		{
//			Mat resizedMat;
//			resize(_loadedImages[index], resizedMat, cv::Size(sacledwidth, scaledheight));
//			scaledloadedImages.push_back(resizedMat);
//		}				
//		generateHogTrainingFeature(scaledloadedImages, trainScale);
//		
//		stringstream ss;
//		ss << trainScale;		
//
//		string  filename = "new"+ss.str() +"MOTOSIDE";
//		svm = ml::SVM::create();
//		svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
//
//
//
//		Ptr<ml::TrainData> data = ml::TrainData::create(_trainingMat, ml::SampleTypes::ROW_SAMPLE, labelMat);
//		svm->trainAuto(data);
//
//		svm->save(filename + ".xml");
//		writeTrainingParameter(filename);		
//	}
//	std::cout << "svm training is done!" << std::endl;
//}

void FeatureExtractor::runTrainProcessOpenCV()
{	
	
	cout << "run OpenCV version" << endl;
	
	if (cv::useOptimized()) 
	{
		cout << "OpenCV optimized CPU open" << endl;
	}	
	if (cv::checkHardwareSupport(CV_SSE2))
	{
		cout << "OpenCV open CV_SSE2" << endl;
	}
	if (cv::checkHardwareSupport(CV_NEON))
	{
		cout << "OpenCV open CV_NEON" << endl;
	}
	if (cuda::getCudaEnabledDeviceCount()) 
	{
		cout << "OpenCV with cuda" << endl;
	}
	//_descriptor = new HOGDescriptor(_hogFeature->winSize, _hogFeature->blockSize, _hogFeature->blockStride, _hogFeature->cellSize, _hogFeature->nbins, _hogFeature->derivAperture, _hogFeature->winSigma, _hogFeature->histogramNormType, _hogFeature->L2HysThreshold, _hogFeature->gammaCorrection, _hogFeature->nlevels, _hogFeature->signedGradient);
	_descriptor->setVersion(true);	 
	std::cout << "Start to computing orientation and magnitude" << std::endl;
	for (int i = 0;i<PosTrainingData.size(); i++) 
	{	
		
		vector<float> descriptorValue;
		
		double t = (double)cv::getTickCount();

		/*UMat U = PosTrainingData[i].getUMat(1);
		InputArray temp = U;
		if (temp.dims() <= 2 && temp.type() == CV_8UC1&&temp.isUMat()) 
		{
			cout << "OCL" << endl;
		}
		else
		{
			cout << "NO OCL" << endl;
		}			
		_descriptor->compute(temp, descriptorValue);*/

		_descriptor->compute(PosTrainingData[i], descriptorValue);				
		t = (double)cv::getTickCount() - t;
		//cout << "Spend Time:"<<t/getTickFrequency() <<"s"<< endl;		
		descriptorValueList.push_back(descriptorValue);		
	}
	int row = descriptorValueList.size();
	int col = descriptorValueList[0].size();
	
	Mat positive(row, col, CV_32F);
	/*------------
	|descriptorValue
	|descriptorValue
	|descriptorValue
	|descriptorValue
	----------------
	*/
	for (int i = 0; i < row; i++)
	{
		memcpy(&(positive.data[col*i * sizeof(float)]), descriptorValueList[i].data(), col * sizeof(float));
	}
	
	descriptorValueList.clear();

	for (int i = 0; i<NegTrainingData.size(); i++)
	{
		vector<float> descriptorValue;		
		_descriptor->compute(NegTrainingData[i], descriptorValue);			
		descriptorValueList.push_back(descriptorValue);
	}
	row = descriptorValueList.size();
	col = descriptorValueList[0].size();
	Mat negative(row, col, CV_32F);
	for (int i = 0; i < row; i++)
	{
		memcpy(&(negative.data[col*i * sizeof(float)]), descriptorValueList[i].data(), col * sizeof(float));
	}

	Mat trainingDataMat(positive.rows + negative.rows, positive.cols, CV_32FC1);
	//�N���˥���Ʃ�JmixData
	memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
	//�N�t�˥���Ʃ�JmixData
	memcpy(&trainingDataMat.data[sizeof(float)*positive.rows*positive.cols], negative.data, sizeof(float)*negative.rows*negative.cols);
	Mat dataProperty(positive.rows + negative.rows, 1, CV_32SC1, Scalar(-1.0));
	dataProperty.rowRange(0, positive.rows) = Scalar(1.0);

	std::cout << "Start to training" << std::endl;

	positive.release();
	negative.release();

	string filename = "����05251_auto";
	svm = ml::SVM::create();
	svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
	Ptr<ml::TrainData> data = ml::TrainData::create(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);
	svm->trainAuto(data);
	
	svm->save(filename + ".xml");
	writeTrainingParameter(filename);
	//	_descriptor->save();	
	std::cout << "svm training is done!" << std::endl;	
}

void FeatureExtractor::generateHogTrainingFeature(vector<Mat>& inputsImages, float scale)
{
	int d = _hogFeature->winSize.width;
	int X_Blocks_Num = (_hogFeature->winSize.width* scale) / _hogFeature->cellSize.width;
	int Y_Blocks_Num = (_hogFeature->winSize.height* scale) / _hogFeature->cellSize.height;
	int cellsInBlock = 2;
	int blockStep = cellsInBlock* cellsInBlock* _hogFeature->nbins;

	//(X_Blocks_Num - 1) * (Y_Blocks_Num - 1) * blockStep ���Ȭ�����
	//_trainingMat = Mat(inputsImages.size(), (X_Blocks_Num - 1) * (Y_Blocks_Num - 1) * blockStep, CV_32FC1);	
	_trainingMat = Mat(inputsImages.size(), (X_Blocks_Num - 1) * (Y_Blocks_Num - 1) * blockStep, CV_32FC1);
	cout << "Calculating Hog Feature..." << endl;
	HogFeature mHog(_hogFeature);
	for (int index = 0; index < inputsImages.size(); index++)
	{
		int innerIndex = 0;

		vector<float> descriptors;
		mHog.calculateHogFeature(inputsImages[index], descriptors);
		
		descriptorValueList.push_back(descriptors);
		cout << "\rNow Current : " << index << flush;

		for (int j = 0; j < descriptors.size(); j++)
		{
			_trainingMat.at<float>(index, innerIndex++) = descriptors.at(j);
		}
	}
	cout << endl;	
}

Mat FeatureExtractor::fillTrainingLabel()
{
	int num_total_files = PosTrainingData.size() + NegTrainingData.size();
	float* label = new float[num_total_files];
	for (int i = 0; i < PosTrainingData.size(); i++)
	{
		label[i] = 1.0;
	}

	for (int i = NegTrainingData.size(); i < num_total_files; i++)
	{
		label[i] = -1.0;
	}
	return Mat(num_total_files, 1, CV_32FC1, label);
}

void FeatureExtractor::writeTrainingParameter(string fileName)
{
	fstream fp;
	fp.open(fileName + ".txt", ios::out);//�}���ɮ�
	if (!fp)
	{//�p�G�}���ɮץ��ѡAfp��0�F���\�Afp���D0
		cout << "Fail to open file: " << fileName << endl;
	}
	int space = 25;
	cout << "SVM Para: " << setw(space) << endl;
	cout << "Type: " << svm->getType() << setw(space) << endl;
	cout << "Kernel type: " << svm->getKernelType() << setw(space) << endl;
	cout << "C: " << svm->getC() << setw(space) << endl;
	cout << "Degree: " << svm->getDegree() << setw(space) << endl;
	cout << "Nu: " << svm->getNu() << setw(space) << endl;
	cout << "Gamm: " << svm->getGamma() << setw(space) << endl;
	cout << "TermCriteria_epsilon: " << svm->getTermCriteria().epsilon << setw(space) << endl;
	cout << "TermCriteria_maxCount: " << svm->getTermCriteria().maxCount << setw(space) << endl;
	cout << "TermCriteria_type: " << svm->getTermCriteria().type << setw(space) << endl;
	cout << "----------------------------------------" << endl;

	fp << "SVM Para: " << setw(space) << endl;
	fp << "Type: " << svm->getType() << setw(space) << endl;
	fp << "Kerne type: " << svm->getKernelType() << setw(space) << endl;
	fp << "C: " << svm->getC() << setw(space) << endl;
	fp << "Degree: " << svm->getDegree() << setw(space) << endl;
	fp << "Nu: " << svm->getNu() << setw(space) << endl;
	fp << "Gamm: " << svm->getGamma() << setw(space) << endl;
	fp << "TermCriteria_epsilon: " << svm->getTermCriteria().epsilon << setw(space) << endl;
	fp << "TermCriteria_maxCount: " << svm->getTermCriteria().maxCount << setw(space) << endl;
	fp << "TermCriteria_type: " << svm->getTermCriteria().type << setw(space) << endl;
	fp << "----------------------------------------" << endl;
	fp << "HOG Para: " << setw(space) << endl;
	fp << "winSize:" << _hogFeature->winSize.width << "," << _hogFeature->winSize.height << setw(space) << endl;
	fp << "blockSize:" << _hogFeature->blockSize.width << "," << _hogFeature->blockSize.height << setw(space) << endl;
	fp << "blockStride:" << _hogFeature->blockStride.width << "," << _hogFeature->blockStride.height << setw(space) << endl;
	fp << "cellSize:" << _hogFeature->cellSize.width << ", " << _hogFeature->cellSize.height << setw(space) << endl;
	fp << "nbins:" << _hogFeature->nbins << setw(space) << endl;
	fp << "derivAperture:" << _hogFeature->derivAperture << setw(space) << endl;
	fp << "winSigma:" << _hogFeature->winSigma << setw(space) << endl;
	fp << "histogramNormType:" << _hogFeature->histogramNormType << setw(space) << endl;
	fp << "L2HysThreshold:" << _hogFeature->L2HysThreshold << setw(space) << endl;
	fp << "gammaCorrection:" << _hogFeature->gammaCorrection << setw(space) << endl;
	fp << "nlevels:" << _hogFeature->nlevels << setw(space) << endl;
	fp << "signedGradient :" << _hogFeature->signedGradient << setw(space) << endl;
	fp << "HOGDescriptorSize:" << _descriptor->getDescriptorSize() << setw(space) << endl;
	fp << "���˥���:" << PosTrainingData.size() << setw(space) << endl;
	fp << "�t�˥���:" << NegTrainingData.size() << setw(space) << endl;

	fp.close();//�����ɮ�	
}

void FeatureExtractor::runTrainProcess()
{	
	cout << "run defineBySelf version" << endl;
	vector<Mat> tempNegVector = NegTrainingData;
	NegTrainingData.clear();
	for (int negIndex = 0; negIndex < tempNegVector.size(); negIndex++)
	{
		if (tempNegVector[negIndex].rows > _hogFeature->winSize.height && tempNegVector[negIndex].cols > _hogFeature->winSize.width)
		{
			for (int heightIndex = 0; heightIndex <= tempNegVector[negIndex].rows - _hogFeature->winSize.height; heightIndex += (_hogFeature->winSize.height / 2))
			{
				for (int widthIndex = 0; widthIndex <= tempNegVector[negIndex].cols - _hogFeature->winSize.width; widthIndex += (_hogFeature->winSize.width / 2))
				{
					Mat frame_roi = tempNegVector[negIndex](Range(heightIndex, heightIndex + _hogFeature->winSize.height), Range(widthIndex, widthIndex + _hogFeature->winSize.width));
					NegTrainingData.push_back(frame_roi);
				}
			}
		}
		else
		{
			NegTrainingData.push_back(tempNegVector[negIndex]);
		}
	}


	float trainScaleStart = 1;
	//float trainScaleTarget = 0.7;
	float trainScaleTarget = 1;
	float trainScaleStep = 0.1;



	for (float trainScale = trainScaleStart; trainScale >= trainScaleTarget; trainScale -= trainScaleStep)
	{
		cout << "Now Training " << trainScale << endl;
		vector<Mat> scaledloadedImages;
		int scaledheight = PosTrainingData[0].cols*trainScale;
		int sacledwidth = PosTrainingData[0].rows*trainScale;
		for (int index = 0; index < PosTrainingData.size(); index++)
		{
			Mat resizedMat;
			resize(PosTrainingData[index], resizedMat, cv::Size(sacledwidth, scaledheight));
			scaledloadedImages.push_back(resizedMat);
		}
		for (int index = 0; index < NegTrainingData.size(); index++)
		{
			Mat resizedMat;
			resize(NegTrainingData[index], resizedMat, cv::Size(sacledwidth, scaledheight));
			scaledloadedImages.push_back(resizedMat);
		}
		HogFeature mHog(_hogFeature);

		int row, col;
		for (int i = 0; i<PosTrainingData.size(); i++)
		{			
			/*vector<float> descriptorValue;
			clock_t t;
			t = clock();
			mHog.calculateHogFeature(PosTrainingData[i], descriptorValue);
			t = clock() - t;
			printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);
			descriptorValueList.push_back(descriptorValue);*/

			
			vector<float> descriptorValue;
			//double t = (double)cv::getTickCount();
			mHog.calculateHogFeature(PosTrainingData[i], descriptorValue);						
			//t = (double)cv::getTickCount() - t;
			//cout << "Spend Time:" << t << endl;
			descriptorValueList.push_back(descriptorValue);
		}
		row = descriptorValueList.size();
		col = descriptorValueList[0].size();

		Mat positive(row, col, CV_32F);
		for (int i = 0; i < row; i++)
		{
			memcpy(&(positive.data[col*i * sizeof(float)]), descriptorValueList[i].data(), col * sizeof(float));
		}

		descriptorValueList.clear();

		for (int i = 0; i<NegTrainingData.size(); i++)
		{
			vector<float> descriptorValue;			
			mHog.calculateHogFeature(NegTrainingData[i], descriptorValue);
			descriptorValueList.push_back(descriptorValue);
		}
		row = descriptorValueList.size();
		col = descriptorValueList[0].size();
		Mat negative(row, col, CV_32F);
		for (int i = 0; i < row; i++)
		{
			memcpy(&(negative.data[col*i * sizeof(float)]), descriptorValueList[i].data(), col * sizeof(float));
		}

		Mat trainingDataMat(positive.rows + negative.rows, positive.cols, CV_32FC1);
		//�N���˥���Ʃ�JmixData
		memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
		//�N�t�˥���Ʃ�JmixData
		memcpy(&trainingDataMat.data[sizeof(float)*positive.rows*positive.cols], negative.data, sizeof(float)*negative.rows*negative.cols);
		Mat dataProperty(positive.rows + negative.rows, 1, CV_32SC1, Scalar(-1.0));
		dataProperty.rowRange(0, positive.rows) = Scalar(1.0);

		stringstream ss;
		ss << trainScale;

		string  filename = "new" + ss.str() + "MOTOSIDE";
		svm = ml::SVM::create();
		svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
		Ptr<ml::TrainData> data = ml::TrainData::create(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);
		cout << "Start to training" << endl;
		svm->trainAuto(data);

		svm->save(filename + ".xml");
		writeTrainingParameter(filename);
	}
	std::cout << "svm training is done!" << std::endl;
}