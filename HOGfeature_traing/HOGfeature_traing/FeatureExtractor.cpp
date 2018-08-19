#include "FeatureExtractor.h"

FeatureExtractor::FeatureExtractor()
{	
	_hogFeature = new HogFeatureParameter();

	/*dir = "ROI/";
	_winSize = Size(72, 104);*/

	
	/*dir = "機車正背面/";
	_hogFeature->winSize = Size(48, 104);*/
	
	
	/*dir = "機車側面全身/";
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
	string path = "pos/" + dir+"/";
	ifstream ifs(path + "000positive.txt", ios::in);
	//vector<vector<Point>> locationList;
	vector<vector<float>> descriptorValueList;

	//fstream fp;
	//fp.open("正背面pos.txt", ios::out);//開啟檔案
	//if (!fp)
	//{//如果開啟檔案失敗，fp為0；成功，fp為非0
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
	//fp.close();//關閉檔案	

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
	string path = "neg\\"+dir+"\\";
	
	ifstream ifs(path + "\\000negative.txt", ios::in);
	//vector<vector<Point>> locationList;
	vector<vector<float>> descriptorValueList;

	//fstream fp;
	//fp.open("正背面neg.txt", ios::out);//開啟檔案
	//if (!fp)
	//{//如果開啟檔案失敗，fp為0；成功，fp為非0
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
	
			vector<float> descriptorValue;		
			_descriptor->compute(img, descriptorValue);
			descriptorValueList.push_back(descriptorValue);
			
		}
	}
	//fp.close();//關閉檔案	
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

	/*vector<float> features_new;
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
	cv::moveWindow("Visualization after", 750, 500);*/
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

void FeatureExtractor::loadTrainingData()
{		
	//load pos data
	string posPath = "pos/" + dir + "/";
	ifstream ifs(posPath+"000positive.txt", ios::in);
	//ifstream ifs("pos/moto/000positive.txt", ios::in);
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
			Mat img = imread(posPath +fileName, IMREAD_GRAYSCALE);
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				throw exception("load POS data fail");				
			}					
			//cout << fileName << endl;			
			PosTrainingData.push_back(img);			
			flip(img, img, 1);			
			PosTrainingData.push_back(img);
			waitKey(1);
		}
	}
	ifs.close();
	cout << "load POS data success" << endl;

	//load neg data
	string negPath = "neg/" + dir + "/";
	ifs.open(negPath + "000negative.txt", ios::in);
	if (!ifs.is_open())
	{
		throw exception("open NEG data fail");
	}
	int i = 0;
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		i++;
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt"&&i%10==0)
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
	//將正樣本資料放入mixData
	memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
	//將負樣本資料放入mixData
	memcpy(&trainingDataMat.data[sizeof(float)*positive.rows*positive.cols], negative.data, sizeof(float)*negative.rows*negative.cols);
	Mat dataProperty(positive.rows + negative.rows, 1, CV_32SC1, Scalar(-1.0));
	dataProperty.rowRange(0, positive.rows) = Scalar(1.0);

	std::cout << "memcpy finish!" << std::endl;

	positive.release();
	negative.release();


	/*string  filename = "正背面0123C_SVC_LINEAR";
	Ptr<ml::SVM> svm = ml::SVM::create();
	svm->setType(cv::ml::SVM::Types::C_SVC);
	svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
	svm->setTermCriteria(cv::TermCriteria(CV_TERMCRIT_ITER, 100000, 1e-6));
	svm->train(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);	*/


	string  filename = "側面0308_auto_LINEAR_gamm_new2";
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
	loadTrainingData();
	
	unsigned long start = clock();
	runTrainProcessOpenCV();
	unsigned long end = clock();
	
	printf("total Spending time=%1.3f seconds\n", (end - start) / 1000.0);
	
	//runTrainProcess();	
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
	//_descriptor->setVersion(true);	 
	_descriptor->setVersion(false);
	

	std::cout << "Start to computing orientation and magnitude" << std::endl;
	
	for (int i = 0;i<PosTrainingData.size(); i++) 
	{					
		//double t = (double)cv::getTickCount();
		vector<float> descriptorValue;
		_descriptor->compute(PosTrainingData[i], descriptorValue);								
		//t = (double)cv::getTickCount() - t;
		//cout << "Spend Time:"<<t/getTickFrequency() <<"s"<< endl;	
		descriptorValueList.push_back(descriptorValue);		
	}
	int row = descriptorValueList.size();//樣本數
	int col = descriptorValueList[0].size();//單一樣本維度
	
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
		if (dir.compare("行人")==0) 
		{
			for (int height = 0; height <NegTrainingData[i].rows - _hogFeature->winSize.height; height += _hogFeature->winSize.height / 2)
			{
				for (int width = 0; width<NegTrainingData[i].cols - _hogFeature->winSize.width; width += _hogFeature->winSize.width / 2)
				{
					vector<float> descriptorValue;
					Rect roi = Rect(width, height, _hogFeature->winSize.width, _hogFeature->winSize.height);
					_descriptor->compute(NegTrainingData[i](roi), descriptorValue);
					descriptorValueList.push_back(descriptorValue);
				}
			}
		}
		else 
		{
			vector<float> descriptorValue;
			_descriptor->compute(NegTrainingData[i], descriptorValue);
			descriptorValueList.push_back(descriptorValue);
		}				
	}

	row = descriptorValueList.size();
	col = descriptorValueList[0].size();
	Mat negative(row, col, CV_32F);
	for (int i = 0; i < row; i++)
	{
		memcpy(&(negative.data[col*i * sizeof(float)]), descriptorValueList[i].data(), col * sizeof(float));
	}

	Mat trainingDataMat(positive.rows + negative.rows, positive.cols, CV_32FC1);
	//將正樣本資料放入mixData
	memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
	//將負樣本資料放入mixData
	memcpy(&trainingDataMat.data[sizeof(float)*positive.rows*positive.cols], negative.data, sizeof(float)*negative.rows*negative.cols);
	Mat dataProperty(positive.rows + negative.rows, 1, CV_32SC1, Scalar(-1.0));
	dataProperty.rowRange(0, positive.rows) = Scalar(1.0);

	std::cout << "Start to training" << std::endl;

	positive.release();
	negative.release();

	string filename = "helmet0818_auto";
	svm = ml::SVM::create();
	svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
	Ptr<ml::TrainData> data = ml::TrainData::create(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);
	svm->trainAuto(data);
	
	svm->save(filename + ".xml");
	writeTrainingParameter(filename);
		//_descriptor->save(filename);
	std::cout << "svm training is done!" << std::endl;	
}

void FeatureExtractor::generateHogTrainingFeature(vector<Mat>& inputsImages, float scale)
{
	int d = _hogFeature->winSize.width;
	int X_Blocks_Num = (_hogFeature->winSize.width* scale) / _hogFeature->cellSize.width;
	int Y_Blocks_Num = (_hogFeature->winSize.height* scale) / _hogFeature->cellSize.height;
	int cellsInBlock = 2;
	int blockStep = cellsInBlock* cellsInBlock* _hogFeature->nbins;

	//(X_Blocks_Num - 1) * (Y_Blocks_Num - 1) * blockStep 的值為維度
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

void FeatureExtractor::flipImage()
{
	//load pos data
	string posPath = "pos/" + dir + "/";
	string destination = "pos/行人加入水平翻轉/";
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
			Mat img = imread(posPath + fileName);
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				throw exception("load POS data fail");
			}
			//cout << fileName << endl;	
			imwrite(destination+"+" +fileName, img);
			flip(img, img, 1);
			imwrite(destination +"-" + fileName, img);
		}
	}
	ifs.close();
}

void FeatureExtractor::writeTrainingParameter(string fileName)
{
	fstream fp;
	fp.open(fileName + ".txt", ios::out);//開啟檔案
	if (!fp)
	{//如果開啟檔案失敗，fp為0；成功，fp為非0
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
	fp << "正樣本數:" << PosTrainingData.size() << setw(space) << endl;
	fp << "負樣本數:" << NegTrainingData.size() << setw(space) << endl;

	fp.close();//關閉檔案	
}

void FeatureExtractor::KFoldCrossValidation(int k)
{
	vector<vector<vector<float>>> allPosData;	
	vector<vector<Mat>> allPosDataMat;
	vector<vector<float>> NegData;
	for (int i = 0; i < k; i++)
	{
		vector<vector<float>> data;
		allPosData.push_back(data);
		vector<Mat> dataMat;
		allPosDataMat.push_back(dataMat);
	}

	//load pos data
	string posPath = "pos/"+dir+"加入水平翻轉/";
	
	ifstream ifs(posPath + "+000positive.txt", ios::in);
	int count = 0;	

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
			//cout << fileName << endl;
			vector<float> descriptorValue;
			allPosDataMat[count].push_back(img);

			_descriptor->compute(img, descriptorValue);			
			allPosData[count].push_back(descriptorValue);			
			count++;
			if (count == k)
				count = 0;			
		}
	}
	ifs.close();

	//cout << "Sum of data:" <<<< endl;
	for (int i = 0; i < k; i++)
	{		
		cout <<"Data "<<i<<"'s sum:"<< allPosData[i].size() << endl;
	}
	
	cout << "load POS data and divide data success" << endl;

	//load neg data
	string negPath = "neg/" + dir + "/";
	ifs.open(negPath + "000negative.txt", ios::in);
	if (!ifs.is_open())
	{
		throw exception("open NEG data fail");
	}
	int i = 0;
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		i++;
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt"&&i % 10 == 0)
		{
			Mat img = imread(negPath + fileName, IMREAD_GRAYSCALE);
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				throw exception("load NEG data fail");
			}

			vector<float> descriptorValue;

			if (dir == "行人") 
			{			
				for (int height = 0; height <img.rows - _hogFeature->winSize.height; height += _hogFeature->winSize.height / 2)
				{
					for (int width = 0; width<img.cols - _hogFeature->winSize.width; width += _hogFeature->winSize.width / 2)
					{
						Rect roi = Rect(width, height, _hogFeature->winSize.width, _hogFeature->winSize.height);
						_descriptor->compute(img(roi), descriptorValue);
						NegData.push_back(descriptorValue);
					}
				}
			}
			else 
			{
				_descriptor->compute(img, descriptorValue);
				NegData.push_back(descriptorValue);
			}			
		}
	}
	ifs.close();
	cout << "load NEG data success" << endl;
	

	//i作為測試資料 i以外的為訓練資料
	for (int i = 0; i < k; i++)
	{
		vector<vector<float>> trainData;
		//vector<vector<float>> testData;
		vector<Mat> testData;

		//分割測試集與訓練集
		int row = 0;
		for (int j = 0; j < k; j++)
		{
			if (j == i) 
			{
				testData = allPosDataMat[j];
				continue;
			}
				
			row += allPosData[j].size();
			trainData.insert(trainData.end(), allPosData[j].begin(), allPosData[j].end());
			
		}
		stringstream s;
		s << i;

		string featureName = "KFOLD//" + dir + "//"+s.str()+".xml";
		
		PrimalSVM* primalSVM;
		try {
			primalSVM = new PrimalSVM(featureName);
			cout << "start verify"<<featureName << endl;
			fstream fp;
			fp.open("KFOLD//" + dir + "//"+ dir +"KFOLD.txt", ios::app);//開啟檔案

			vector<float> hogVector;
			primalSVM->getSupportVector(hogVector);			
			_descriptor->setSVMDetector(hogVector);
			double sumOfMiss = 0.0;
			for (int j = 0; j<testData.size(); j++)
			{
				vector< Point > foundLocations;
				_descriptor->detect(testData[j], foundLocations);
				if (foundLocations.size() == 0) 
				{
					sumOfMiss+=1.0;
				}				
			}
			double accuracy = 1.0 - sumOfMiss / (double)testData.size();
			fp << i << "accuracy :" << accuracy << endl;
			cout <<i<<"accuracy :"<< accuracy << endl;
		}
		catch(Exception e)
		{
			
			//cout << trainData .size()<< endl;
			//cout << testData.size() << endl;
			int col = trainData[0].size();

			Mat positive(row, col, CV_32F);
			for (int j = 0; j < row; j++)
			{
				memcpy(&(positive.data[col*j * sizeof(float)]), trainData[j].data(), col * sizeof(float));
			}

			row = NegData.size();
			col = NegData[0].size();
			Mat negative(row, col, CV_32F);
			for (int i = 0; i < row; i++)
			{
				memcpy(&(negative.data[col*i * sizeof(float)]), NegData[i].data(), col * sizeof(float));
			}

			Mat trainingDataMat(positive.rows + negative.rows, positive.cols, CV_32FC1);
			//將正樣本資料放入mixData
			memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
			//將負樣本資料放入mixData
			memcpy(&trainingDataMat.data[sizeof(float)*positive.rows*positive.cols], negative.data, sizeof(float)*negative.rows*negative.cols);
			Mat dataProperty(positive.rows + negative.rows, 1, CV_32SC1, Scalar(-1.0));
			dataProperty.rowRange(0, positive.rows) = Scalar(1.0);
			std::cout << "Start to training" << std::endl;

			positive.release();
			negative.release();
			stringstream s;
			s << i;
			string filename = "KFOLD//" + dir + "//" + s.str();
			svm = ml::SVM::create();
			svm->setKernel(cv::ml::SVM::KernelTypes::LINEAR);
			Ptr<ml::TrainData> data = ml::TrainData::create(trainingDataMat, ml::SampleTypes::ROW_SAMPLE, dataProperty);
			svm->trainAuto(data);

			svm->save(filename + ".xml");
			writeTrainingParameter(filename);			
			std::cout << "svm training " << i << " is done!" << std::endl;
		}		
	}	
	
	std::cout << "svm training all done!" << std::endl;
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
			vector<float> descriptorValue;
			clock_t t;
			t = clock();
			mHog.calculateHogFeature(PosTrainingData[i], descriptorValue);
			t = clock() - t;
			printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);
			descriptorValueList.push_back(descriptorValue);
			imshow("test", PosTrainingData[i]);
			waitKey(0);
			/*vector<float> descriptorValue;
			double t = (double)cv::getTickCount();
			mHog.calculateHogFeature(PosTrainingData[i], descriptorValue);						
			t = (double)cv::getTickCount() - t;
			cout << "Spend Time:" << t << endl;
			descriptorValueList.push_back(descriptorValue);*/
			//cin.get();
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
		//將正樣本資料放入mixData
		memcpy(trainingDataMat.data, positive.data, sizeof(float)*positive.rows*positive.cols);
		//將負樣本資料放入mixData
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

void FeatureExtractor::verifyNewHOG()
{
	PrimalSVM* primalSVM;	
	
	primalSVM = new PrimalSVM("pedestrianFeature.xml");
	//primalSVM = new PrimalSVM("ped0626_auto.xml");
		
	vector<float> hogVector;
	primalSVM->getSupportVector(hogVector);
	_descriptor->setVersion(false);
	_descriptor->setCache(false);
	_descriptor->setSVMDetector(hogVector);

	//load pos data
	string path = "testPic/";
	ifstream ifs("compare/"+path + "000positive.txt", ios::in);
	//ifstream ifs("pos/moto/000positive.txt", ios::in);
	if (!ifs.is_open())
	{
		throw exception("open POS data fail");
	}
	int count = 0;
	double sum = 0.0;
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt")
		{
			Mat img = imread("compare/" + path + "/" + fileName);
			Mat grayImg;
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				throw exception("load POS data fail");
			}

			//cout << fileName << endl;			
			cvtColor(img, grayImg, CV_BGR2GRAY);
			vector<Rect> result;
			unsigned long start = clock();
			_descriptor->detectMultiScale(grayImg, result, 0.0, Size(8, 8), Size(), 1.05, 2.0,false);
			unsigned long end = clock();
			count++;
			sum += ((end - start) / 1000.0);
			printf("total time=%1.3f seconds\n", (end - start) / 1000.0);
			
			for (int i = 0; i<result.size(); i++)
			{
				cv::rectangle(img, result[i], cv::Scalar(255, 0, 0), 2);
			}			;
			//cv::imshow("image", img);
			//cv::waitKey(0);
			cv::imwrite("compare/preresult/"+fileName+".jpg", img);			
		}
	}
	ifs.close();	

}

void FeatureExtractor::verifyHelmet()
{
	enum detectType
	{
		HOGSVM,
		HOUGH
	};
		
	int type = HOGSVM;

	PrimalSVM* primalSVM;
	primalSVM = new PrimalSVM("helmet0818_auto.xml");
	
	vector<float> hogVector;
	primalSVM->getSupportVector(hogVector);
	_descriptor->setVersion(false);
	_descriptor->setCache(false);
	_descriptor->setSVMDetector(hogVector);


	//string path = "pos/機車正背面";
	//ifstream ifs(path + "/000positive.txt", ios::in);	
	//if (!ifs.is_open())
	//{
	//	throw exception("open POS data fail");
	//}
	//int count = 0;
	//double sum = 0.0;
	//while (!ifs.eof())
	//{
	//	string fileName;
	//	getline(ifs, fileName);
	//	if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt")
	//	{
	//		Mat img = imread(path + "/" + fileName);
	//		Mat grayImg;
	//		if (img.empty())
	//		{
	//			std::cout << fileName << " error" << std::endl;
	//			throw exception("load POS data fail");
	//		}
	//		//cout << fileName << endl;		

	//		

	//		resize(img, img, Size(), 1.4, 1.4);
	//		cvtColor(img, grayImg, CV_BGR2GRAY);

	//		Rect roi = checkROI(grayImg,"正背面");
	//		
	//		if(type==HOUGH)
	//		{
	//			Head head;				
	//			unsigned long start = clock();
	//			detectedHeadHoughCircles(grayImg(roi), &head);
	//			unsigned long end = clock();
	//			count++;
	//			sum += ((end - start) / 1000.0);
	//			circle(img(roi), head.center, head.radius, Scalar(0, 255, 0),3);
	//			printf("total time=%1.3f seconds\n", (end - start) / 1000.0);
	//		//	cv::rectangle(img, roi, cv::Scalar(255, 255, 0), 1);
	//			cv::imwrite("helmetResult/hough/" + fileName + ".jpg", img);
	//		}
	//		else  if(type == HOGSVM)
	//		{
	//			vector<Rect> result;
	//			vector<double> resultWeight;
	//			unsigned long start = clock();
	//			_descriptor->detectMultiScale(grayImg(roi), result, resultWeight, 0.8, Size(1, 1), Size(), 1.01, 2.0, false);				
	//			int Max = INT32_MIN;
	//			int maxIndex = -1;
	//			for (int i = 0; i<resultWeight.size(); i++)
	//			{
	//				if (resultWeight[i] > Max)
	//				{
	//					Max = resultWeight[i];
	//					maxIndex = i;
	//				}
	//			}
	//			
	//			unsigned long end = clock();
	//			count++;
	//			sum += ((end - start) / 1000.0);
	//			printf("total time=%1.3f seconds\n", (end - start) / 1000.0);				
	//			if (maxIndex != -1)
	//			{
	//				cv::rectangle(img(roi), result[maxIndex], cv::Scalar(255, 0, 0), 2);
	//			}
	//		//	cv::rectangle(img, roi, cv::Scalar(255, 255, 0), 1);
	//			//cv::imshow("image", img);
	//			//cv::waitKey(0);
	//			cv::imwrite("helmetResult/hogsvm/" + fileName + ".jpg", img);
	//		}			
	//	}
	//}
	//ifs.close();

	string path2 = "pos/機車側面全身";
	ifstream ifs2(path2 + "/000positive.txt", ios::in);
	if (!ifs2.is_open())
	{
		throw exception("open POS data fail");
	}
	
	while (!ifs2.eof())
	{
		string fileName;
		getline(ifs2, fileName);
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt")
		{
			Mat img = imread(path2 + "/" + fileName);
			Mat grayImg;
			if (img.empty())
			{
				cout << fileName << " error" << endl;
				throw exception("load POS data fail");
			}
			//cout << fileName << endl;		
			resize(img, img, Size(), 1.6, 1.6);
			cvtColor(img, grayImg, CV_BGR2GRAY);

			Rect roi = checkROI(grayImg, "側面");

			if (type == HOUGH)
			{
				Head head;				
				unsigned long start = clock();
				detectedHeadHoughCircles(grayImg(roi), &head);
				unsigned long end = clock();
				/*count++;
				sum += ((end - start) / 1000.0);*/
				circle(img(roi), head.center, head.radius, Scalar(0, 255, 0),3);
				printf("total time=%1.3f seconds\n", (end - start) / 1000.0);
			//	cv::rectangle(img, roi, cv::Scalar(255, 255, 0), 1);
				cv::imwrite("helmetResult/hough2/" + fileName + ".jpg", img);
			}
			else  if(type == HOGSVM)
			{
				vector<Rect> result;
				vector<double> resultWeight;
				unsigned long start = clock();
				_descriptor->detectMultiScale(grayImg(roi), result, resultWeight, 0.8, Size(1, 1), Size(), 1.01, 2.0, false);				
				int Max = INT32_MIN;
				int maxIndex = -1;
				for (int i = 0; i<resultWeight.size(); i++)
				{
					if (resultWeight[i] > Max)
					{
						Max = resultWeight[i];
						maxIndex = i;
					}
				}
				
				unsigned long end = clock();
				/*count++;
				sum += ((end - start) / 1000.0);*/
				printf("total time=%1.3f seconds\n", (end - start) / 1000.0);				
				if (maxIndex != -1)
				{
					cv::rectangle(img(roi), result[maxIndex], cv::Scalar(255, 0, 0), 2);
				}
			//	cv::rectangle(img, roi, cv::Scalar(255, 255, 0), 1);
				//cv::imshow("image", img);
				//cv::waitKey(0);
				cv::imwrite("helmetResult/hogsvm2/" + fileName + ".jpg", img);
			}			
		}
	}
	ifs2.close();

}

Rect FeatureExtractor::checkROI(Mat frame,string type)
{	
	int x;
	int y;
	int width;
	int height;

	if (type.compare("正背面") == 0)
	{
		x = 0;
		y = 0;
		width = frame.cols;
		height = frame.rows / 3;
	}
	if (type.compare("側面") == 0) 
	{
		x = frame.cols / 5;
		y = 0;
		width = frame.cols * 3 / 5;
		height = frame.rows*3/ 10;
	}	
	return Rect(x, y, width, height);
}

bool FeatureExtractor::detectedHeadHoughCircles(Mat grayFrame, Head * head)
{
	Mat sample = grayFrame;
	vector<cv::Vec3f> circles;
	double dp = 1.0;
	double minDist = sample.rows;
	if (sample.rows < sample.cols)
	{
		minDist = sample.cols;
	}
	double param1 = 100;
	double param2 = 15;
	int minRadius = sample.cols / 5;
	int maxRadius = sample.cols / 3;
	if (sample.cols > sample.rows)
	{
		minRadius = sample.rows / 5;
		maxRadius = sample.rows / 3;
	}
	HoughCircles(sample, circles, CV_HOUGH_GRADIENT, dp, minDist, param1, param2, minRadius, maxRadius);
	if (circles.size() != 0)
	{
		cv::Point center(cvRound(circles[0][0]), cvRound(circles[0][1]));
		int radius = cvRound(circles[0][2]);
		head->center = center;
		head->radius = radius;
		return true;
	}

	Mat sample2 = grayFrame;
	Mat dst;
	Sobel(sample2, dst, CV_8U, 1, 1);
	Mat abs_dst;
	convertScaleAbs(dst, abs_dst);
	Mat add_dst;
	addWeighted(sample2, 1, abs_dst, -1, CV_8U, add_dst);
	double dp2 = 1;
	double minDist2 = sample2.rows;
	if (sample2.rows < sample2.cols)
	{
		minDist2 = sample2.cols;
	}
	double param12 = 50;
	double param22 = 5;
	int minRadius2 = sample.cols / 5;
	int maxRadius2 = sample.cols / 3;
	if (sample.cols > sample.rows)
	{
		minRadius2 = sample.rows / 5;
		maxRadius2 = sample.rows / 3;
	}
	vector<cv::Vec3f> circles2;
	HoughCircles(add_dst, circles2, CV_HOUGH_GRADIENT, dp2, minDist2, param12, param22, minRadius2, maxRadius2);

	if (circles2.size() != 0)
	{
		cv::Point center(cvRound(circles2[0][0]), cvRound(circles2[0][1]));
		int radius = cvRound(circles2[0][2]);
		head->center = center;
		head->radius = radius;
		return true;
	}
	return false;
}