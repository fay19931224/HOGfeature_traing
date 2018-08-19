#include "HogFeature.h"
#include <time.h>
# define M_PI           3.14159265

HogFeature::HogFeature()
{
}

HogFeature::HogFeature(HogFeatureParameter * hogFeatureParameter)
{
	_hogFeatureParameter = hogFeatureParameter;
	_radiusToDegree = (180 / M_PI);
}

HogFeature::~HogFeature()
{
}

int HogFeature::pixelSubtract(uchar firstPixel, uchar secPixel)
{
	int temp = firstPixel - secPixel;
	if (temp <= 10)
	{
		return 0;
	}
	else
	{
		return temp;
	}
}

void HogFeature::calculateResiduale_X(Mat &inputMat, Mat &y_Residuals)
{
	for (int y = 0; y < inputMat.rows; y++)
	{
		for (int x = 0; x < inputMat.cols; x++)
		{
			if (y == 0) //first row
			{
				y_Residuals.at<int>(y, x) = inputMat.at<uchar>(y + 1, x) - inputMat.at<uchar>(y, x);
				//y_Residuals.at<int>(y, x) = pixelSubtract(inputMat.at<uchar>(y + 1, x), inputMat.at<uchar>(y, x));
			}
			else if (y == (inputMat.rows - 1)) //last row
			{
				y_Residuals.at<int>(y, x) = inputMat.at<uchar>(y, x) - inputMat.at<uchar>(y - 1, x);
				//y_Residuals.at<int>(y, x) = pixelSubtract(inputMat.at<uchar>(y, x), inputMat.at<uchar>(y - 1, x));
			}
			else
			{
				y_Residuals.at<int>(y, x) = inputMat.at<uchar>(y + 1, x) - inputMat.at<uchar>(y - 1, x);
				//y_Residuals.at<int>(y, x) = pixelSubtract(inputMat.at<uchar>(y + 1, x), inputMat.at<uchar>(y - 1, x));
			}
		}
	}
}

void HogFeature::calculateResiduale_Y(Mat &inputMat, Mat &x_Residuals)
{
	for (int y = 0; y < inputMat.rows; y++)
	{
		for (int x = 0; x < inputMat.cols; x++)
		{
			if (x == 0) //first column
			{
				x_Residuals.at<int>(y, x) = inputMat.at<uchar>(y, x + 1) - inputMat.at<uchar>(y, x);
				//x_Residuals.at<int>(y, x) = pixelSubtract(inputMat.at<uchar>(y, x + 1), inputMat.at<uchar>(y, x));
			}
			else if (x == (inputMat.cols - 1)) //last column
			{
				x_Residuals.at<int>(y, x) = inputMat.at<uchar>(y, x) - inputMat.at<uchar>(y, x - 1);
				//x_Residuals.at<int>(y, x) = pixelSubtract(inputMat.at<uchar>(y, x), inputMat.at<uchar>(y, x - 1));
			}
			else
			{
				x_Residuals.at<int>(y, x) = inputMat.at<uchar>(y, x + 1) - inputMat.at<uchar>(y, x - 1);
				//x_Residuals.at<int>(y, x) = pixelSubtract(inputMat.at<uchar>(y, x + 1), inputMat.at<uchar>(y, x - 1));
			}
		}
	}
}

Mat HogFeature::calculateResidualedEdgeImage(Mat &inputMat)
{
	Mat x_Residuals(inputMat.rows, inputMat.cols, CV_32FC1, cvScalar(0));
	Mat y_Residuals(inputMat.rows, inputMat.cols, CV_32FC1, cvScalar(0));
	calculateResiduale_X(inputMat, x_Residuals);
	calculateResiduale_Y(inputMat, y_Residuals);
	Mat result(inputMat.rows, inputMat.cols, CV_32FC1, cvScalar(0));

	for (int y = 0; y < inputMat.rows; y++)
	{
		for (int x = 0; x < inputMat.cols; x++)
		{
			result.at<int>(y, x) = squareRoot(pow(y_Residuals.at<int>(y, x), 2) + pow(x_Residuals.at<int>(y, x), 2));
		}
	}

	return result;
}

//Calcaulate each cell orientation
void HogFeature::calculateOrientation(Mat &x_Residuals, Mat &y_Residuals, vector<float> &hogFeatureVector, int &cell_X, int &cell_Y, int &vectorIndex)
{	
	int total = 0;
	float tanOrientation0 = 0, tanOrientation1 = 0, tanOrientation2 = 0, tanOrientation3 = 0, tanOrientation4 = 0, tanOrientation5 = 0, tanOrientation6 = 0, tanOrientation7 = 0, tanOrientation8 = 0;
	float magnitude0 = 0, magnitude1 = 0, magnitude2 = 0, magnitude3 = 0, magnitude4 = 0, magnitude5 = 0, magnitude6 = 0, magnitude7 = 0, magnitude8 = 0;
	float magnitudeSum = 0;
	
	for (int innerY = 0; innerY < _hogFeatureParameter->cellSize.height; innerY++)
	{
		for (int innerX = 0; innerX <_hogFeatureParameter->cellSize.width; innerX++)
		{
			int xValue = x_Residuals.at<int>(cell_Y + innerY, cell_X + innerX);
			int yValue = y_Residuals.at<int>(cell_Y + innerY, cell_X + innerX);
			if (xValue != 0 && yValue != 0)
			{				
				int tan = myATan2(yValue, xValue) * _radiusToDegree;
				float magnitude = sqrtf(xValue * xValue + yValue * yValue);
				magnitudeSum += magnitude;

				if (tan < 0)
				{
					tan += 180;
				}
				//沒有將值依照權重放入?
				switch (tan / 20)
				{
				case 0:
					tanOrientation0++;
					magnitude0 += magnitude;
					total++;
					break;
				case 1:
					tanOrientation1++;
					magnitude1 += magnitude;
					total++;
					break;
				case 2:
					tanOrientation2++;
					magnitude2 += magnitude;
					total++;
					break;
				case 3:
					tanOrientation3++;
					magnitude3 += magnitude;
					total++;
					break;
				case 4:
					tanOrientation4++;
					magnitude4 += magnitude;
					total++;
					break;
				case 5:
					tanOrientation5++;
					magnitude5 += magnitude;
					total++;
					break;
				case 6:
					tanOrientation6++;
					magnitude6 += magnitude;
					total++;
					break;
				case 7:
					tanOrientation7++;
					magnitude7 += magnitude;
					total++;
					break;
				case 8:
					tanOrientation8++;
					magnitude8 += magnitude;
					total++;
					break;				
				}
			}
			else if (xValue == 0 && yValue != 0)
			{
				tanOrientation4++;
				float magnitude = sqrtf(xValue * xValue + yValue * yValue);
				magnitude4 += magnitude;
				magnitudeSum += magnitude;
				total++;
			}
			else if (xValue != 0 && yValue == 0)
			{
				tanOrientation0++;
				float magnitude = sqrtf(xValue * xValue + yValue * yValue);
				magnitude0 += magnitude;
				magnitudeSum += magnitude;
				total++;
			}
			else
			{
				tanOrientation0++;
				tanOrientation1++;
				tanOrientation2++;
				tanOrientation3++;
				tanOrientation4++;
				tanOrientation5++;
				tanOrientation6++;
				tanOrientation7++;
				tanOrientation8++;
				total++;

				magnitude0 = 0;
				magnitude1 = 0;
				magnitude2 = 0;
				magnitude3 = 0;
				magnitude4 = 0;
				magnitude5 = 0;
				magnitude6 = 0;
				magnitude7 = 0;
				magnitude8 = 0;
				magnitudeSum = 1;				
			}
		}
	}
	
	
	if (total != 0)
	{		
		hogFeatureVector[vectorIndex] = (magnitude0 / magnitudeSum);
		hogFeatureVector[vectorIndex + 1] = (magnitude1 / magnitudeSum);
		hogFeatureVector[vectorIndex + 2] = (magnitude2 / magnitudeSum);
		hogFeatureVector[vectorIndex + 3] = (magnitude3 / magnitudeSum);
		hogFeatureVector[vectorIndex + 4] = (magnitude4 / magnitudeSum);
		hogFeatureVector[vectorIndex + 5] = (magnitude5 / magnitudeSum);
		hogFeatureVector[vectorIndex + 6] = (magnitude6 / magnitudeSum);
		hogFeatureVector[vectorIndex + 7] = (magnitude7 / magnitudeSum);
		hogFeatureVector[vectorIndex + 8] = (magnitude8 / magnitudeSum);
	}
	else
	{
		for (int i = 0; i < _hogFeatureParameter->nbins; i++)
		{
			hogFeatureVector.push_back(0);
		}
	}
}

void HogFeature::normalizeHogFeature(Mat &inputMat, vector<float>& hogFeatureVector, vector<float> &normalizedFeature)
{
	//vector<float> normalizedHogFeature;
	//int size = _cellsInBlock * _cellsInBlock * _orientations;
	//float L2HysThreshold = 0.2;
	int normalizedFeatureIndex = 0;
	int widthStep = (inputMat.cols / _hogFeatureParameter->cellSize.width);	

	//每一次處理一個block
	for (int blockY = 0; blockY <= inputMat.rows / _hogFeatureParameter->cellSize.height - cellsInBlock; blockY++)
	{
		for (int blockX = 0; blockX <= inputMat.cols / _hogFeatureParameter->cellSize.width - cellsInBlock; blockX++)
		{			
			//歸一化範圍
			float blockSum = 0;
			for (int innerBlockY = 0; innerBlockY < cellsInBlock; innerBlockY++)
			{
				for (int innerBlockX = 0; innerBlockX < cellsInBlock; innerBlockX++)
				{					
					//cout << ((blockY + innerBlockY) * widthStep + (blockX + innerBlockX))*_hogFeatureParameter->nbins << endl;					
					for (int index = 0; index < _hogFeatureParameter->nbins; index++)
					{						
						float tempValue = hogFeatureVector.at(((blockY + innerBlockY) * widthStep + (blockX + innerBlockX)) *_hogFeatureParameter->nbins + index);
						blockSum += tempValue * tempValue;
					}					
				}
			}
			//cout << endl;
			float scale = 1.f / (sqrtf(blockSum) + 1e-3f);//改良區塊?

			//只對block裡的4個cell做歸一化?
			for (int innerBlockY = 0; innerBlockY < cellsInBlock; innerBlockY++)
			{
				for (int innerBlockX = 0; innerBlockX < cellsInBlock; innerBlockX++)
				{
					//cout <<"nor::"<< ((blockY + innerBlockY) * widthStep + (blockX + innerBlockX)) * _hogFeatureParameter->nbins << endl;
					for (int index = 0; index < _hogFeatureParameter->nbins; index++)
					{												
						normalizedFeature[normalizedFeatureIndex] = hogFeatureVector.at(((blockY + innerBlockY) * widthStep + (blockX + innerBlockX)) * _hogFeatureParameter->nbins + index) * scale;
						normalizedFeatureIndex++;
					}
				}
			}
		}
	}
}

void HogFeature::calculateHogFeature(Mat &inputMat, vector<float> &returnHogFeature)
{
	//Tools utility;
	Mat x_Residuals(inputMat.rows, inputMat.cols, CV_32FC1);
	Mat y_Residuals(inputMat.rows, inputMat.cols, CV_32FC1);
	calculateResiduale_X(inputMat, x_Residuals);
	calculateResiduale_Y(inputMat, y_Residuals);
	
	//hogFeatureVector的大小為block數乘上bin數，並用來存放每個block裡的直方圖資料，存放每個cell的梯度向量直方圖
	vector<float> hogFeatureVector((inputMat.rows / _hogFeatureParameter->cellSize.height) * (inputMat.cols / _hogFeatureParameter->cellSize.width) * _hogFeatureParameter->nbins, float());	
	int vectorIndex = 0;
	

	//計算每個cell的梯度向量
	for (int globalY = 0; globalY <= (inputMat.rows - _hogFeatureParameter->cellSize.height); globalY += _hogFeatureParameter->cellSize.height)
	{
		for (int globalX = 0; globalX <= (inputMat.cols - _hogFeatureParameter->cellSize.width); globalX += _hogFeatureParameter->cellSize.width)
		{
			//calculateOrientation求論文中的M0
			calculateOrientation(x_Residuals, y_Residuals, hogFeatureVector, globalX, globalY, vectorIndex);
			vectorIndex += _hogFeatureParameter->nbins;		
		}
	}
	
	int dim = ((inputMat.cols / _hogFeatureParameter->cellSize.width) - 1) * ((inputMat.rows / _hogFeatureParameter->cellSize.height) - 1) * cellsInBlock * cellsInBlock * _hogFeatureParameter->nbins;
	returnHogFeature.resize(dim);	
	normalizeHogFeature(inputMat, hogFeatureVector, returnHogFeature);
	
}

//?
void HogFeature::getWindowBlocks(Mat &frame, vector<float> &hogFeatures, int &index_X, int &index_Y, vector<float> &returnFeature, float &trainScale)
{
	int X_Blocks_Num = (_hogFeatureParameter->winSize.width* trainScale) / _hogFeatureParameter->cellSize.width;
	int Y_Blocks_Num = (_hogFeatureParameter->winSize.height * trainScale) / _hogFeatureParameter->cellSize.height;
	int stratX = index_X / _hogFeatureParameter->cellSize.width;
	int startY = index_Y / _hogFeatureParameter->cellSize.height;
	
	int widthStep = (frame.cols / _hogFeatureParameter->cellSize.width);
	int blockStep = cellsInBlock * cellsInBlock * _hogFeatureParameter->nbins;

	int returnFeatureIndex = 0;
	returnFeature.resize((X_Blocks_Num - 1) * (Y_Blocks_Num - 1) * blockStep);
	
	for (int yBlockPosition = startY; yBlockPosition <= startY + Y_Blocks_Num - cellsInBlock; yBlockPosition++)
	{
		for (int xBlockPosition = stratX; xBlockPosition <= stratX + X_Blocks_Num - cellsInBlock; xBlockPosition++)
		{
			for (int innerIndex = 0; innerIndex < blockStep; innerIndex++)
			{
				returnFeature[returnFeatureIndex] = hogFeatures[((yBlockPosition * (widthStep - 1) + xBlockPosition) * blockStep) + innerIndex];
				returnFeatureIndex++;
				//returnFeature.push_back(hogFeatures.at(((yBlockPosition * (widthStep - 1) + xBlockPosition) * _programAttributes.cellsInBlock * _programAttributes.cellsInBlock * _programAttributes.orientations) + innerIndex));
			}
		}
	}
}

float HogFeature::squareRoot(float x)
{
	unsigned int i = *(unsigned int*)&x;

	// adjust bias
	i += 127 << 23;
	// approximation of square root
	i >>= 1;

	return *(float*)&i;
}

float HogFeature::myATan2(int y, int x)
{	
	//此處為加速hog的重點吧
	float coeff_1 = M_PI / 4;
	float coeff_2 = 3 * coeff_1;
	float abs_y = abs(y);
	float angle;
	if (x >= 0)
	{
		float r = (x - abs_y) / (x + abs_y);
		angle = coeff_1 - coeff_1 * r;
	}
	else
	{
		float r = (x + abs_y) / (abs_y - x);
		angle = coeff_2 - coeff_1 * r;
	}
	return y < 0 ? -angle : angle;
}
