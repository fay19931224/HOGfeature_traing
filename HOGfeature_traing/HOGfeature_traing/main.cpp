#include <iostream>
#include <cstdlib> /* �üƬ������ */
#include <ctime>   /* �ɶ�������� */
#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <time.h>
#include "FeatureExtractor.h"
#include "Trainer.h"


void randomProduceSample()
{
	
	srand(time(NULL));
	int cnt = 0;
	
	
	ifstream ifs("sample/000.txt", ios::in);
	if (!ifs.is_open())
	{
		cout << "fail" << endl;
		system("PAUSE");	
	}
	while (!ifs.eof())
	{
		string fileName;
		getline(ifs, fileName);
		if (fileName != "" && fileName.substr(fileName.size() - 3, 3) != "txt") 
		{			
			Mat img = imread("sample/" + fileName, CV_LOAD_IMAGE_UNCHANGED);			
			if (img.empty())
			{
				cout << fileName << "error" << endl;
				system("PAUSE");
				continue;
			}
			/* ���w�üƽd�� */
			int min_x = 0; int max_x = img.cols - 1 - 72;
			int min_y = 0; int max_y = img.rows - 1 - 88;

			/* ���� [min , max] ����ƶü� */
			int x, y;
			for (int i = 0; i < 20; i++) {
				x = rand() % (max_x - min_x + 1) + min_x;
				y = rand() % (max_y - min_y + 1) + min_y;

				//std::cout << "x = " << x + 72 << ", y = " << y + 88 << std::endl;
				Mat ROI = img(cv::Rect(x, y, 72, 88));				
				stringstream ss;
				ss << cnt;
				imwrite("output\\" + ss.str() + ".jpg", ROI);
				cnt++;
			}
		}
	}	
}

int main()
{	
	//string trainingImformation = "HOG�Ѽ���/��������.txt";
	string trainingImformation = "HOG�Ѽ���/�������I��.txt";
	FeatureExtractor* extractor = new FeatureExtractor(trainingImformation);	
	Trainer trainer(extractor);

	extractor->ShowHOGFeature("pos/�������I��/videoData1202_2  (60).jpg");
	//extractor->ShowHOGFeature("pos/������������/videoData1202_6 (30).jpg");
	trainer.runTrainProcess(); 
	
	system("pause");
}
