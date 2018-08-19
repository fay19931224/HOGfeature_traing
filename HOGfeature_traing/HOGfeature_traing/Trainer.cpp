#include "Trainer.h"

Trainer::Trainer(FeatureExtractor* featureExtractor)
{
	_featureExtractor = featureExtractor;
}


Trainer::~Trainer()
{
}

void Trainer::runTrainProcess()
{	
	//_featureExtractor->run();
	_featureExtractor->verifyHelmet();
	//_featureExtractor->verifyNewHOG();
	//_featureExtractor->flipImage();
	//_featureExtractor->KFoldCrossValidationForTrain(10);
}
