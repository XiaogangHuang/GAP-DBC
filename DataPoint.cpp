#include "DataPoint.h"

//Default constructor
DataPoint::DataPoint()
{
}

//constructor
DataPoint::DataPoint(int c) :clusterId(-1)
{
}

vector<int>& DataPoint::PCore()
{
	return this->pcore;
}

void DataPoint::SetPCore(int c)
{
	this->pcore.push_back(c);
}

int DataPoint::GetClusterId()
{
	return this->clusterId;
}

void DataPoint::SetClusterId(int clusterId)
{
	this->clusterId = clusterId;
}
