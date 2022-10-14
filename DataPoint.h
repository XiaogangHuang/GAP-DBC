#include <vector>
using namespace std;

class DataPoint
{
private:
	vector<int> pcore;
	int clusterId;                    //the cluster ID
public:
	DataPoint();
	DataPoint(int c);
	vector<int>& PCore();
	void SetPCore(int c);
	int GetClusterId();
	void SetClusterId(int classId);
};

