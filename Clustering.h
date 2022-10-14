#include "Toolfunction.h"

//Cluster analysis type.
class Clustering
{
private:
    vector< DataPoint > dadaSets;
    array2ddouble dataSets;
    kdtree2* kdtree;
    vector< vector< pair< int, int > > > ADJ;
    vector< vector< int > > PCOR;
    vector< vector< int > > NOI;
    vector< double > MAXDIST;
    DisjSet clust_struct;
    State* state;
    CORE_STAT* cores;
    int* index;
    double radius;
    double percent;
    int dataNum;
    int dataDim;
    int minPTs;
    int clusterId;
    int queries_num;
    CYW_TIMER timer;
public:
    Clustering() {}                //Default constructor.
    ~Clustering() {
        delete kdtree;
        free(state);
        free(cores);
        free(index); 
    }
    bool Init(const char* fileName, double radius, int minPTs, double p);    //Initialization operation
    bool DoDBSCANRecursive();

    void SetArrivalPoints_kd(int dpId);
    bool CheckCorePt_kd(double d, int cnt, int core1, int core2, vector<int>& mark);
    void BuildGraph_kd();
    void FindConnectedCom();
    bool StopCondition();
    void SelectObjects();
    void RangeQuery_kd();
    void MergeAndUpdate();
    void Merge_kd();
    bool Merge_fast_kd(int a, int b);
    bool Filter_kd(int cbs_1, int cbs_2, int* visited_cbs1,
        int* visited_cbs2, int& query);
    void AssignPoints_kd();

    double distance_euclidean(const vector<point_coord_type>& a,
        const vector<point_coord_type>& b);
    void PrintMessage();
    bool WriteToFile(const char* fileName);    //save results
};
