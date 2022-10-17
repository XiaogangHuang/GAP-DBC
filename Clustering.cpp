#pragma warning(disable:4996)
#include "Clustering.h"

bool Clustering::Init(const char* fileName, double radius, int minPTs, double p) {
    this->radius = radius;        //set the radius
    this->percent = p;
    this->minPTs = minPTs;        //set the minPTs
    printf("Reading data...\n");

    double* raw_data = read_data(fileName, " ", &dataDim, &dataNum);
    state = (State*)malloc(sizeof(State) * dataNum);
    for (size_t i = 0; i < dataNum; i++)
    {
        state[i].visited = 0;
        state[i].type = 0;
        state[i].nei_count = 0;
    }
    index = (int*)malloc(dataNum * sizeof(int));
    dataSets.resize(dataNum);
    for (int i = 0; i < dataNum; i++) {
        DataPoint tempDP(-1);
        dadaSets.push_back(tempDP); 
        dataSets[i].resize(dataDim - 1);
        for (size_t j = 0; j < dataDim - 1; j++)
            dataSets[i][j] = raw_data[dataDim * i + j];
    }
    printf("n = %d  dim = %d\n", dataSets.size(), dataSets[0].size());

    timer.start_my_timer();
    kdtree = new kdtree2(dataSets);
    //kdtree = new kdtree2(dataSets, false);
    timer.stop_my_timer();
    printf("kd-tree build time = %.4f\n", timer.get_temp_timer());
    timer.init_my_timer();
    free(raw_data);

    return true;
}

bool Clustering::DoDBSCANRecursive() {
    timer.start_my_timer();
    vector<int> pos(dataNum);
    for (int i = 0; i < dataNum; i++)
        pos[i] = i;
    random_device rd;
    uint32_t seed = rd();
    //seed = 313360971;
    cout << "random seed = " << seed << endl;
    mt19937 rng(seed);
    shuffle(pos.begin(), pos.end(), rng);
    for (auto& it : pos)
    {
        if (state[it].type == 0)
        {
            queries_num++;
            SetArrivalPoints_kd(it);
        }
    }
    timer.stop_my_timer();
    printf("Step 1 runs in %.4fs\n", timer.get_temp_timer());
    printf("There are %zd Core Balls\n", PCOR.size());
    printf("There are %zd Non-core Balls\n", NOI.size());
    
    //Finding connected components
    if (PCOR.size() > 0)
    {
        timer.start_my_timer();
        cores = (CORE_STAT*)malloc(PCOR.size() * sizeof(CORE_STAT));
        for (size_t i = 0; i < PCOR.size(); i++)
            cores[i].scanned = 0;
        clust_struct.reSize(PCOR.size());
        printf("***********Initializing***********\n");
        BuildGraph_kd();
        for (size_t i = 0; i < dataNum; i++)
        {
            if (state[i].type > 0)
            {
                state[i].visited = 1;
            }
        }
        timer.stop_my_timer();
        printf("Number of clusters: %d\n", clust_struct.Size());
        printf("Number of queries: %d\n", queries_num);
        printf("Cumulative time: %.4f\n", timer.get_my_timer());
        printf("Starting iteration:\n");
        while (true)
        {
            timer.start_my_timer();
            FindConnectedCom();     //label nodes via their yes connected components
            if (StopCondition())    //Check stopping condition
                break;
            SelectObjects();        //Selecting objects and range queries
            MergeAndUpdate();       //Merging connected components, update graph 
            timer.stop_my_timer();
            printf("*************Refining*************\n");
            printf("Number of clusters: %d\n", clust_struct.Size());
            printf("Number of queries: %d\n", queries_num);
            printf("Cumulative time: %.4f\n", timer.get_my_timer());
        }
        printf("**********************************\n");
        printf("Merging...\n");
        timer.start_my_timer();
        Merge_kd();                    //Merging density-connected clusters
        timer.stop_my_timer();
        printf("Cumulative time: %.4f\n", timer.get_my_timer());
    }
    printf("Processing outliers...\n");
    timer.start_my_timer();
    AssignPoints_kd();          // Assign points to clusters
    timer.stop_my_timer();
    printf("Number of queries: %d\n", queries_num);
    printf("Cumulative time: %.4f\n", timer.get_my_timer());
    
    return true;
}

void Clustering::SetArrivalPoints_kd(int dpId) {
    kdtree2_result_vector res;
    kdtree->r_nearest_around_point(dpId, 0, radius * radius, res);
    int cout1 = res.size();
    state[dpId].nei_count = cout1;
    vector <int> pcor;
    vector <int> noi;
    if (cout1 >= minPTs)
    {
        state[dpId].type = 1;
    }
    else
    {
        state[dpId].type = 4;
    }
    if (cout1 >= minPTs) {
        int cnt = 0;
        double maxdist = 0.0;
        for (int i = 0; i < cout1; i++)
        {
            res[i].dis = sqrt(res[i].dis);
            if (2 * res[i].dis <= radius)
                cnt++;
        }
        pcor.push_back(dpId);
        for (int i = 0; i < cout1; i++)
        {
            int id = res[i].idx;
            double dist = res[i].dis;
            if (state[id].type <= 0)
            {
                state[id].nei_count++;
                if (dist <= radius * percent)
                {
                    state[id].type = -1;
                    if (state[id].nei_count >= minPTs)
                    {
                        state[id].type = 1;
                    }
                    else if (cnt >= minPTs && 2 * dist <= radius)
                    {
                        state[id].type = 1;
                    }
                    if (dist > maxdist)
                    {
                        maxdist = dist;
                    }
                    pcor.push_back(id);
                }
                else if (state[id].type == -1 && state[id].nei_count >= minPTs)
                {
                    state[id].type = 1;
                }
            }
            else if (state[id].type == 4)
            {
                state[id].type = 3;
            }
        }
        PCOR.push_back(pcor);
        MAXDIST.push_back(maxdist);
        pcor.clear();
    }
    else
    {
        noi.push_back(dpId);
        for (int i = 0; i < cout1; i++)
        {
            int id = res[i].idx;
            if (state[id].type <= 0)
            {
                state[id].nei_count++;
                if (state[id].nei_count >= minPTs && state[id].type == -1)
                    state[id].type = 1;
            }
            noi.push_back(id);
        }
        NOI.push_back(noi);
        noi.clear();
    }
}

bool Clustering::CheckCorePt_kd(double d, int cnt, int core1, int core2, vector<int>& mark)
{
    if (d > 2 * percent * radius)
    {
        return false;
    }
    vector< point_coord_type > midpoint(dataDim - 1);
    for (size_t i = 0; i < dataDim - 1; i++)
    {
        midpoint[i] = (dataSets[PCOR[core1][0]][i] + dataSets[PCOR[core2][0]][i]) / 2;
    }
    double radiustemp = radius * percent;
    double thresh = sqrt(radiustemp * radiustemp - d * d / 4);
    if (thresh < radius - thresh)
    {
        thresh = radius - thresh;
    }
    int pt = -1;
    for (size_t i = 0; i < PCOR[core2].size(); i++)
    {
        if (mark[PCOR[core2][i]] == 1 && distance_euclidean(midpoint,
            dataSets[PCOR[core2][i]]) <= radius - thresh)
        {
            pt = PCOR[core2][i];
            break;
        }
    }
    if (pt == -1)
    {
        return false;
    }
    else if (cnt >= minPTs)
    {
        for (size_t i = 0; i < PCOR[core2].size(); i++)
        {
            if (distance_euclidean(midpoint, dataSets[PCOR[core2][i]])
                <= radius - thresh)
            {
                state[PCOR[core2][i]].type = 1;
            }
        }
        return true;
    }
    for (size_t i = 0; i < PCOR[core2].size(); i++)
    {
        if (mark[PCOR[core2][i]] == 1)
        {
            mark[PCOR[core2][i]] = 2;
        }
        else if (distance_euclidean(midpoint, dataSets[PCOR[core2][i]]) <= thresh)
        {
            cnt++;
        }
    }
    if (cnt >= minPTs)
    {
        for (size_t i = 0; i < PCOR[core2].size(); i++)
        {
            if (distance_euclidean(midpoint, dataSets[PCOR[core2][i]])
                <= radius - thresh)
            {
                state[PCOR[core2][i]].type = 1;
            }
        }
        for (size_t i = 0; i < PCOR[core2].size(); i++)
        {
            if (mark[PCOR[core2][i]] == 2)
            {
                mark[PCOR[core2][i]] = 1;
            }
        }
        return true;
    }
    for (size_t i = 0; i < PCOR[core1].size(); i++)
    {
        if (mark[PCOR[core1][i]] == 1 &&
            distance_euclidean(midpoint, dataSets[PCOR[core1][i]]) <= thresh)
        {
            cnt++;
        }
    }
    for (size_t i = 0; i < PCOR[core2].size(); i++)
    {
        if (mark[PCOR[core2][i]] == 2)
        {
            mark[PCOR[core2][i]] = 1;
        }
    }
    if (cnt >= minPTs)
    {
        for (size_t i = 0; i < PCOR[core2].size(); i++)
        {
            if (distance_euclidean(midpoint, dataSets[PCOR[core2][i]])
                <= radius - thresh)
            {
                state[PCOR[core2][i]].type = 1;
            }
        }
        return true;
    }
    return false;
}

void Clustering::BuildGraph_kd()
{
    for (int jt = 0; jt < PCOR.size(); jt++) {
        for (int p = 0; p < PCOR[jt].size(); p++)
        {
            dadaSets[PCOR[jt][p]].SetPCore(jt);
        }
    }
    for (size_t i = 0; i < dataNum; i++)
    {
        if (state[i].type == 1)
        {
            vector<int>& vec = dadaSets[i].PCore();
            for (size_t j = 1; j < vec.size(); j++)
            {
                int ni = vec[j];
                int nj = vec[j - 1];
                clust_struct.Union(ni, nj);
            }
        }
    }
    array2ddouble pcores;
    pcores.resize(PCOR.size());
    for (int jt = 0; jt < PCOR.size(); jt++) {
        pcores[jt].resize(dataDim - 1);
        for (int p = 0; p < dataDim - 1; p++)
        {
            pcores[jt][p] = dataSets[PCOR[jt][0]][p];
        }
    }
    kdtree2 pctree(pcores);
    kdtree2_result_vector matches;
    matches.reserve(dataNum);
    double radiustemp = radius + 2 * percent * radius;
    vector< int > mark(dataNum, 0);
    ADJ.resize(PCOR.size());
    for (int it = 0; it < PCOR.size(); it++)
    {
        for (int i = 0; i < PCOR[it].size(); i++) {
            mark[PCOR[it][i]] = 1;
        }
        matches.resize(0);
        pctree.r_nearest_around_point(it, 0, radiustemp * radiustemp, matches);
        for (int i = 0; i < matches.size(); i++)
        {
            int id = matches[i].idx;
            double dist = sqrt(matches[i].dis);
            if (it < id && clust_struct.Is_same(it, id) == false
                && dist - radius <= MAXDIST[it] + MAXDIST[id])
            {
                int edge_state = 1;
                if (dist <= radius)
                {
                    edge_state = 0;
                    clust_struct.Union(it, id);
                }
                else if (dist <= 2 * radius)
                {
                    int intersect = 0;
                    for (int p = 0; p < PCOR[id].size(); p++)
                    {
                        if (mark[PCOR[id][p]] == 1)
                        {
                            intersect++;
                            if (state[PCOR[id][p]].type == 1)
                            {
                                intersect = -1;
                                break;
                            }
                        }
                    }
                    if (intersect == -1 ||
                        CheckCorePt_kd(dist, intersect, it, id, mark))
                    {
                        edge_state = 0;
                        clust_struct.Union(it, id);
                    }
                    else if (intersect > 0)
                    {
                        edge_state = 4;
                    }
                }
                ADJ[it].push_back({ id, edge_state });
                ADJ[id].push_back({ it, edge_state });
            }
        }
        for (int i = 0; i < PCOR[it].size(); i++) {
            mark[PCOR[it][i]] = 0;
        }
    }
}

void Clustering::FindConnectedCom()
{
    for (size_t i = 0; i < PCOR.size(); i++)
        cores[i].clust = -1;
    clusterId = 0;
    for (size_t i = 0; i < PCOR.size(); i++)
    {
        int p = clust_struct.Find(i);
        if (cores[p].clust == -1)
        {
            cores[p].clust = clusterId++;
        }
        cores[i].clust = cores[p].clust;
    }
}

bool Clustering::StopCondition()
{
    for (size_t i = 0; i < ADJ.size(); i++)
    {
        for (size_t j = 0; j < ADJ[i].size(); j++)
        {
            int nei = ADJ[i][j].first;
            int edgetype = ADJ[i][j].second;
            if (cores[i].clust != cores[nei].clust && edgetype != 0)
            {
                return false;
            }
        }
    }
    return true;
}

void Clustering::SelectObjects()
{
    for (size_t i = 0; i < PCOR.size(); i++)
    {
        cores[i].un_border = 0;
        if (cores[i].scanned == 0)
        {
            int unprocessed = 0;
            for (size_t j = 0; j < PCOR[i].size(); j++)
            {
                unprocessed += 1 - state[PCOR[i][j]].visited;
            }
            cores[i].un_border = unprocessed;
            if (unprocessed == 0)
                cores[i].scanned = 1;
        }
    }
    for (size_t i = 0; i < dataNum; i++)
    {
        if (state[i].visited == 0)
        {
            index[i] = 0;
        }
        else
        {
            index[i] = 1;
        }
    }
    RangeQuery_kd();
}

void Clustering::RangeQuery_kd()
{
    vector<double> pt(dataDim - 1);
    for (size_t dpId = 0; dpId < dataNum; dpId++)
    {
        if (index[dpId] == 0)
        {
            index[dpId] = 1;
            state[dpId].visited = 1;
            queries_num++;
            kdtree2_result_vector res;
            kdtree->r_nearest_around_point(dpId, 0, radius * radius, res);
            int cout1 = res.size();
            state[dpId].nei_count = cout1;
            if (cout1 >= minPTs) {
                int cnt = 0;
                for (int i = 0; i < cout1; i++)
                {
                    double dist = sqrt(res[i].dis);
                    if (2 * dist <= radius)
                        cnt++;
                }
                state[dpId].type = 1;
                vector<int>& pcore = dadaSets[dpId].PCore();
                for (size_t j = 1; j < pcore.size(); j++)
                {
                    if (cores[pcore[j]].clust != cores[pcore[j - 1]].clust)
                    {
                        clust_struct.Union(pcore[j], pcore[j - 1]);
                    }
                }
                for (int i = 0; i < cout1; i++)
                {
                    int id = res[i].idx;
                    index[id] = 1;
                    if (state[id].type < 0)
                    {
                        state[id].nei_count++;
                        if (state[id].nei_count >= minPTs)
                        {
                            state[id].visited = 1;
                            state[id].type = 1;
                            vector<int>& pcore1 = dadaSets[id].PCore();
                            for (size_t j = 1; j < pcore1.size(); j++)
                            {
                                if (cores[pcore1[j]].clust != cores[pcore1[j - 1]].clust)
                                {
                                    clust_struct.Union(pcore1[j], pcore1[j - 1]);
                                }
                            }
                            clust_struct.Union(pcore[0], pcore1[0]);
                        }
                        else if (cnt >= minPTs)
                        {
                            double dist = sqrt(res[i].dis);
                            if (2 * dist <= radius)
                            {
                                state[id].visited = 1;
                                state[id].type = 1;
                                vector<int>& pcore1 = dadaSets[id].PCore();
                                for (size_t j = 1; j < pcore1.size(); j++)
                                {
                                    if (cores[pcore1[j]].clust != cores[pcore1[j - 1]].clust)
                                    {
                                        clust_struct.Union(pcore1[j], pcore1[j - 1]);
                                    }
                                }
                                clust_struct.Union(pcore[0], pcore1[0]);
                            }
                        }
                    }
                    else if (state[id].type == 4)
                    {
                        state[id].type = 3;
                    }
                }
            }
            else
            {
                state[dpId].type = 4;
                for (int i = 0; i < cout1; i++)
                {
                    int id = res[i].idx;
                    if (state[id].type < 0)
                    {
                        state[id].nei_count++;
                        if (state[id].nei_count >= minPTs)
                        {
                            state[id].visited = 1;
                            state[id].type = 1;
                            index[id] = 1;
                            vector<int>& pcore1 = dadaSets[id].PCore();
                            for (size_t j = 1; j < pcore1.size(); j++)
                            {
                                if (cores[pcore1[j]].clust != cores[pcore1[j - 1]].clust)
                                {
                                    clust_struct.Union(pcore1[j], pcore1[j - 1]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Clustering::MergeAndUpdate()
{
    for (size_t i = 0; i < ADJ.size(); i++)
    {
        for (size_t j = 0; j < ADJ[i].size(); j++)
        {
            int nei = ADJ[i][j].first;
            int& edgetype = ADJ[i][j].second;
            if (edgetype > 0 && cores[i].un_border + cores[nei].un_border == 0)
            {
                edgetype = 0;
            }
        }

        if (cores[i].scanned == 0)
        {
            bool cond = true;
            for (size_t j = 0; j < ADJ[i].size(); j++)
            {
                int nei = ADJ[i][j].first;
                int edgetype = ADJ[i][j].second;
                if (clust_struct.Is_same(i, nei) == false && edgetype != 0)
                {
                    cond = false;
                    break;
                }
            }
            if (cond)
            {
                for (size_t it = 0; it < PCOR[i].size(); it++)
                {
                    state[PCOR[i][it]].visited = 1;
                }
            }
        }
    }
}

void Clustering::Merge_kd()
{
    /*for (size_t i = 0; i < PCOR.size(); i++)
    {
        double maxdist = 0.0;
        for (size_t j = 0; j < PCOR[i].size(); j++)
        {
            if (state[PCOR[i][j]].type == 1)
            {
                double dist = distance_euclidean(dataSets[PCOR[i][0]],
                    dataSets[PCOR[i][j]]);
                if (dist > maxdist)
                {
                    maxdist = dist;
                }
            }
        }
        MAXDIST[i] = maxdist;
    }*/
    for (size_t i = 0; i < ADJ.size(); i++)
    {
        for (size_t j = 0; j < ADJ[i].size(); j++)
        {
            int nei = ADJ[i][j].first;
            if (i < nei && clust_struct.Is_same(i, nei) == false)
            {
                if (Merge_fast_kd(i, nei))
                {
                    clust_struct.Union(i, nei);
                }
            }
        }
    }
}

bool Clustering::Merge_fast_kd(int cbs_1, int cbs_2) {
    int n1 = PCOR[cbs_1].size();
    int* visited_cbs1 = (int*)calloc(n1 + 1, sizeof(int));
    int cnt = n1;
    int query;
    for (size_t i = 0; i < PCOR[cbs_1].size(); i++)
    {
        if (state[PCOR[cbs_1][i]].type != 1)
        {
            visited_cbs1[i] = 1;
            cnt--;
        }
        else if (distance_euclidean(dataSets[PCOR[cbs_2][0]], 
            dataSets[PCOR[cbs_1][i]]) > MAXDIST[cbs_2] + radius)
        {
            visited_cbs1[i] = 1;
            cnt--;
        }
        else
        {
            query = i;
        }
    }
    visited_cbs1[n1] = cnt;
    int n2 = PCOR[cbs_2].size();
    int* visited_cbs2 = (int*)calloc(n2 + 1, sizeof(int));
    cnt = n2;
    for (size_t i = 0; i < PCOR[cbs_2].size(); i++)
    {
        if (state[PCOR[cbs_2][i]].type != 1)
        {
            visited_cbs2[i] = 1;
            cnt--;
        }
        else if (distance_euclidean(dataSets[PCOR[cbs_1][0]], 
            dataSets[PCOR[cbs_2][i]]) > MAXDIST[cbs_1] + radius)
        {
            visited_cbs2[i] = 1;
            cnt--;
        }
    }
    visited_cbs2[n2] = cnt;
    while (visited_cbs1[n1] > 20 && visited_cbs2[n2] > 20)
    {
        if (Filter_kd(cbs_1, cbs_2, visited_cbs1, visited_cbs2, query) ||
            Filter_kd(cbs_2, cbs_1, visited_cbs2, visited_cbs1, query))
        {
            free(visited_cbs1);
            free(visited_cbs2);
            return true;
        }
    }
    for (int i = 0; i < n2; i++) {
        if (visited_cbs2[i] == 0)
        {
            for (int j = 0; j < n1; j++) {
                if (visited_cbs1[j] == 0)
                {
                    double dist = distance_euclidean(dataSets[PCOR[cbs_1][j]],
                        dataSets[PCOR[cbs_2][i]]);
                    if (dist <= radius)
                    {
                        free(visited_cbs1);
                        free(visited_cbs2);
                        return true;
                    }
                }
            }
        }
    }
    free(visited_cbs1);
    free(visited_cbs2);
    return false;
}

bool Clustering::Filter_kd(int cbs_1, int cbs_2, int* visited_cbs1,
    int* visited_cbs2, int& query)
{
    int q = 0;
    double* distTemp = (double*)calloc(PCOR[cbs_2].size(), sizeof(double));
    double minidist = FLT_MAX;
    visited_cbs1[query] = 1;
    visited_cbs1[PCOR[cbs_1].size()] -= 1;
    for (int i = 0; i < PCOR[cbs_2].size(); i++) {
        if (visited_cbs2[i] == 0)
        {
            distTemp[i] = distance_euclidean(dataSets[PCOR[cbs_1][query]],
                dataSets[PCOR[cbs_2][i]]);
            if (distTemp[i] <= radius)
                return true;
            else if (distTemp[i] < minidist)
            {
                minidist = distTemp[i];
                q = i;
            }
        }
    }
    vector<double> point(dataDim - 1);
    double temp = 0;
    for (int i = 0; i < dataDim - 1; i++) {
        double a = dataSets[PCOR[cbs_1][query]][i];
        point[i] = a - dataSets[PCOR[cbs_2][q]][i];
        temp += a * point[i];
    }
    double theta = 0;
    for (int i = 0; i < PCOR[cbs_2].size(); i++) {
        if (visited_cbs2[i] == 0)
        {
            double alpha = temp;
            if (i == q)
            {
                alpha = asin(radius / distTemp[i]);
            }
            else
            {
                for (int j = 0; j < dataDim - 1; j++) {
                    alpha -= dataSets[PCOR[cbs_2][i]][j] * point[j];
                }
                alpha = acos(alpha / (distTemp[i] * minidist)) +
                    asin(radius / distTemp[i]);
            }
            if (alpha > theta)
            {
                theta = alpha;
                if (theta >= MY_PI)
                {
                    theta = MY_PI;
                    break;
                }
            }
        }
    }
    if (theta == MY_PI)
    {
        for (int i = 0; i < PCOR[cbs_1].size(); i++) {
            if (visited_cbs1[i] == 0)
            {
                double dist = distance_euclidean(dataSets[PCOR[cbs_1][query]],
                    dataSets[PCOR[cbs_1][i]]);
                if (dist < minidist - radius)
                {
                    visited_cbs1[i] = 1;
                    --visited_cbs1[PCOR[cbs_1].size()];
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < PCOR[cbs_1].size(); i++) {
            if (visited_cbs1[i] == 0)
            {
                double dist = distance_euclidean(dataSets[PCOR[cbs_1][query]],
                    dataSets[PCOR[cbs_1][i]]);
                double alpha = temp;
                for (int j = 0; j < dataDim - 1; j++) {
                    alpha -= dataSets[PCOR[cbs_1][i]][j] * point[j];
                }
                if (dist < minidist - radius ||
                    acos(alpha / (dist * minidist)) > theta)
                {
                    visited_cbs1[i] = 1;
                    --visited_cbs1[PCOR[cbs_1].size()];
                }
            }
        }
    }
    query = q;
    free(distTemp);
    return false;
}

void Clustering::AssignPoints_kd()
{
    for (size_t i = 0; i < PCOR.size(); i++)
        cores[i].clust = -1;
    clusterId = 0;
    for (size_t i = 0; i < PCOR.size(); i++)
    {
        int p = clust_struct.Find(i);
        if (cores[p].clust == -1)
        {
            cores[p].clust = clusterId++;
        }
        for (size_t j = 0; j < PCOR[i].size(); j++)
        {
            dadaSets[PCOR[i][j]].SetClusterId(cores[p].clust);
        }
    }
    vector<double> pt(dataDim - 1);
    for (int i = 0; i < NOI.size(); i++)
    {
        if (dadaSets[NOI[i][0]].GetClusterId() == -1)
        {
            for (int j = 1; j < NOI[i].size(); j++)
            {
                if (state[NOI[i][j]].type == 1)
                {
                    int c = dadaSets[NOI[i][j]].GetClusterId();
                    dadaSets[NOI[i][0]].SetClusterId(c);
                    break;
                }
            }
        }
        if (dadaSets[NOI[i][0]].GetClusterId() == -1)
        {
            for (int j = 1; j < NOI[i].size(); j++)
            {
                int c = dadaSets[NOI[i][j]].GetClusterId();
                if (state[NOI[i][j]].type == -1)
                {
                    kdtree2_result_vector res;
                    kdtree->r_nearest_around_point(NOI[i][j], 0, radius * radius, res);
                    queries_num++;
                    if (res.size() >= minPTs)
                    {
                        state[NOI[i][j]].type = 1;
                        dadaSets[NOI[i][0]].SetClusterId(c);
                        break;
                    }
                    else
                    {
                        state[NOI[i][j]].type = 3;
                    }
                }
            }
        }
    }
}

/** Euclidean distance between two vector **/
double Clustering::distance_euclidean(const vector<point_coord_type>& a,
    const vector<point_coord_type>& b) {
    int dim = a.size();
    double d = 0.0;
    for (int i = 0; i < dim; ++i) {
        double t = (a[i] - b[i]);
        d += t * t;
    }
    return sqrt(d);
}

void Clustering::PrintMessage()
{
    printf("******************************************\n");
    int* temp = (int*)calloc(clusterId, sizeof(int));
    int noisePts = 0;
    for (int i = 0; i < dataNum; i++)
    {
        if (dadaSets[i].GetClusterId() != -1)
        {
            temp[dadaSets[i].GetClusterId()] += 1;
        }
        else
            noisePts++;
    }
    printf("Number of clusters: %d\n", clusterId);
    printf("Outlier: %d pts\n", noisePts);
    printf("Total running time = %.4fs\n", timer.get_my_timer());
    timer.print();

    free(temp);
}

bool Clustering::WriteToFile(const char* fileName) {
    ofstream of1(fileName);
    for (int i = 0; i < dataNum; i++)
    {
        of1 << dadaSets[i].GetClusterId() << ' ';
    }
    of1 << endl;
    of1.close();
    return true;
}
