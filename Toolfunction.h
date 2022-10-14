#pragma warning(disable:4996)
#pragma once
#include "kdtree2/kdtree2.hpp"
#include "DataPoint.h"
#include "cyw_timer.h"
#include "DisjSet.h"
#include <random>
#include <string.h>


#define FLT_MIN         1.175494351e-38F
#define MY_PI           3.14159265358979323846

struct State {
    int visited;
    int type;   // -1: uborder; 0: untouched; 1: pcore; 3: pborder; 4: pnoise;
    int nei_count;
    State() : visited(0), type(0), nei_count(0) {}
};

struct CORE_STAT {
    int un_border;
    int scanned;
    int clust;
    CORE_STAT() : un_border(0), scanned(0), clust(-1) {}
};

int get_dim(char* s, const char* delims);
double* get_data(char* s, int dim, const char* delims);
void read_data_dim_size(const char* filename, int* data_dim, int* data_size, const char* delims);
double* read_data(const char* filename, const char* delims);
double* read_data(const char* filename, const char* delims, int* dim, int* data_size);