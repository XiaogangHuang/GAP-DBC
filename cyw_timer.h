#include <string.h>
#include <windows.h>
#include <iostream>
#ifndef CYW_TIMER_H_INCLUDED
#define CYW_TIMER_H_INCLUDED
#include <stdio.h>
#include <vector>
class CYW_TIMER
{
public:
    CYW_TIMER() {
        init_my_timer();
    }
    void my_strcat_double(char* buffer, double val) {
        char str_tmp[200];
        sprintf(str_tmp, "%f", val);

        strcat(buffer, str_tmp);
    }
    /**
     * Helper method for computing the current time (w.r.t to an offset).
     *
     *@return System in in microseconds
     */
    long get_system_time_in_microseconds() {
        /* --------RUN IN  LINUX----------
            struct timeval tempo;
            gettimeofday(&tempo, NULL);
            return tempo.tv_sec * 1000000 + tempo.tv_usec;
        */
        LARGE_INTEGER nFreq;
        LARGE_INTEGER t1;
        double dt;
        QueryPerformanceFrequency(&nFreq);
        QueryPerformanceCounter(&t1);
        dt = t1.QuadPart / (double)nFreq.QuadPart;
        return dt * 1000000;
    }

    double get_system_time_in_seconds() {
        /* --------RUN IN  LINUX----------
            struct timeval tempo;
            gettimeofday(&tempo, NULL);
            return tempo.tv_sec * 1000000 + tempo.tv_usec;
        */
        LARGE_INTEGER nFreq;
        LARGE_INTEGER t1;
        double dt;
        QueryPerformanceFrequency(&nFreq);
        QueryPerformanceCounter(&t1);
        dt = t1.QuadPart / (double)nFreq.QuadPart;
        return dt;
    }

    /**
     * Initializes a timer
     *
     * @param *timer Pointer to timer struct instance
     */
    void init_my_timer() {
        vec_time.clear();
        start_time = 0.0;
        elapsed_time = 0.0;
        elapsed_time_total = 0.0;
    }

    /**
     * Starts a given timer
     *
     * @param *timer Pointer to timer struct instance
     */
    void start_my_timer() {
        start_time = get_system_time_in_seconds();
    }

    /**
     * Stops a given timer
     *
     * @param *timer Pointer to timer struct instance
     */
    void stop_my_timer() {
        elapsed_time = (double)get_system_time_in_seconds() - start_time;
        elapsed_time_total += elapsed_time;
        vec_time.push_back(elapsed_time_total);
    }

    /**
     * Returns the time measured by a given timer
     *
     * @param *timer Pointer to timer struct instance
     * @return Passed time in seconds
     */
    double get_temp_timer() {
        return elapsed_time;
    }

    double get_my_timer() {
        return elapsed_time_total;
    }

    void print(char* out_info) {
        std::cout << out_info;
        print();
    }

    void print() {
        printf("Cumulative runtime: ");
        for (size_t i = 0; i < vec_time.size(); i++)
        {
            printf("%.4f ", vec_time[i]);
        }
        printf("\n");
    }

    void strcat_to_buffer(char* buffer) {
        strcat(buffer, "running time: ");
        my_strcat_double(buffer, elapsed_time_total);
    }

    void strcat_to_buffer(char* out_info, char* buffer) {
        strcat(buffer, out_info);
        my_strcat_double(buffer, elapsed_time_total);
    }
private:
    std::vector<double> vec_time;
    double start_time;
    double elapsed_time;
    double elapsed_time_total;
};


#endif // CYW_TIMER_H_INCLUDED
