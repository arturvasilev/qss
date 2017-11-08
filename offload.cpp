#include <iostream>
#include <unistd.h>

#include <omp.h>
#include <nlopt.hpp>

#include "graph.hpp"

int main(int argc, char ** argv)
{
    using namespace std;

    if(argc == 0)
    {
        cout << "Enter name of config-file" << endl;
        return 1;
    }

    __attribute__((target(mic:0))) int myStaticInt = 13; // local var - available on MIC

    int tag;

#pragma offload target(mic:0) signal(&tag) inout(myStaticInt)
{
    cout << "\tLogical cores on mic: " << omp_get_num_procs() << endl;
    sleep(1);
    myStaticInt++;
    cout << "\tmain Variable: " << myStaticInt << endl;
}

myStaticInt++;

// #pragma offload_wait target(mic:0) wait(&tag)
// {
//     cout << "\tLogical cores on host: " << omp_get_num_procs() << endl;
//     cout << "\tmain Variable: " << myStaticInt << endl;
// }

    return 0;
}

