#include <iostream>
#include <omp.h>

int main() {
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        std::cout << "Hello World from thread " << threadID << std::endl;
    }
    return 0;
}
