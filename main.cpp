#include "guidedfilter.h"
#include <iostream>
#include <chrono>
#include <opencv2/opencv.hpp>

using namespace std;

int main(int argc, char* argv[])
{
    int r = 8, numThreads=atoi(argv[1]);
    double eps = 0.02 * 0.02 * 255 * 255;
    cv::Mat I = cv::imread("../cave-noflash.png", cv::IMREAD_COLOR);
    cv::Mat p = cv::imread("../cave-flash.png", cv::IMREAD_COLOR);
    // cv::cvtColor(I, I, cv::COLOR_BGR2GRAY);
    // cv::cvtColor(p, p, cv::COLOR_BGR2GRAY);
    cv::Mat q(p.cols, p.rows, p.type());
    AuxMono<uint8_t> *am = new AuxMono<uint8_t>(I.cols, I.rows, p.channels(), numThreads);
    AuxColor<uint8_t> *ac = new AuxColor<uint8_t>(I.cols, I.rows, p.channels(), I.channels(), numThreads);
    auto t1 = std::chrono::high_resolution_clock::now();
    if (I.channels() == 1)
        guidedFilterMono<uint8_t>(I.data, p.data, r, eps, q.data, am);
    else
        guidedFilterColor<uint8_t>(I.data, p.data, r, eps, q.data, ac);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms\n";
    cv::imwrite("../output.png", q);
    delete am, ac;
	return 0;
}