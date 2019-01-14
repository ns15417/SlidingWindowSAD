#include <iostream>
#include <opencv2/opencv.hpp>
#include <random>
#include <ctime>
#include <vector>
#include<stdint.h>

#include"SAD.h"

using namespace cv;
using namespace std;

int main(int argc, char *argv[]) {

    if (3 != argc)
    {
        fprintf(stderr, "Usage: %s ref_image tgt_image\n", argv[0]);
        return -1;
    }

    const char *ref_file = argv[1];
    const char *tgt_file = argv[2];

    Mat ref_img = imread(ref_file, 0);  //8 bit png
    Mat tgt_img = imread(tgt_file, 0);  //8 bit png
    if (ref_img.empty() || tgt_img.empty())
    {
        return -1;
    }

    Mat ref = ref_img.clone();
    Mat tgt = tgt_img.clone();

    MATCH_PARAM para = {ref_img.cols, ref_img.rows, 1, 5, -127, 128, 0, 0, 0, 600, 900, 800 }; //x方向视差
    SAD sad_match(para);

    int w = para.width / para.scale;
    int h = para.height / para.scale;
    
    //计算cost
    sad_match.computerCost(ref, tgt);

    //获得最小代价函数和视差
    Mat1w m_minCost = Mat::zeros(Size(w, h), CV_16UC1);
    Mat1s m_disparity = Mat::zeros(Size(w, h), CV_16SC1);

    sad_match.getMinCostandDisparity(m_minCost, m_disparity);
    m_disparity = m_disparity + 50;
    imwrite("m_disparity.png", m_disparity);
    imwrite("m_minCost.png", m_minCost);

    return 0;

}