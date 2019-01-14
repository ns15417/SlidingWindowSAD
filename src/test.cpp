
#include <iostream>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"


#include "sliding_window_sample_code.hpp"

using namespace cv;
using namespace std;

int main(int argc, char* argv[])
{ 
    if(argc!=3)
    {
        cout<<"Usage: ./main  leftimg  rightimg!"<<endl;
        return -1;
    }else
    {
        Mat leftimage = imread(argv[1],0);
        Mat rightimage = imread(argv[2],0);


        int block_size = 5;
        int min_disp = -127;
        int max_disp = 128;
        int disp_range = 256;
        int width = leftimage.cols;
        int height = rightimage.rows;


        Mat tgtimage;
        Mat refimage;

        copyMakeBorder(leftimage, tgtimage,2,2,129,130, BORDER_CONSTANT, Scalar(0));
        copyMakeBorder(rightimage, refimage,2,2,129,130, BORDER_CONSTANT, Scalar(0));

        imwrite("leftimage.png",tgtimage);
        imwrite("rightimage.png",tgtimage);

        DispSlidingWindow  dispSlidingWindow;


        dispSlidingWindow.setParams(block_size, min_disp,
                    disp_range, width, height, max_disp);

        
        DEImageU8 ext_tgt;
        ext_tgt.height = tgtimage.rows;
        ext_tgt.width = tgtimage.cols;
        ext_tgt.data = tgtimage.data;

        DEImageU8 ext_ref;
        ext_ref.height = refimage.rows;
        ext_ref.width = refimage.cols;
        ext_ref.data = refimage.data;


        int cols = leftimage.cols;
        int rows = leftimage.rows;

        Mat dispImg(rows,cols, CV_16SC1,Scalar(0));
        DEImageS16 disp;
        disp.height = dispImg.rows;
        disp.width = dispImg.cols;
        disp.data = dispImg.ptr<int16_t>();
        
        Mat costImg(rows,cols, CV_16UC1,Scalar(0));
        DEImageU16 cost;
        cost.height = costImg.rows;
        cost.width = costImg.cols;
        cost.data = costImg.ptr<uint16_t>();

        dispSlidingWindow.SlidingWindowMatch(ext_tgt,ext_ref,disp,cost);

        cv::imwrite("costImg.png",costImg);

        //调整视差范围
        Mat dispoutImg(rows,cols, CV_16SC1,Scalar(0));
        dispoutImg = abs(dispImg-127);
        cv::imwrite("dispout.png",dispoutImg);


        return 0;
    }
}