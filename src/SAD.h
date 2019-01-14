
#include"iostream"
#include"opencv2/opencv.hpp"
#include"iomanip"

using namespace std;
using namespace cv;
#define COST_MAX 65535

struct MATCH_PARAM{
    //图像大小
    int width;
    int height;

    //比例
    int scale;

    //匹配窗口大小
    int win_size;//25x25

    //搜索范围
    int disp_xmin;
    int disp_xmax;
    int disp_ymin;
    int disp_ymax;

    int disp_direction; //视差方向 0表示x方向, 1表示y方向

    //区域内点阈值
    int dot_num_min;
    int dot_num_max;

    //SAD阈值
    int SAD_th;
};

class SAD
{

public:

    SAD(MATCH_PARAM para);
    ~SAD();

    //计算SAD
    Mat computerSAD(const Mat &ref, const Mat &tgt); 

    // 计算代价函数
    void computerCost(const Mat &ref, const Mat &tgt); 

    //计算最小代价和视差（整像素视差）
    void getMinCostandDisparity(Mat1w &minCost, Mat1s &disparity);

private:
    MATCH_PARAM match_para;
    Mat *costVol;
    int DSR; //视差搜索区域的面积
    int DSR_width;//视差搜索区域的宽度
};
