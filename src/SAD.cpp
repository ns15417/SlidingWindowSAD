#include"SAD.h"

using namespace std;
using namespace cv;

SAD::SAD(MATCH_PARAM para):match_para(para)
{    
    DSR = (match_para.disp_xmax - match_para.disp_xmin + 1) * (match_para.disp_ymax - match_para.disp_ymin + 1);
    DSR_width = match_para.disp_xmax - match_para.disp_xmin + 1;
    
    int cost_width = match_para.width / match_para.scale;
    int cost_height = match_para.height / match_para.scale;

    costVol = new Mat[DSR];
    
    for (int mIdx = 0; mIdx < DSR; mIdx++) {
        costVol[mIdx] = Mat(Size(cost_width, cost_height), CV_16UC1, Scalar::all(COST_MAX));
    }
}

SAD::~SAD()
{
    if (costVol != nullptr)
    {
        delete [] costVol;
        costVol = nullptr;
    }
}

void SAD::computerCost(const Mat &ref, const Mat &tgt)
{
    int height = ref.rows;
    int width = ref.cols;

    int winSize = match_para.win_size;
    int half_size = winSize / 2;

    Mat ref_region(Size(winSize, winSize), CV_8U, Scalar::all(0));
    Mat tgt_region(Size(winSize, winSize), CV_8U, Scalar::all(0));

    //以tgt图像为基准在ref上的匹配点
    for (int j = half_size; j < height - half_size; j = j + match_para.scale)
    {
        int index_j = j / match_para.scale;
        for (int i = half_size; i< width - half_size; i = i + match_para.scale)
        {
            //先求以当前点为中心的窗口内灰度和，避免全黑或者全白的区域
            tgt_region = tgt(Rect(i - half_size, j - half_size, winSize, winSize));
            Scalar sum_region = sum(tgt_region);
        
            //这里的阈值暂时还未找到合适的值
            //if (sum_region(0) > match_para.dot_num_max || sum_region(0) < match_para.dot_num_min)
            //{
            //    continue;
            //}

            //在搜索范围内搜索
            int ref_start_x = i + match_para.disp_xmin;
            int ref_end_x = i + match_para.disp_xmax;
            int ref_start_y = j + match_para.disp_ymin;
            int ref_end_y = j + match_para.disp_ymax;
            int d = 0;
            int index_i = i / match_para.scale;

            //在ref的区域内搜索
            for (int ref_y = ref_start_y; ref_y <= ref_end_y; ref_y++)
            {            
                for (int ref_x = ref_start_x; ref_x <= ref_end_x; ref_x++)
                {
                    ushort* cost = (ushort*)costVol[d].ptr<ushort>(index_j);
                    d++;

                    if (ref_x - half_size < 0 || ref_y - half_size < 0 || ref_x + half_size >= (ref.cols -1) || ref_y + half_size >= (ref.rows -1))
                    {
                        cost[index_i] = 65535;
                        continue;
                    }

                    ref_region = ref(Rect(ref_x - half_size, ref_y - half_size, winSize, winSize));

                    Mat Dif;
                    absdiff(tgt_region, ref_region, Dif);
                    Scalar ADD = sum(Dif);
                    cost[index_i] = ADD[0];
                     
                }
            }
        }   
    }

}

void SAD::getMinCostandDisparity(Mat1w &minCost, Mat1s &disparity)
{
    if (minCost.empty() || disparity.empty())
    {
        return;
    }

    for (int y = 0; y < minCost.rows - 0; y++) {
        ushort *min_cost = minCost.ptr<ushort>(y);
        short *disp = disparity.ptr<short>(y);

        for (int x = 0; x < minCost.cols - 0; x++) {
            ushort mincost = COST_MAX;
            int minDis = 0;

            for (int d = 0; d < DSR; d++) {
                const ushort *costData = costVol[d].ptr<ushort>(y);

                if ((costData[x] < mincost) && (minDis != 200)) 
                {
                    mincost = costData[x];
                    minDis = d;
                }
            }
            vector<cv::Point> points;
            for (int d = 0; d < DSR; d++) 
            {
                const ushort *costData = costVol[d].ptr<ushort>(y);
                ushort cost = costData[x];
                cv::Point tmp_P;
                tmp_P.x = d * 10;
                tmp_P.y = cost;
                points.push_back(tmp_P);         
            }

            //当最小代价满足一定的阈值条件才是有效的
            //if (mincost < match_para.SAD_th)
            {
                min_cost[x] = mincost;
                if(SAD::match_para.disp_direction == 0)  //视差为 x 方向
                {
                    disp[x] = -(match_para.disp_xmin + minDis % DSR_width);
                }else   //视差为y方向
                {
                    disp[x] = -(match_para.disp_ymin + minDis / DSR_width);
                }
                
            }
        }
    }
}

//简单计算SAD匹配
Mat SAD::computerSAD(const Mat &ref, const Mat &tgt)
{
    int Height = ref.rows;
    int Width = ref.cols;
    int winSize = match_para.win_size;
    int half_size = winSize / 2;

    Mat Kernel_L(Size(winSize, winSize), CV_8U, Scalar::all(0));
    Mat Kernel_R(Size(winSize, winSize), CV_8U, Scalar::all(0));
    Mat Disparity(Height, Width, CV_8U, Scalar(0)); //视差图

    for (int i = 0; i<Width - winSize; i++)  //左图从DSR开始遍历
    {
        for (int j = 0; j<Height - winSize; j++)
        {
            Kernel_L = ref(Rect(i, j, winSize, winSize));
            Mat MM(1, DSR, CV_32F, Scalar(0)); //
            for (int k = 0; k<DSR; k++)
            {
                int x = i - k;
                if (x >= 0)
                {
                    Kernel_R = tgt(Rect(x, j, winSize, winSize));

                    Mat Dif;
                    absdiff(Kernel_L, Kernel_R, Dif);//
                    Scalar ADD = sum(Dif);
                    float a = ADD[0];
                    MM.at<float>(k) = a;
                }
            }

            Point minLoc;
            minMaxLoc(MM, NULL, NULL, &minLoc, NULL);
            int loc = minLoc.x;
            Disparity.at<char>(j, i) = loc * 16;
        }
    }
    return Disparity;
}