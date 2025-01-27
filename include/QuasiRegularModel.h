#pragma once
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

#include "structs.h"


class QuasiRegularModel {
public:
    QuasiRegularModel();
    QuasiRegularModel(const int& imageOption,const int& iterationTime, const int& colorChoice);

    PointInfo getPointInfo(const double& x, const double& y);
    int getLevels();
    std::string getImageName();

    const int GRID_SIZE = 512;
    const double CELL_SIZE = 40.0 / GRID_SIZE;
    const int imageSize = 512;
    const double offset = -20.0;
    const double cellSize = this->CELL_SIZE;
    int imageOption;
    int iterationTime;
    int colorChoice;
    const std::vector<std::vector<cv::Scalar>> Scheme = {
        {
            cv::Scalar(160, 10, 180), cv::Scalar(10, 100, 200), cv::Scalar(50, 50, 150), cv::Scalar(250, 0, 190),
            cv::Scalar(0, 180, 100), cv::Scalar(250, 0, 250), cv::Scalar(250, 0, 0), cv::Scalar(255, 0, 0),
            cv::Scalar(90, 90, 90), cv::Scalar(50, 150, 220), cv::Scalar{0, 50, 0}, cv::Scalar(0, 0, 70),
            cv::Scalar(0, 0, 150), cv::Scalar(0, 180, 180), cv::Scalar(0, 100, 100), cv::Scalar(160, 0, 160),
            cv::Scalar(250, 5, 250), cv::Scalar(200, 20, 180), cv::Scalar(160, 10, 180)
        }
    };

    const std::vector<std::vector<double>> valueScheme = {
        {
            -9, -5, -4, -3, -2, -1, -0.5,
            -0.1, 0, 0.1, 0.2, 0.3, 0.5, 1,
            2.2, 3, 5, 7, 9, 11
        }
    };
};
