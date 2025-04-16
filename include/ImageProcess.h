#pragma once
#include "QuasiRegularModel.h"
#include "topology.h"

class ImageProcess {
private:
	QuasiRegularModel& image;
    
    // version 1 : functions
    static std::vector<Point> marching_squares_points(const std::vector<std::vector<double>>& grid, double level);
    static Point interpolate(double x1, double y1, double x2, double y2, double v1, double v2, double level);
    static int get_case(double v1, double v2, double v3, double v4, double level);
    static std::vector<std::vector<Point>> connectPoints(const std::vector<Point>& points);
    static double distance(const Point& p1, const Point& p2);
    static void tryConnectCurves(std::vector<std::vector<Point>>& curves);
    void getTopology(const std::vector<std::vector<std::vector<Point>>>& curves, const std::vector<cv::Mat>& images);
    static double findClosestPoints(const std::vector<Point>& set1, const std::vector<Point>& set2);
    void sortTopology(std::vector<std::vector<std::vector<cv::Point>>> contourss,
                      std::vector<std::vector<cv::Vec4i>> hierarchys);



public:

    explicit ImageProcess(QuasiRegularModel& model) : image(model) { };
    // generate original image
    void to_image();
    void to_linesImage();

    void to_train_image();
};

