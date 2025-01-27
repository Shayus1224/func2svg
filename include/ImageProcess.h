#pragma once
#include "QuasiRegularModel.h"
#include "topology.h"

class ImageProcess {
private:
	QuasiRegularModel& image;
    
    // functions
    std::vector<Point> marching_squares_points(const std::vector<std::vector<double>>& grid, double level);
    Point interpolate(double x1, double y1, double x2, double y2, double v1, double v2, double level);
    int get_case(double v1, double v2, double v3, double v4, double level);
    std::vector<std::vector<Point>> connectPoints(const std::vector<Point>& points);
    double distance(const Point& p1, const Point& p2);
    void tryConnectCurves(std::vector<std::vector<Point>>& curves);
    void getTopology(const std::vector<std::vector<std::vector<Point>>>& curves, const std::vector<cv::Mat>& images);
    double findClosestPoints(const std::vector<Point>& set1, const std::vector<Point>& set2);

public:

    ImageProcess(QuasiRegularModel& model) : image(model) { };
    // generate original image
    void to_image();
    void to_linesImage();

};

