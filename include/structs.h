#pragma once
#include <Eigen/Dense>
#include <iostream>

struct Point {
    double x;
    double y;
    double h;
    bool visited;
    Point()
    {
        this->x = 0;
        this->y = 0;
        this->h = 0;
        this->visited = 0;
    }

    Point(const double& x, const double& y)
    {
        this->x = x;
        this->y = y;
    }
    bool operator==(const Point& other) const
    {
        return std::abs(x - other.x) < 1e-6 && std::abs(y - other.y) < 1e-6;
    }

    void operator=(const Point& other)
    {
        this->x = other.x;
        this->y = other.y;
        this->h = other.h;
    }
    bool operator!=(const Point& other) const
    {
        return other.x != x || other.y != y;
    }
    bool operator<(const Point& other) const
    {
        if (x == other.x) return y < other.y;
        return x < other.x;
    }
};

struct PointInfo {
    double value;
    int colorNum;
    cv::Scalar Color;
};

struct Segment {
    Point p1, p2;
};