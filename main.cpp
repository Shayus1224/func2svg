#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include "alglib/optimization.h"
#include <iostream>

#include "QuasiRegularModel.h"
#include "ImageProcess.h"
#include "ImageRestore.h"

int main()
{

    size_t t = 5;
    std::cout << "start process  " << t << " choice" << std::endl;
    QuasiRegularModel qrm(t, 5, 0);
    ImageProcess ip(qrm);

    ip.to_image();
    ip.to_linesImage();

    //qrm.to_pointImage();
    // qrm.to_linesImage();
    // qrm.to_svg();


    return 0;
}