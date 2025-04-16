#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include "alglib/optimization.h"
#include <iostream>

#include "QuasiRegularModel.h"
#include "ImageProcess.h"
#include "ImageRestore.h"

int main()
{
for(int t=0;t<111;t++) {
//    int t = 3;
    std::cout << "start process  " << t << " choice" << std::endl;
    QuasiRegularModel qrm(t, 5, 0);
    ImageProcess ip(qrm);
    ip.to_train_image();
//    ip.to_image();
//    ip.to_linesImage();
    //ip.to_linesImage_new();
}

    return 0;
}