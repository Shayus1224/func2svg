#include "QuasiRegularModel.h"

// default generate a picture
QuasiRegularModel::QuasiRegularModel()
{
	this->imageOption = 0;
	this->iterationTime = 5;
	this->colorChoice = 0;
}

// costume generate a picture
QuasiRegularModel::QuasiRegularModel(const int& imageOption,const int& iterationTime, const int& colorChoice)
{
	this->imageOption = imageOption;
	this->iterationTime = iterationTime;
	this->colorChoice = colorChoice;
}

double getPixelValue(const int& i, const int& opt, const int& q, const double& x, const double& y)
{
	double pi = 3.1415926f;
	switch (opt) {
	case 0:
		return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));
	case 1:
		return cos(i * x * cos(2 * pi * i / q) + i * y * sin(2 * pi * i / q));
	case 2:
		return i * cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));
	case 3:
		return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q) + i);
	case 4:
		return 2 * cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));
	case 5:
		return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q + 0.5));
	default:
		return 0;
	}
}

int findColorAndPosition(double x, const std::vector<double>& values)
{
	// 遍历 values 找到 x 所处的区间
	for (size_t i = 0; i < values.size(); ++i) {
		if (x < values[i]) {
			// 如果 x 小于当前分界点，返回区间索引
			return i ;
		}
	}
	// 如果 x 大于等于 values 的最大值，返回区间索引
	return static_cast<int>(values.size());
}


PointInfo QuasiRegularModel::getPointInfo(const double& x, const double& y)
{
	PointInfo target = {};

	for (int i = 0; i < this->iterationTime; ++i) {
		target.value += getPixelValue(i, this->imageOption, this->iterationTime, x, y);
	}
	target.colorNum = findColorAndPosition(target.value, this->valueScheme[this->colorChoice]);
	target.Color = this->Scheme[this->colorChoice][target.colorNum];

	return target;
}

int QuasiRegularModel::getLevels()
{
	return this->Scheme[colorChoice].size();
}

std::string QuasiRegularModel::getImageName()
{
	return std::to_string(imageOption) + std::to_string(iterationTime) + std::to_string(colorChoice);
}