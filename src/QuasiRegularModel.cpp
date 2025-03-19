#include "QuasiRegularModel.h"

/**
 * @brief 准规则斑图初始化 0，5，0
 *
 * 根据规定的第0种图像函数，迭代5次，第一种配色的方法生成准规则斑图
 *
 * @param imageOption 图像生成函数选择
 * @param iterationTime 函数迭代次数
 * @param colorChoice 配色选择
 */
QuasiRegularModel::QuasiRegularModel() : imageOption(0), iterationTime(5), colorChoice(0) {}

/**
 * @brief 准规则斑图初始化
 *
 * 根据规定的第imageOption种图像函数，迭代iterationTime次，第colorChoice配色的方法生成准规则斑图
 *
 * @param imageOption 图像生成函数选择
 * @param iterationTime 函数迭代次数
 * @param colorChoice 配色选择
 */
QuasiRegularModel::QuasiRegularModel(const int& imageOption,const int& iterationTime, const int& colorChoice)
{
	this->imageOption = imageOption;
	this->iterationTime = iterationTime;
	this->colorChoice = colorChoice;
}

/**
 * @brief 计算图像中某个点第i次迭代的高度值
 *
 * 根据给定的坐标 (x, y) 计算该点的高度值
 *
 * @param i 第i次迭代
 * @param opt 图像函数的选择
 * @param q 总计迭代次数
 * @param x 点的横坐标
 * @param y 点的纵坐标
 * @return double 某点在第i次迭代的高度值
 */
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

/**
 * @brief 根据给定的高度值 x，在一组颜色值 values 中找到对应的颜色位置
 *
 * 根根据给定的高度值 x，在一组颜色值 values 中找到对应的颜色位置
 *
 * @param x 高度值
 * @param values 配色
 * @return int 配色索引
 */
int findColorAndPosition(double x, const std::vector<double>& values)
{
    if (values.empty()) {
        throw std::invalid_argument("values is empty");
    }

    int colorNum = -1;
    for (size_t i = 0; i < values.size(); ++i) {
        if (x <= values[i]) {
            colorNum = static_cast<int>(i);
            break;
        }
    }

    // 检查是否找到有效的颜色
    if (colorNum == -1) {
        throw std::out_of_range("value is greater than all elements in values");
    }

    return colorNum;
}

/**
 * @brief 计算指定点的颜色信息
 *
 * 根据给定的坐标 (x, y) 计算该点的强度值，并根据预设的颜色方案确定颜色。
 *
 * @param x 点的横坐标
 * @param y 点的纵坐标
 * @return std::tuple<int, double> 包含颜色索引和该点的强度值
 */
PointInfo QuasiRegularModel::getPointInfo(const double& x, const double& y)
{
	PointInfo target = {};

	for (int i = 0; i < this->iterationTime; ++i) {
		target.value += getPixelValue(i, this->imageOption, this->iterationTime, x, y);
	}
    try {
        target.colorNum = findColorAndPosition(target.value, this->valueScheme[this->colorChoice]);
    } catch (const std::invalid_argument& e) {
        // 处理 valueScheme 为空的情况
        throw std::runtime_error("Value scheme is empty or invalid");
    }
    target.Color = this->Scheme[this->colorChoice][target.colorNum];
	return target;
}

/**
 * @brief 返回图像方案中颜色级别的数量
 *
 * 返回图像方案中颜色级别的数量
 *
 * @return int 返回图像方案中颜色级别的数量
 */
[[maybe_unused]] int QuasiRegularModel::getLevels()
{
	return (int)this->Scheme[colorChoice].size();
}

/**
 * @brief 返回图像名称（根据图像选项、迭代次数、颜色选择组合生成）
 *
 * 返返回图像名称（根据图像选项、迭代次数、颜色选择组合生成）
 *
 * @return string 返回图像名称（根据图像选项、迭代次数、颜色选择组合生成）
 */
std::string QuasiRegularModel::getImageName() const
{
	return std::to_string(imageOption) + std::to_string(iterationTime) + std::to_string(colorChoice);
}