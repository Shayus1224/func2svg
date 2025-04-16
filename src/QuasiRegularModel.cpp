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
double getPixelValue(const int& i, const int& opt, const int& q, const double& x, const double& y) {
    double pi = 3.1415926f;
    switch (opt) {
        case 0:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));
        case 1:
            return cos(i * x * cos(2 * pi * i / q) + i * y * sin(2 * pi * i / q));//???? i ??????
        case 2:
            return i * cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));//???? i ??????
        case 3:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q) + i);//???? i ???????
        case 4:
            return 2 * cos(x * cos(2 * pi * i / q) +
                           y * sin(2 * pi * i / q));//????????????  a = 2??b = c =1??d =e= 0
        case 5:
            return cos(x * cos(2 * pi * i / q) +
                       y * sin(2 * pi * i / q + 0.5));//????????????  a = 1??b = c =1??d =0.5??e= 0
        case 6:
            return cos(3.5 * x * cos(2 * pi * i / q) +
                       1.5 * y * sin(2 * pi * i / q));//????????????? a = 1?? b =3.5??c = 1.5?? e = 0
        case 7:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   1;//????????????? a = 1?? b =c = 1?? e = 1
        case 8:
            return cos(3.5 * x * cos(2 * pi * i / q) + 1.5 * y * sin(2 * pi * i / q) + 0.5) +
                   1;//????????????? a = 1?? b =3.5??c = 1.5??d =0.5??e = 1
        case 9:
            return abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)));//???????????????????任
        case 10:
            return cos(x * abs(cos(2 * pi * i / q)) + y * abs(sin(2 * pi * i / q)));//???????????????????任
        case 11:
            return abs(cos(x * abs(cos(2 * pi * i / q)) + y * abs(sin(2 * pi * i / q))));//???????????????????任
        case 12:
            return abs(i * cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)));//???????仯???????????任
        case 13:
            return sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));//???????????????????任???
        case 14:
            return tan(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));//???????????????????任????
        case 15:
            return 1 / cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));//???????????????????任????
        case 16:
            return cos(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)));//???????????????????任???
        case 17:
            return sin(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)));//???????????????????任?塣
        case 18:
            return tan(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)));//???????????????????任????
        case 19:
            return cos(abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//???????????????????任??.
        case 20:
            return tan(abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//???????????????????任???
        case 21:
            return cos(x * sin(2 * pi * i / q) + y * cos(2 * pi * i / q));//?????????????????任???
        case 22:
            return cos(x / cos(2 * pi * i / q) + y / sin(2 * pi * i / q));//?????????????????任????
        case 23:
            return cos(sin(x * cos(2 * pi * i / q)) + cos(y * sin(2 * pi * i / q)));//?????????????????任????
        case 24:
            return cos(
                    cos(x) * cos(2 * pi * i / q) + (sin(y) * sin(2 * pi * i / q)));//?????????????????????任? ??
        case 25:
            return cos(tan(x) * cos(2 * pi * i / q) +
                       (1 / tan(y) * sin(2 * pi * i / q)));//?????????????????????任?? ??
        case 26:
            return abs(
                    cos(sin(x * cos(2 * pi * i / q)) + cos(y * sin(2 * pi * i / q))));//???????????????????任? ??
        case 27:
            return cos(abs(sin(x * cos(2 * pi * i / q)) +
                           cos(y * sin(2 * pi * i / q))));//???????????????????任?? ??
        case 28:
            return tan(cos(x) * cos(2 * pi * i / q) + sin(y) * sin(2 * pi * i / q));//???????????????????任?? ??
        case 29:
            return abs(cos(x) * cos(2 * pi * i / q) + sin(y) * sin(2 * pi * i / q));//???????????????????任?? ??
        case 30:
            return cos(tan(cos(x) * cos(2 * pi * i / q)) +
                       tan(sin(y) * sin(2 * pi * i / q)));//???????????????????任?? ??
        case 31:
            return pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 2);//?????????????????任???
        case 32:
            return pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//?????????????????任????
        case 33:
            return pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 15);//?????????????????任????
        case 34:
            return pow(abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))), 0.5);//?????????????????任???
        case 35:
            return pow(tan(pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 2)),
                       3);//?????????????????任??(1)??
        case 36:
            return pow(tan(pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3)),
                       3);//?????????????????任??(2)??
        case 37:
            return pow(abs(pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3)),
                       0.5);//?????????????????任??(3)??
        case 38:
            return cos(
                    x * pow(cos(2 * pi * i / q), 2) + y * pow(sin(2 * pi * i / q), 2));//?????????????????任???
        case 39:
            return cos(
                    x * pow(cos(2 * pi * i / q), 3) + y * pow(sin(2 * pi * i / q), 3));//?????????????????任????
        case 40:
            return cos(x * pow(cos(2 * pi * i / q), 15) +
                       y * pow(sin(2 * pi * i / q), 15));//?????????????????任????
        case 41:
            return cos(x * pow(abs(cos(2 * pi * i / q)), 0.5) +
                       y * pow(abs(sin(2 * pi * i / q)), 0.5));//?????????????????任???
        case 42:
            return cos(
                    pow(abs(x), 0.5) * cos(2 * pi * i / q) +
                    pow(abs(y), 0.5) * sin(2 * pi * i / q));//???????????????????任???
        case 43:
            return cos(
                    pow(abs(x), 0.75) * cos(2 * pi * i / q) +
                    pow(abs(y), 0.75) * sin(2 * pi * i / q));//???????????????????任????
        case 44:
            return cos(
                    pow(abs(x), 1.25) * cos(2 * pi * i / q) +
                    pow(abs(y), 1.25) * sin(2 * pi * i / q));//???????????????????任????
        case 45:
            return cos(sin(pow(abs(x), 0.5)) * cos(2 * pi * i / q) +
                       cos(pow(abs(y), 0.5)) * sin(2 * pi * i / q));//???????????????????任???
        case 46:
            return pow(cos(x * pow(cos(2 * pi * i / q), 2) + y * pow(sin(2 * pi * i / q), 2)),
                       3);//?????????????????任???
        case 47:
            return pow(cos(pow(abs(x), 0.5) * cos(2 * pi * i / q) + pow(abs(y), 0.5) * sin(2 * pi * i / q)),
                       2);//?????????????????任????
        case 48:
            return cos(pow(abs(x), 0.5) * pow(cos(2 * pi * i / q), 2) +
                       pow(abs(y), 0.5) * pow(sin(2 * pi * i / q), 2));//?????????????????任????
        case 49:
            return pow(cos(pow(abs(x), 0.5) * pow(cos(2 * pi * i / q), 2) +
                           pow(abs(y), 0.5) * pow(sin(2 * pi * i / q), 2)), 3);//?????????????????任???
        case 50:
            return cos(sin(pow(abs(x), 0.5)) * pow(cos(2 * pi * i / q), 2) +
                       cos(pow(abs(y), 0.5)) * pow(sin(2 * pi * i / q), 2));//?????????????????任?塣
        case 51:
            return exp(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)));//???????????????????任?
        case 52:
            return exp(abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//???????????????????任??
        case 53:
            return cos(x * exp(cos(2 * pi * i / q)) + y * exp(sin(2 * pi * i / q)));//???????????????????任?
        case 54:
            return cos(x * exp(pow(abs(cos(2 * pi * i / q)), 0.5)) +
                       y * exp(pow(abs(sin(2 * pi * i / q)), 0.5)));//???????????????????任??
        case 55:
            return exp(cos(x * exp(cos(2 * pi * i / q)) + y * exp(sin(2 * pi * i / q))));//???????????????????任?
        case 56:
            return exp(pow(abs(cos(x * exp(cos(2 * pi * i / q)) + y * exp(sin(2 * pi * i / q)))),
                           0.5));//???????????????????任??
        case 57:
            return log(abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//???????????????????任?
        case 58:
            return log(
                    abs(tan(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)))));//???????????????????任??
        case 59:
            return cos(x * log(abs(cos(2 * pi * i / q))) +
                       y * log(abs(sin(2 * pi * i / q))));//???????????????????任?
        case 60:
            return cos(log(abs(x)) * cos(2 * pi * i / q) +
                       log(abs(y)) * sin(2 * pi * i / q));//???????????????????任??
        case 61:
            return log(abs(cos(log(abs(x)) * cos(2 * pi * i / q) +
                               log(abs(y)) * sin(2 * pi * i / q))));//???????????????????任?
        case 62:
            return log(abs(cos(x * log(abs(cos(2 * pi * i / q))) +
                               y * log(abs(sin(2 * pi * i / q))))));//???????????????????任??
        case 63:
            return cos(log(abs(x)) * log(abs(cos(2 * pi * i / q))) +
                       log(abs(y)) * log(abs(sin(2 * pi * i / q))));//???????????????????任??
        case 64:
            return log(abs(cos(log(abs(x)) * log(abs(cos(2 * pi * i / q))) +
                               log(abs(y)) * log(abs(sin(2 * pi * i / q))))));//???????????????????任??
        case 65:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   (sin(x) + cos(y)) / 5;//??????????????????
        case 66:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   cos((sin(x) + y)) / 5;//???????????????????
        case 67:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) + (x + y) / 25;//???????????????
        case 68:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   (pow(x, 2) + pow(y, 2)) / 200;//????????????????
        case 69:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   pow(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 2);//????????????????
        case 70:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   pow(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 2) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//?????????????????
        case 71:
            return tan(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 2);//?????????????????
        case 72:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 100) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 101);//?????????????????
        case 73:
            return 5 * pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 5);//?????????????????
        case 74:
            return 5 * pow(abs(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))), 5);//?????????????????
        case 75:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   exp(cos(x) + sin(y));//??????????????????
        case 76:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   log(x * x + y * y) / 5;//??????????????????
        case 77:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                   log(abs((sin(x) + cos(y)))) / 5;//???????????????????
        case 78:
            return abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   pow(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//????????????????任???????????
        case 79:
            return cos(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   pow(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//???????????????任???????????
        case 80:
            return pow(cos(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))), 3) +
                   sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q));//???????????????任????????????
        case 81:
            return tan(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//???????????????任????????????
        case 82:
            return exp(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//???????????????任???????????
        case 83:
            return exp(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//???????????????任????????????
        case 84:
            return log(abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)))) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//???????????????任???????????
        case 85:
            return log(abs(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)))) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3);//???????????????任????????????
        case 86:
            return log(abs(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)))) +
                   pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)),
                       100);//???????????????任????????????
        case 87:
            return abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   log(abs(cos(
                           x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//????????????????任?????????????
        case 88:
            return sin(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   log(abs(sin(
                           x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//???????????????任?????????????
        case 89:
            return exp(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) +
                   log(abs(cos(
                           x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//???????????????任?????????????
        case 90:
            return sqrt(exp(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)))) +
                   log(abs(cos(
                           x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//???????????????任??????????????
        case 91:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) * cos(x) *
                   cos(y);//??????????????????任?
        case 92:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) * sin(x) *
                   cos(y);//??????????????????任??
        case 93:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) * sin(x - y) *
                   cos(x + y);//??????????????????任??
        case 94:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) * exp(sin(x)) *
                   exp(cos(y));//???????????????任?
        case 95:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) * exp(sin(x)) * exp(cos(y)) *
                   exp(sin(y)) *
                   exp(cos(x));//???????????????任??
        case 96:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) * exp(sin(x)) * exp(cos(y)) *
                   log(abs(sin(x) * cos(y)));//???????????????任??
        case 97:
            return exp(cos(x * sin(cos(2 * pi * i / q)) + y * sin(sin(2 * pi * i / q)))) +
                   sqrt(abs(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))));//??????????????任?
        case 98:
            return pow(tan(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) *
                           pow(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 3)),
                       4);//??????????????任??
        case 99:
            return exp(pow(cos(x * cos(cos(2 * pi * i / q)) + y * sin(sin(2 * pi * i / q))),
                           q));//??????????????????????任????
        case 100:
            return exp(pow(cos(x * cos(cos(2 * pi * i / q)) + y * sin(sin(2 * pi * i / q))),
                           i));//??????????????????????任????
        case 101:
            return cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) * cos(2 * pi * i / q) *
                   sin(2 * pi * i / q) *
                   (-1);//????????????????任
        case 102:
            return (-(exp(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) *
                      cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) +
                      2 * sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) * cos(2 * pi * i / q) *
                    sin(2 * pi * i / q));//???????????????????????任
        case 103:
            return (-sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) *
                    ((2 * x * pi * i / pow(q, 2)) * sin(2 * pi * i / q) -
                     2 * y * pi * i / pow(q, 2) * cos(2 * pi * i / q)));//????????q??任
        case 104:
            return exp(sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q))) *
                   cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) *
                   ((2 * x * pi * i / pow(q, 2)) * sin(2 * pi * i / q) -
                    2 * y * pi * i / pow(q, 2) * cos(2 * pi * i / q)) -
                   3 * pow(cos(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)), 2) *
                   ((2 * x * pi * i / pow(q, 2)) * sin(2 * pi * i / q) -
                    2 * y * pi * i / pow(q, 2) * cos(2 * pi * i / q));    //?????????????????q??任
        case 105:
            return -sin(x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q)) *
                   (2 * x * pi / q * sin(2 * pi * i / q) - 2 * y * pi / q * cos(2 * pi * i / q));//????????i??任
        case 106:
            return exp(sin(x * sin(2 * pi * i / q) + y * cos(2 * pi * i / q))) *
                   cos(x * sin(2 * pi * i / q) + y * cos(2 * pi * i / q)) *
                   ((2 * x * pi / q) * cos(2 * pi * i / q) - 2 * y * pi / q * sin(2 * pi * i / q)) +
                   4 * pow(cos(x * sin(2 * pi * i / q) + y * cos(2 * pi * i / q)), 3) *
                   (2 * x * pi / q * cos(2 * pi * i / q) -
                    2 * y * pi / q * sin(2 * pi * i / q));//?????????????????i??任
        case 107: {
            float t = 2 * pi * i / q;
            float Cost = 1 - t * t / 2 + pow(t, 4) / 24 - pow(t, 6) / 720 + pow(t, 8) / 40320 -
                         pow(t, 10) / 3628800 +
                         pow(t, 12) / 3628800 / 11 / 12 - pow(t, 14) / 3628800 / 11 / 12 / 13 / 14;
            float Sint = t - t * t * t / 6 + pow(t, 5) / 120 - pow(t, 7) / 5040 + pow(t, 9) / 362880 -
                         pow(t, 11) / 3628800 / 11 + pow(t, 13) / 3628800 / 11 / 12 / 13 -
                         pow(t, 15) / 3628800 / 11 / 12 / 13 / 14 / 15;
            return cos(x * Cost + y * Sint);//???????
        }
        case 108: {
            float t = 2 * pi * i / q;
            float Cost = 1 - t * t / 2 + pow(t, 4) / 24 - pow(t, 6) / 720 + pow(t, 8) / 40320 -
                         pow(t, 10) / 3628800 +
                         pow(t, 12) / 3628800 / 11 / 12;
            float Sint = t - t * t * t / 6 + pow(t, 5) / 120 - pow(t, 7) / 5040 + pow(t, 9) / 362880 -
                         pow(t, 11) / 3628800 / 11 + pow(t, 13) / 3628800 / 11 / 12 / 13;
            return cos(x * Cost + y * Sint);//????????
        }
        case 109: {
            float t = 2 * pi * i / q;
            float Cost = 1 - t * t / 2 + pow(t, 4) / 24 - pow(t, 6) / 720 + pow(t, 8) / 40320;
            float Sint = t - t * t * t / 6 + pow(t, 5) / 120 - pow(t, 7) / 5040 + pow(t, 9) / 362880;
            return cos(x * Cost + y * Sint);//????????
        }
        case 110: {
            float t = x * cos(2 * pi * i / q) + y * sin(2 * pi * i / q);
            float Cost = 1 - t * t / 2 + pow(t, 4) / 24 - pow(t, 6) / 720 + pow(t, 8) / 40320 -
                         pow(t, 10) / 3628800 +
                         pow(t, 12) / 3628800 / 11 / 12 - pow(t, 14) / 3628800 / 11 / 12 / 13 / 14;
            return Cost;//????????
        }
        case 111: {
            float t = 2 * pi * i / q;
            float Cost = 1 - t * t / 2 + pow(t, 4) / 24 - pow(t, 6) / 710 + pow(t, 8) / 40320 -
                         pow(t, 10) / 3628800 +
                         pow(t, 12) / 3628800 / 11 / 12 - pow(t, 14) / 3628800 / 11 / 12 / 13 / 14;
            float Sint = t - t * t * t / 6 + pow(t, 5) / 121 - pow(t, 7) / 5040 + pow(t, 9) / 362880 -
                         pow(t, 11) / 3628800 / 11 + pow(t, 13) / 3628800 / 11 / 12 / 13 -
                         pow(t, 15) / 3628800 / 11 / 12 / 13 / 14 / 15;
            return cos(x * Cost + y * Sint);//????????
        }
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
        colorNum = values[values.size()-1];
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
	return std::to_string(imageOption)+ "_" + std::to_string(iterationTime)+ "_" + std::to_string(colorChoice);
}