#include <stack>
#include "ImageProcess.h"
#include <alglib/interpolation.h>
#include <graphviz/gvc.h>
/**
 * @brief 用于生成四个点的状态值
 *
 * 通过将每个点与 level 值进行比较，计算出一个二进制状态（四个点的状态用一个整数表示，按位进行“或”操作）。
 *
 * @param v1 v1点高度值
 * @param v2 v2点高度值
 * @param v3 v3点高度值
 * @param v4 v4点高度值
 * @param level 比较高度值
 * @return int 该square对应的case值
 */
int ImageProcess::get_case(double v1, double v2, double v3, double v4, double level)
{
	return (v1 >= level) | ((v2 >= level) << 1) | ((v3 >= level) << 2) | ((v4 >= level) << 3);
}

/**
 * @brief 计算两个点之间的插值点
 *
 * 计算两个点之间的插值点
 *
 * @param x1 point1 x
 * @param y1 point1 y
 * @param x2 point2 x
 * @param y2 point2 y
 * @param v1 point1 高度值
 * @param v2 point2 高度值
 * @param level 比较高度值
 * @return Point 拟合插值点
 */
Point ImageProcess::interpolate(double x1, double y1, double x2, double y2, double v1, double v2, double level)
{
	double t = (level - v1) / (v2 - v1);
	return { x1 + t * (x2 - x1), y1 + t * (y2 - y1) };
}

/**
 * @brief 计算两个点之间的欧式距
 *
 * 计算两个点之间的欧式距
 *
 * @param p1 point1
 * @param p2 point2
 * @return double 距离值
 */
double ImageProcess::distance(const Point& p1, const Point& p2)
{
	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief 寻找最近点距离
 *
 * 通过遍历set中的每一个点，找出最小距离
 *
 * @param set1
 * @param set2
 * @return double 最近点距离
 */
double ImageProcess::findClosestPoints(const std::vector<Point>& set1, const std::vector<Point>& set2)
{
	double minDistance = std::numeric_limits<double>::max();
	// std::pair<Point, func2_Point> closestPair;

	for (const auto& point1 : set1) {
		for (const auto& point2 : set2) {
			double dis = distance(point1, point2);
			if (dis < minDistance) {
				minDistance = dis;
				//                closestPair = std::make_pair(point1, point2);
			}
		}
	}

	return minDistance;
}

/**
 * @brief 尝试将曲线连接起来
 *
 * 通过判断两条曲线之间最近点的距离是否足够接近来合并曲线
 *
 * @param curves
 */
void ImageProcess::tryConnectCurves(std::vector<std::vector<Point>>& curves)
{
	bool merged;
	do {
		merged = false;
		for (int i = 0; i < curves.size(); ++i) {
			for (int j = i + 1; j < curves.size(); ++j) {
				if (distance(curves[i].back(), curves[j].front()) < 3) {
					curves[i].insert(curves[i].end(), curves[j].begin(), curves[j].end());
					curves.erase(curves.begin() + j);
					merged = true;
					break;
				} else if (distance(curves[i].front(), curves[j].back()) < 3) {
					curves[j].insert(curves[j].end(), curves[i].begin(), curves[i].end());
					curves.erase(curves.begin() + i);
					merged = true;
					break;
				} else if (distance(curves[i].back(), curves[j].back()) < 3) {
					std::reverse(curves[j].begin(), curves[j].end());
					curves[i].insert(curves[i].end(), curves[j].begin(), curves[j].end());
					curves.erase(curves.begin() + j);
					merged = true;
					break;
				} else if (distance(curves[i].front(), curves[j].front()) < 3) {
					std::reverse(curves[i].begin(), curves[i].end());
					curves[i].insert(curves[i].end(), curves[j].begin(), curves[j].end());
					curves.erase(curves.begin() + j);

					merged = true;
					break;
				}
			}
			if (merged) break;
		}
	} while (merged);
}

/**
 * @brief 用marching square方法提出边界曲线
 *
 * 该函数遍历整个网格，计算每个格子四个顶点的状态，并根据状态决定该格子中的哪些边应该被插值。
 * 插值使用的是 interpolate 函数来计算边的中点。
 * 最后返回所有边界点。
 *
 * @param grid
 * @param level
 * @return
 */
std::vector<Point> ImageProcess::marching_squares_points(const std::vector<std::vector<double>>& grid, double level)
{
	std::vector<Point> res;
	auto rows = grid.size();
    auto cols = grid[0].size();

	for (auto i = 0; i < rows - 1; ++i) {
		for (auto j = 0; j < cols - 1; ++j) {
			// 获取当前格子四个顶点的值和坐标
			double v0 = grid[i][j], v1 = grid[i][j + 1];
			double v2 = grid[i + 1][j + 1], v3 = grid[i + 1][j];
			Point p0 = { j, i }, p1 = { j + 1, i }, p2 = { j + 1, i + 1 }, p3 = { j, i + 1 };

			// 计算边中点
			Point m0 = interpolate(p0.x, p0.y, p1.x, p1.y, v0, v1, level);
			Point m1 = interpolate(p1.x, p1.y, p2.x, p2.y, v1, v2, level);
			Point m2 = interpolate(p2.x, p2.y, p3.x, p3.y, v2, v3, level);
			Point m3 = interpolate(p3.x, p3.y, p0.x, p0.y, v3, v0, level);

			// 根据 case 确定边界点
			switch (get_case(v0, v1, v2, v3, level)) {
			case 1:
			case 14: res.insert(res.end(), { m3, m0 }); break;
			case 2:
			case 13: res.insert(res.end(), { m0, m1 }); break;
			case 3:
			case 12: res.insert(res.end(), { m3, m1 }); break;
			case 4:
			case 11: res.insert(res.end(), { m1, m2 }); break;
			case 5: res.insert(res.end(), { m0, m1, m2, m3 }); break;
			case 6:
			case 9: res.insert(res.end(), { m0, m2 }); break;
			case 7:
			case 8: res.insert(res.end(), { m2, m3 }); break;
			case 10: res.insert(res.end(), { m0, m3, m1, m2 }); break;
			}
		}
	}

	return res;
}

/**
 * @brief 将提取出的边界点进行连接，生成完整的曲线
 * @param points
 * @return std::vector<std::vector<Point>>
 */
std::vector<std::vector<Point>> ImageProcess::connectPoints(const std::vector<Point>& points)
{
	auto n = points.size();
	std::vector<bool> visited(n, false);
	std::vector<std::vector<Point>> curves;

	auto compareDist = [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
		return a.first > b.first;
    };

	for (int i = 0; i < n; ++i) {
		if (!visited[i]) {
			std::vector<Point>curve;
			std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, decltype(compareDist)> pq(compareDist);
			pq.emplace( 0, i );

			while (!pq.empty()) {
				auto dist = pq.top();
				auto current = dist.second;
				pq.pop();

				if (visited[current]) continue;
				visited[current] = true;
				curve.push_back(points[current]);

				for (int j = 0; j < n; ++j) {
					if (!visited[j]) {
						double dis = distance(points[current], points[j]);
						if (dis < 3) {
							pq.emplace( dis, j );
						}
					}
				}
			}

			curves.push_back(curve);
		}
	}

	tryConnectCurves(curves);

	return curves;
}

/**
 * @brief 生成图像的拓扑结构
 *
 * 遍历每一层图像，提取出颜色点，并将这些点进行连接，生成节点。
 * 通过分析相邻层之间的关系，创建节点之间的边，并将图像拓扑结构连接起来。
 * 拓扑最终会经过整理，去除冗余边，形成最终的拓扑结构。
 *
 * @param curves
 * @param images
 */
void ImageProcess::getTopology(const std::vector<std::vector<std::vector<Point>>>& curves, const std::vector<cv::Mat>& images)
{
	std::cout << "开始整理拓扑" << std::endl;
	Topology topology;
	// auto colorScheme = getColorValues();
	std::vector<std::vector<std::vector<Point>>> colorDots_s;
	// 初始化所有拓扑节点
	for (int t = 0; t < images.size() - 1; ++t) {
		auto image_t = images[t];
		cv::Mat grayImg;
		// 将彩色图像转换为灰度图像
		cv::cvtColor(image_t, grayImg, cv::COLOR_BGR2GRAY);
		// 遍历图像
		std::vector<Point>colorPoint;
		for (int i = 0; i < image_t.rows; ++i) {
			const auto* row_ptr = image_t.ptr<uchar>(i); // 获取第i行的指针
			for (int j = 0; j < image_t.cols; ++j) {
				if (row_ptr[j] == 0) { // 假设条件是像素值等于255
					colorPoint.emplace_back(j, i); // 存储满足条件的点的位置
				}
			}
		}
		auto colorDots = connectPoints(colorPoint);
		colorDots_s.emplace_back(colorDots);
		for (int i = 0; i < colorDots.size(); ++i) {
			int id = ((t + 1) * 10000) + i;
			cv::Scalar color = image.Scheme[image.colorChoice][i];
			topology.addNode(id, color);
		}
		std::cout << t << "层节点建立完成, 现在有" << topology.nodes.size() << std::endl;
	}
	// 整理所有的拓扑，一个vertex + 一个edge
	for (int t = 0; t < images.size() - 1; ++t) {
		auto colorDots = colorDots_s[t];
		auto lines = curves[t];
		for (int i = 0; i < colorDots.size(); ++i) {
			// 1. 找出距离很近的边
			for (int k = 0; k < lines.size(); ++k) {
				auto dis = findClosestPoints(colorDots[i], lines[k]);
				if (dis < 1) {
					int dotId = ((t + 1) * 10000) + i;
					int edgeId = ((t + 1) * 10000) + k;
					// topology.addNode(dotId, getColor(colorScheme[t]));
					topology.addEdge(dotId, 0, edgeId, lines[k]);
				}
			}
			lines = curves[t + 1];
			for (int k = 0; k < lines.size(); ++k) {
				auto dis = findClosestPoints(colorDots[i], lines[k]);
				if (dis < 1) {
					int dotId = ((t + 1) * 10000) + i;
					int edgeId = ((t + 2) * 10000) + k;
					// topology.addNode(dotId, getColor(colorScheme[t]));
					topology.addEdge(dotId, 0, edgeId, lines[k]);
				}
			}
		}
		std::cout << t << "层拓扑建立完成, 现在有" << topology.edges.size() << std::endl;
	}

	// 将散开的拓扑连接起来
	topology.connectTopo();
	topology.removeRedundantEdges();

	int level = 0;
    cv::Mat image_ee(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(255, 255, 255));

    for (auto& l : topology.edges) {
        auto points = l.outerLine;
        std::set<Point> uniquePoints(points.begin(), points.end());
        points.assign(uniquePoints.begin(), uniquePoints.end());

        auto size = points.size();
        std::cout << "第" << level++ << "个有" << size << "个点\n";
        // 将数据点转换为alglib格式
        cv::Mat image_e(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(255, 255, 255));
        std::ofstream outfile("line_" + std::to_string(level) + "point.txt");

        alglib::real_1d_array x_array, y_array, t_array;

        x_array.setlength((int) size);
        y_array.setlength((int) size);
        t_array.setlength((int) size);
        if (outfile.is_open()) {
            for (int i = 0; i < points.size(); ++i) {
                int x = (int) points[i].x;
                int y = (int) points[i].y;

                image_e.at<cv::Vec3b>(y, x)[0] = 0;
                image_e.at<cv::Vec3b>(y, x)[1] = 0;
                image_e.at<cv::Vec3b>(y, x)[2] = 0;
                image_ee.at<cv::Vec3b>(y, x)[0] = 0;
                image_ee.at<cv::Vec3b>(y, x)[1] = 0;
                image_ee.at<cv::Vec3b>(y, x)[2] = 0;
                outfile << x << " , " << y << std::endl;
                x_array[i] = x;
                y_array[i] = y;
                t_array[i] = i;
            }
        } else {
            std::cerr << "无法打开文件进行追加   " << std::to_string(level) << std::endl;
        }
        outfile.close();
        std::string file_name = "line_" + std::to_string(level) + ".jpg";
        cv::imwrite(file_name, image_e);
        try {

            alglib::spline1dinterpolant res1,res2;
            alglib::spline1dbuildcubic(t_array, x_array, res1);
            alglib::spline1dbuildcubic(t_array, y_array, res2);
            l.line1 = res1;
            l.line2 = res2;

        } catch (alglib::ap_error alglib_exception) {

            printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());

        }
    }
    cv::imwrite("all_lines.jpg", image_ee);
//	topology.printTopology();
//	std::string name = image.getImageName();
//	topology.generateDotFile(name);
}

/**
 * 将图像从像素值转换为色块图像
 */
void ImageProcess::to_image() {
    cv::Mat bezierFittedImage(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(255, 255, 255));
    // 获取图像指针

    for (int i = 0; i < image.imageSize; ++i) {
        auto *rowPtr = bezierFittedImage.ptr<cv::Vec3b>(i);
        for (int j = 0; j < image.imageSize; ++j) {
            // 计算像素点的位置和颜色
            double x = image.offset + image.cellSize * i;
            double y = image.offset + image.cellSize * j;
            auto pi = image.getPointInfo(x, y);

            rowPtr[j] = cv::Vec3b((uchar) pi.Color[2], (uchar) pi.Color[1], (uchar) pi.Color[0]);
        }
    }
    std::string n = image.getImageName() + "_originalImage.jpg";
    std::cout << "a original image is generated as " + n << std::endl;
    cv::imwrite(n, bezierFittedImage);

}

void ImageProcess::to_linesImage()
{
	std::vector<std::vector<double>> grid;
	for (int i = 0; i < image.imageSize; ++i) {
		std::vector<double> line;
		for (int j = 0; j < image.imageSize; ++j) {
			double x = image.offset + image.cellSize * i;
			double y = image.offset + image.cellSize * j;
			auto pi = image.getPointInfo(x, y);
			line.emplace_back(pi.value);
		}
		grid.emplace_back(line);
	}

	std::vector<std::vector<Point>> points;
	auto levels = image.valueScheme[image.colorChoice];
	for (double n : levels) {
		auto t = marching_squares_points(grid, n);
		points.push_back(t);
	}

	std::vector<std::vector<std::vector<Point>>> curves;
	for (const auto& func2_Point : points) {
		auto curve = connectPoints(func2_Point);
		curves.emplace_back(curve);
	}
    cv::Mat bezierFittedImage(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0, 0, 0));
    int count = 0;
    for(const auto& curve : curves){

        cv::Mat im(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0, 0, 0));
        for(const auto& c: curve) {
            for (const auto &cc: c) {
                bezierFittedImage.at<cv::Vec3b>(cc.y, cc.x) = cv::Vec3b(255, 255, 255);
                im.at<cv::Vec3b>(cc.y, cc.x) = cv::Vec3b(255, 255, 255);
            }
        }
        // cv::imwrite(std::to_string(levels[count]) + "_test.jpg", bezierFittedImage);
    }
    // cv::imwrite("test.jpg", bezierFittedImage);
    // info: 使用marching square提取出的边界线 curves 无问题
	// 生成辅助图片
	std::vector<cv::Mat> images;
	for (int a = 0; a < curves.size(); a++ ) { // a 表示层数
		double m, n;
		auto vScheme = image.valueScheme[image.colorChoice];  // k
        auto cScheme = image.Scheme[image.colorChoice];// k+1
        // double target_m = colorScheme[a+1];
		if ((a - 1) < 0) { // a = 0 的情况
			m = -std::numeric_limits<double>::max();
			n = vScheme[a];
		} else if (a+1 == vScheme.size()) { // a + 1 = k
			m = vScheme[a];
			n = std::numeric_limits<double>::max();
		} else { // a [1, k-2]
			m = vScheme[a];
			n = vScheme[a + 1];
		}
		// auto color = getColor(colorScheme[a]);
		cv::Mat image_e(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0, 0, 0));

		for (int i = 0; i < image_e.rows; ++i) {
			// 获取指向第 i 行的指针
			auto* rowPtr = image_e.ptr<cv::Vec3b>(i);
			for (int j = 0; j < image_e.cols; ++j) {
				// 获取(i,j)像素值
				cv::Vec3b& theColor = rowPtr[j];
				auto h = grid[i][j];
				if ( h >= m && h < n) {
					theColor = cv::Vec3b(255, 255, 255);
				}
                // 处理图像四周信息，连接未封闭的曲线。
                if( i == 0 || j == 0 || i == image_e.rows-1 || j == image_e.cols-1 ){
                    if ( h >= m && h < n) {
                        theColor = cv::Vec3b(255, 255, 255);
                    }
                }
			}
		}



		for (const auto& c : curves[a]) {
			for (auto cc : c) {
				cv::circle(image_e, cv::Point((int)cc.x, (int)cc.y), 1, cv::Scalar(255, 255, 255), -1);
			}
		}
//		for (const auto& c : curves[a + 1]) {
//			for (auto cc : c) {
//				cv::circle(image_e, cv::Point((int)cc.x, (int)cc.y), 1, cv::Scalar(0, 0, 0), -1);
//			}
//		}
        cv::imwrite(std::to_string(a) + "_fill_image.jpg", image_e);
		images.emplace_back(image_e);
	}
    // findContours
    std::vector<std::vector<std::vector<cv::Point>>> contourss;
    std::vector<std::vector<cv::Vec4i>> hierarchys;
    int aa=0;
    for(const auto& im : images){

        std::vector<std::vector<cv::Point>> contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::Mat im_2,im_t;
        cv::cvtColor(im,im_2,cv::COLOR_BGR2GRAY);
        cv::threshold(im_2,im_t,144,255,cv::THRESH_BINARY);
        cv::findContours(im_t, contours,hierarchy,cv::RETR_TREE,cv::CHAIN_APPROX_TC89_L1);
        contourss.emplace_back(contours);
        hierarchys.emplace_back(hierarchy);
        // return 1. contours[i] 表示第i个轮廓
        cv::Mat image_e(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0, 0, 0));
        for(const auto& c:contours){
            auto d = cv::contourArea(c, true);
            if(d<=0) {
                for (const auto &p: c) {
                    cv::circle(image_e, p, 1, cv::Scalar(255, 0, 0), -1);
                }
            }else{
                for (const auto &p: c) {
                    cv::circle(image_e, p, 1, cv::Scalar(0, 255, 0), -1);
                }
            }
        }
        cv::imwrite(std::to_string(aa++) + "_line_image.jpg", image_e);
    }

    sortTopology(contourss, hierarchys);
	// getTopology(curves, images);
}
/***
 *  input :
 *   1. contourss
 *   2. hierarchys
 *   3. fatherLineNum
 *   4. thisLineNum
 *   5. newContour
 *   6. newHierarchy
 *   7. level
 */
void regroup(std::vector<cv::Vec4i>& new_hierarchy,
             std::vector<std::pair<std::vector<cv::Point>,int>>& new_contours,
             int fatherLine, int thisLine,int n,
             const std::vector<std::vector<std::vector<cv::Point>>>& contourss,
             const std::vector<std::vector<cv::Vec4i>>& hierarchys){
    // 重组拓扑结构。
    // 获取父轮廓和子轮廓
//    auto fatherCurve = contourss[n][fatherLine];
//    auto childCurve = contourss[n + 1][thisLine];
//    auto fatherTopo = hierarchys[n][fatherLine]; // [Next, Previous, First_Child, Parent]
//    auto childTopo = hierarchys[n + 1][thisLine];
//    int fatherID = n * 10000 + fatherLine;
//    int childID = (n+1) * 10000 + thisLine;

//    int fatherIdx,childIdx;
//    for(const auto& nc : new_contours){
//        if  (nc.second == fatherID) {
//            fatherIdx = nc.second;
//        }else if(nc.second == childID){
//
//        }
//    }
//
//    // 将父轮廓和子轮廓添加到新的轮廓集合中
//    int fatherIdx = new_contours.size();
//    new_contours.push_back(fatherCurve);
//    int childIdx = new_contours.size();
//    new_contours.push_back(childCurve);
//
//    // 更新层次结构
//    cv::Vec4i newFatherTopo = {-1, -1, childIdx, -1}; // 父轮廓的第一个子轮廓指向子轮廓
//    cv::Vec4i newChildTopo = {-1, -1, -1, fatherIdx}; // 子轮廓的父轮廓指向父轮廓
//
//    // 如果父轮廓在原始层次结构中已有子轮廓，需处理
//    if (fatherTopo[2] != -1) {
//        // 将原始第一个子轮廓的索引放入新层次结构
//        int origFirstChildIdx = new_contours.size();
//        new_contours.push_back(contourss[n][fatherTopo[2]]);
//        newFatherTopo[2] = origFirstChildIdx; // 更新父轮廓的第一个子轮廓
//        new_hierarchy.push_back({-1, -1, -1, fatherIdx}); // 原始子轮廓的父轮廓指向 fatherIdx
//    }
//
//    // 添加新的层次信息
//    new_hierarchy.push_back(newFatherTopo);
//    new_hierarchy.push_back(newChildTopo);
}



void ImageProcess::sortTopology(std::vector<std::vector<std::vector<cv::Point>>> contourss,
                                std::vector<std::vector<cv::Vec4i>> hierarchys) {
    std::vector<cv::Vec4i> new_hierarchy;
    std::vector<std::vector<cv::Point>> new_contours;
    cv::Mat image_all(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0, 0, 0));

    // contourss.size();
    // test for 3，4 level
    for(auto n=0;n<contourss.size()-1;++n) {

        auto hi_2 = hierarchys[n];
        auto hi_3 = hierarchys[n + 1];
        auto co_2 = contourss[n];
        auto co_3 = contourss[n + 1];

        cv::Mat image_e(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0, 0, 0));
        int help_count = 1;
        for (const auto &c: co_2) {
            for (const auto &p: c) {
                auto r = static_cast<unsigned char>(help_count % 256);
                auto g = static_cast<unsigned char>((help_count / 256) % 256);
                auto b = static_cast<unsigned char>((help_count / (256 * 256)) % 256);
                cv::circle(image_e, p, 1, cv::Scalar(r, g, b), -1);
                cv::circle(image_all, p, 1, cv::Scalar(255, 255, 255), -1);

            }
            help_count++;
        }
        help_count = 1;
        for (const auto &c: co_3) {
            bool visited = false;
            for (const auto &p: c) {
                if (!visited && (image_e.at<cv::Vec3b>(p.y, p.x) != cv::Vec3b(0, 0, 0))) {
                    auto color = image_e.at<cv::Vec3b>(p.y, p.x);
                    auto fatherLine = color[0] + (color[1] * 256) + (color[2] * 256 * 256);
                    auto thisLine = help_count;
                    //regroup(new_hierarchy, new_contours, fatherLine - 1, thisLine - 1, n, contourss, hierarchys);
                    visited = true;
                }
                cv::circle(image_all, p, 1, cv::Scalar(255, 255, 255), -1);
            }
            help_count++;
        }
        cv::imwrite(image.getImageName() + "test_line_image.jpg", image_all);
    }

    // 检查拓扑内容和结构是否正确
    // 检查内容
    cv::Mat image_check(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0, 0, 0));
    for(const auto& l : new_contours ){
        for(const auto& p : l){
            cv::circle(image_check, p, 1, cv::Scalar(255,255,255), -1);
        }
    }
    cv::imwrite(image.getImageName()+"testcheck_line_image.jpg", image_check);
    int a = 0;
    std::cin >> a;
}


// 根据findContours函数给出的内外圈信息，仅边界信息处理有误。借助之前写的，根据边界值是否在曲线内部，封闭边界不连续的曲线。

/***
 * TODO: 根据上下层的连接关系，获取拓扑关系。
 *  1. 选择3，4两层，借助contours的结构，匹配曲线，并形成新的拓扑关系。
 *  2. 可以形成拓扑关系。
 *  3.1 遍历时，如果有image_e.at<cv::Vec3b>(p.y, p.x) != cv::Vec3b(0, 0, 0)的情况：
 *      3.1.1 说明该绘制曲线为两层level共享的一条边界，也就是下层的内边界，上层的外边界。
 *      3.1.2 找到该曲线在下层的父结点,如果没有父节点，该节点为父节点。
 *      3.1.3 找到该曲线在上层的子节点，如果没有子节点，该节点为叶子节点。
 *      3.1.4 将以上节点形成爷父子节点在拓扑中存储。
 */


//void ImageProcess::sortTopology(std::vector<std::vector<std::vector<cv::Point>>> contourss,
//                  std::vector<std::vector<cv::Vec4i>> hierarchys) {
//    std::vector<std::vector<cv::Point>> new_contours;
//    std::vector<cv::Vec4i> new_hierarchy;
//
//    // 辅助函数：将一层轮廓及其内部拓扑关系加入全局结构
//    auto addLayer = [&](int layerIdx) {
//        auto& contours = contourss[layerIdx];
//        auto& hierarchy = hierarchys[layerIdx];
//        int startIdx = new_contours.size();
//
//        // 添加轮廓并初始化层次结构
//        for (const auto& contour : contours) {
//            new_contours.push_back(contour);
//            new_hierarchy.push_back(cv::Vec4i(-1, -1, -1, -1)); // 初始为空
//        }
//
//        // 映射层内拓扑关系到全局索引
//        for (int i = 0; i < contours.size(); ++i) {
//            int globalIdx = startIdx + i;
//            auto& topo = hierarchy[i];
//            // 设置 Next, Previous, First_Child, Parent
//            if (topo[0] != -1) new_hierarchy[globalIdx][0] = startIdx + topo[0]; // Next
//            if (topo[1] != -1) new_hierarchy[globalIdx][1] = startIdx + topo[1]; // Previous
//            if (topo[2] != -1) new_hierarchy[globalIdx][2] = startIdx + topo[2]; // First_Child
//            if (topo[3] != -1) new_hierarchy[globalIdx][3] = startIdx + topo[3]; // Parent
//        }
//        return startIdx;
//    };
//
//    // 辅助函数：设置层间父子关系
//    auto setParentChild = [&](int fatherIdx, int childIdx) {
//        // 将 childIdx 设为 fatherIdx 的子轮廓
//        if (new_hierarchy[fatherIdx][2] == -1) {
//            new_hierarchy[fatherIdx][2] = childIdx; // First_Child
//        } else {
//            // 如果已有子轮廓，追加到兄弟链
//            int lastChild = new_hierarchy[fatherIdx][2];
//            while (new_hierarchy[lastChild][0] != -1) {
//                lastChild = new_hierarchy[lastChild][0];
//            }
//            new_hierarchy[lastChild][0] = childIdx;  // Next
//            new_hierarchy[childIdx][1] = lastChild;  // Previous
//        }
//        new_hierarchy[childIdx][3] = fatherIdx; // Parent
//    };
//
//    // Step 1: 处理每一层并保留层内拓扑关系
//    std::vector<int> layerStartIndices;
//    for (int n = 0; n < contourss.size(); ++n) {
//        int startIdx = addLayer(n);
//        layerStartIndices.push_back(startIdx);
//    }
//
//    // Step 2: 建立层间关系
//    for (int n = 0; n < contourss.size() - 1; ++n) {
//        int start_n = layerStartIndices[n];
//        int end_n = (n + 1 < contourss.size()) ? layerStartIndices[n + 1] : new_contours.size();
//        int start_n1 = layerStartIndices[n + 1];
//        int end_n1 = (n + 2 < contourss.size()) ? layerStartIndices[n + 2] : new_contours.size();
//
//        for (int i = start_n; i < end_n; ++i) {
//            for (int j = start_n1; j < end_n1; ++j) {
//                // 检查第 n+1 层轮廓 j 是否在第 n 层轮廓 i 内部
//                if (cv::pointPolygonTest(new_contours[i], new_contours[j][0], false) > 0) {
//                    setParentChild(i, j);
//                }
//            }
//        }
//    }
//
//    // 结果：new_contours 和 new_hierarchy 包含完整的拓扑信息
//    generateTopologyGraph(new_hierarchy);
//}

void ImageProcess::to_train_image()
{
        auto levels = image.valueScheme[0];
        // 获取图像指针
        for(int level = 0 ; level < levels.size()-1; ++ level) {
            cv::Mat bezierFittedImage(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(0,0,0));

            auto h1 = levels[level];
            auto h2 = levels[level + 1];

            for (int i = 0; i < image.imageSize; ++i) {
                auto *rowPtr = bezierFittedImage.ptr<cv::Vec3b>(i);
                for (int j = 0; j < image.imageSize; ++j) {
                    // 计算像素点的位置和颜色
                    double x = image.offset + image.cellSize * i;
                    double y = image.offset + image.cellSize * j;
                    auto pi = image.getPointInfo(x, y);
                    if (pi.value >= h1 && pi.value < h2) {
                        rowPtr[j] = cv::Vec3b(255,255,255);
                    }
                }
            }
            //std::string outputDir = "D:/coding/train_data/";// path/to/your/folder/"; // 必须以斜杠结尾
            std::string n = image.getImageName() + "_" + std::to_string(level) + "_originalImage.jpg";
            std::cout << "a original image is generated as " + n << std::endl;
            cv::imwrite(n, bezierFittedImage);
        }
    }

