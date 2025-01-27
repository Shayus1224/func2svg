#include <stack>
#include "ImageProcess.h"
#include <alglib/interpolation.h>

int ImageProcess::get_case(double v1, double v2, double v3, double v4, double level)
{
	return (v1 >= level) | ((v2 >= level) << 1) | ((v3 >= level) << 2) | ((v4 >= level) << 3);
}

Point ImageProcess::interpolate(double x1, double y1, double x2, double y2, double v1, double v2, double level)
{
	double t = (level - v1) / (v2 - v1);
	return { x1 + t * (x2 - x1), y1 + t * (y2 - y1) };
}

double ImageProcess::distance(const Point& p1, const Point& p2)
{
	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	return std::sqrt(dx * dx + dy * dy);
}

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

std::vector<Point> ImageProcess::marching_squares_points(const std::vector<std::vector<double>>& grid, double level)
{
	std::vector<Point> res;
	int rows = grid.size();
	int cols = grid[0].size();

	for (double i = 0; i < rows - 1; ++i) {
		for (double j = 0; j < cols - 1; ++j) {
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

std::vector<std::vector<Point>> ImageProcess::connectPoints(const std::vector<Point>& points)
{
	int n = points.size();
	std::vector<bool> visited(n, false);
	std::vector<std::vector<Point>> curves;

	auto compareDist = [&points](const std::pair<double, int>& a, const std::pair<double, int>& b) {
		return a.first > b.first;
		};

	for (int i = 0; i < n; ++i) {
		if (!visited[i]) {
			std::vector<Point>curve;
			std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, decltype(compareDist)> pq(compareDist);
			pq.push({ 0, i });

			while (!pq.empty()) {
				auto dist = pq.top();
				auto current = dist.second;
				pq.pop();

				if (visited[current]) continue;
				visited[current] = true;
				curve.push_back(points[current]);

				for (int j = 0; j < n; ++j) {
					if (!visited[j]) {
						double dist = distance(points[current], points[j]);
						if (dist < 5) {
							pq.push({ dist, j });
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

void ImageProcess::getTopology(const std::vector<std::vector<std::vector<Point>>>& curves, const std::vector<cv::Mat>& images)
{
	std::cout << "开始整理拓扑" << std::endl;
	Topology topology;
	// auto colorScheme = getColorValues();
	std::vector<std::vector<std::vector<Point>>> colorDotss;
	// 初始化所有拓扑节点
	for (int t = 0; t < images.size() - 1; ++t) {
		auto imaget = images[t];
		cv::Mat grayImg;
		// 将彩色图像转换为灰度图像
		cv::cvtColor(imaget, grayImg, cv::COLOR_BGR2GRAY);
		// 遍历图像
		std::vector<Point>colorPoint;
		for (int i = 0; i < imaget.rows; ++i) {
			const auto* row_ptr = imaget.ptr<uchar>(i); // 获取第i行的指针
			for (int j = 0; j < imaget.cols; ++j) {
				if (row_ptr[j] == 0) { // 假设条件是像素值等于255
					colorPoint.emplace_back(j, i); // 存储满足条件的点的位置
				}
			}
		}
		auto colorDots = connectPoints(colorPoint);
		colorDotss.emplace_back(colorDots);
		for (int i = 0; i < colorDots.size(); ++i) {
			int id = ((t + 1) * 10000) + i;
			cv::Scalar color = image.Scheme[image.colorChoice][i];
			topology.addNode(id, color);
		}
		std::cout << t << "层节点建立完成, 现在有" << topology.nodes.size() << std::endl;
	}
	// 整理所有的拓扑，一个vertex + 一个edge
	for (int t = 0; t < images.size() - 1; ++t) {
		auto colorDots = colorDotss[t];
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
	for (auto& l : topology.edges) {
		auto points = l.outerLine;
		std::set<Point> uniquePoints(points.begin(), points.end());
		points.assign(uniquePoints.begin(), uniquePoints.end());

		int size = points.size();
		std::cout << "第" << level++ << "个有" << size << "个点\n";
		// 将数据点转换为alglib格式
		cv::Mat imagee(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(255, 255, 255));
		// std::ofstream outfile("line_" + std::to_string(level) + "point.txt");

		alglib::real_1d_array x_array, y_array;

		x_array.setlength(size);
		y_array.setlength(size);
		// if (outfile.is_open()) {
			for (int i = 0; i < points.size(); ++i) {
				auto x = points[i].x;
				auto y = points[i].y;

				imagee.at<cv::Vec3b>(y, x)[0] = 0;
				imagee.at<cv::Vec3b>(y, x)[1] = 0;
				imagee.at<cv::Vec3b>(y, x)[2] = 0;
				// outfile << x << " , " << y << std::endl;
				x_array[i] = x;
				y_array[i] = y;
			}
		// } else {
		// 	std::cerr << "无法打开文件进行追加   " << std::to_string(level) << std::endl;
		// }
		// outfile.close();
		// std::string file_name = "line_" + std::to_string(level) + ".jpg";
		// cv::imwrite(file_name, imagee);
		try {

			alglib::spline1dinterpolant res1;
			alglib::spline1dbuildcubic(x_array, y_array, res1);
			l.line1 = res1;

		} catch (alglib::ap_error alglib_exception) {

			printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());

		}
	}

	topology.printTopology();
	std::string name = image.getImageName();
	topology.generateDotFile(name);
}

void ImageProcess::to_image()
{
	cv::Mat bezierFittedImage(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(255, 255, 255));
	// 获取图像指针
	for (int i = 0; i < image.imageSize; ++i) {
		auto* rowPtr = bezierFittedImage.ptr<cv::Vec3b>(i);
		for (int j = 0; j < image.imageSize; ++j) {
			// 计算像素点的位置和颜色
			double x = image.offset + image.cellSize * i;
			double y = image.offset + image.cellSize * j;
			auto pi = image.getPointInfo(x, y);

			rowPtr[j] = cv::Vec3b((uchar)pi.Color[2], (uchar)pi.Color[1], (uchar)pi.Color[0]);
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
	// 生成辅助图片
	std::vector<cv::Mat> images;
	for (int a = 0; a < curves.size() - 1; ++a) {
		double m, n;
		auto colorScheme = image.valueScheme[image.colorChoice];
		if ((a - 1) < 0) {
			m = -std::numeric_limits<double>::max();
			n = colorScheme[a];
		} else if ((a + 1) == colorScheme.size()) {
			m = colorScheme[a];
			n = std::numeric_limits<double>::max();
		} else {
			m = colorScheme[a];
			n = colorScheme[a + 1];
		}
		// auto color = getColor(colorScheme[a]);
		cv::Mat imagee(image.imageSize, image.imageSize, CV_8UC3, cv::Scalar(255, 255, 255));

		for (int i = 0; i < imagee.rows; ++i) {
			// 获取指向第 i 行的指针
			cv::Vec3b* rowPtr = imagee.ptr<cv::Vec3b>(i);
			for (int j = 0; j < imagee.cols; ++j) {
				// 获取像素值
				cv::Vec3b& theColor = rowPtr[j];
				auto h = grid[i][j];
				if (m < h && h < n) {
					theColor = cv::Vec3b(0, 0, 0);
				}
			}
		}
		for (const auto& c : curves[a]) {
			for (auto cc : c) {
				cv::circle(imagee, cv::Point(cc.x, cc.y), 1, cv::Scalar(0, 0, 0), -1);
			}
		}
		for (const auto& c : curves[a + 1]) {
			for (auto cc : c) {
				cv::circle(imagee, cv::Point(cc.x, cc.y), 1, cv::Scalar(0, 0, 0), -1);
			}
		}
		images.emplace_back(imagee);
	}

	getTopology(curves, images);
}

