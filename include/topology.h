#pragma once
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <alglib/ap.h>
#include <alglib/interpolation.h>

#include "QuasiRegularModel.h"
#include "ImageProcess.h"

struct spline_info {
    int n;
    double* knots;
    double* a, * b, * c, * d;
};

struct Node {
//    // 截面对应的高度值
//    double level;
//    // 在截面上曲线的标号
//    int id;
//    // 邻居节点（上下层）同层不需要连接
//    std::unordered_set<int> neighbors;
    int id;
    cv::Scalar color;
    std::unordered_set<int> neighbors;
};
struct Edge {
    int startNode;
    int endNode;
    int id;
    std::vector<Point> outerLine; // 目前用离散点表示，之后用曲线拟合
    alglib::spline1dinterpolant line1, line2; // 拟合后的曲线信息
    spline_info l1, l2;
    Edge(int s, int e, int id, std::vector<Point> l) : startNode(s), endNode(e), id(id), outerLine(std::move(l))
    {
    }
};


class Topology {
public:

    std::unordered_map<int, Node> nodes; // 节点映射，使用节点id作为键
    std::vector<Edge> edges; // 边的列表

    // 添加节点
    // return 1:添加成功，return 2，添加失败，已存在。
    bool addNode(int id, cv::Scalar color)
    {
        if (nodes.find(id) != nodes.end()) {
            std::cerr << "Node with id " << id << " already exists." << std::endl;
            return 0;
        }
        nodes[id] = { id, std::move(color), {} };
        return 1;
    }

    // 添加边
    void addEdge(int startNode, int endNode, int edgeId, const std::vector<Point>& outerLine)
    {
        if (nodes.find(startNode) == nodes.end()) {
            std::cerr << "One node not found." << std::endl;
            return;
        }
        edges.emplace_back(startNode, endNode, edgeId, outerLine);
        // nodes[startNode].neighbors.push_back(endNode);
    }
    // 连接拓扑
    void connectTopo()
    {
        for (auto& e1 : this->edges) {
            for (auto& e2 : this->edges) {
                if (e1.id == e2.id && e1.startNode != e2.startNode) {
                    nodes[e1.startNode].neighbors.insert(e2.startNode);
                    e1.endNode = e2.startNode;
                }
            }
        }
    }


    // 删除多余的拓扑连接线
    void removeRedundantEdges()
    {
        try {
            std::unordered_set<int> uniqueEdgeIds;
            auto it = edges.begin();
            while (it != edges.end()) {
                if (uniqueEdgeIds.find(it->id) != uniqueEdgeIds.end()) {
                    // 找到重复的边，删除
                    it = edges.erase(it);
                } else {
                    uniqueEdgeIds.insert(it->id);
                    ++it;
                }
            }
        } catch (const std::exception& e) {
            // 捕获并处理所有std::exception派生的异常
            std::cerr << "An error occurred: " << e.what() << std::endl;
        } catch (...) {
            // 捕获并处理所有其他类型的异常
            std::cerr << "An unknown error occurred." << std::endl;
        }
    }
    // 获取节点映射，使用节点id作为键
    const std::unordered_map<int, Node>& getNodes() const
    {
        return nodes;
    }

    // 获取边的列表
    const std::vector<Edge>& getEdges() const
    {
        return edges;
    }
    // 打印拓扑结构
    void printTopology() const
    {
        std::ofstream outfile("test.txt");
        std::ofstream outfile_line("lines.txt");
        if (!outfile.is_open()) {
            std::cerr << "Error opening file for writing!" << std::endl;
            return;
        }
        outfile << "Nodes:\n";
        for (const auto& in : nodes) {
            outfile << "Node ID: " << in.first << ", Color: (" << in.second.color[0] << ", " << in.second.color[1] << ", " << in.second.color[2] << "), Neighbors: ";
            for (const auto& neighbor : in.second.neighbors) {
                outfile << neighbor << " ";
            }
            outfile << "\n";
        }
        outfile << "\nEdges:\n";
        for (const auto& edge : edges) {
            outfile << "Edge ID: " << edge.id << ", Start Node: " << edge.startNode << ", End Node: " << edge.endNode << std::endl; //  << "    ";

            std::string serializedSpline;
            auto& spline1 = edge.line1;
            alglib::spline1dserialize(spline1, serializedSpline);
			std::cout << serializedSpline << std::endl;
            outfile_line<< "Edge ID: " << edge.id << std::endl << serializedSpline << std::endl << std::endl;

            auto& spline2 = edge.line2;
            alglib::spline1dserialize(spline2, serializedSpline);
            std::cout << serializedSpline << std::endl;
            outfile_line<< "Edge ID: " << edge.id << std::endl << serializedSpline << std::endl << std::endl;

        }
        outfile.close();
        outfile_line.close();
    }

    void generateDotFile(std::string file_loc)
    {
        std::string filename = file_loc + ".dot";
        std::ofstream dotFile(filename);
        if (!dotFile.is_open()) {
            std::cerr << "Unable to open file: " << filename << std::endl;
            return;
        }

        // 写入 dot 格式内容
        dotFile << "graph G {" << std::endl;
        
        // 写入节点
        for (const auto& node : nodes) {
            dotFile << "  " << node.first << " [style=filled, fillcolor=\"#"
                << std::hex << node.second.color[0] << node.second.color[1] << node.second.color[2] << "\"];" << std::endl;
        }

        // 写入边
        for (const auto& edge : edges) {
            dotFile << "  " << edge.startNode << " -- " << edge.endNode << ";" << std::endl;
        }

        dotFile << "}" << std::endl;
        dotFile.close();
        std::cout << "Dot file generated: " << filename << std::endl;
    }
};