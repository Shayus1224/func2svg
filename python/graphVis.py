import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

# 解析节点数据
def parse_nodes(file_content):
    nodes = []
    for line in file_content:
        if line.startswith("Node ID:"):
            parts = line.split(",")
            node_id = int(parts[0].split(":")[1].strip())
            color = tuple(map(int, parts[1].split("(")[1].split(")")[0].split(",")))
            neighbors = list(map(int, parts[2].split(":")[1].strip().split())) if "Neighbors" in line else []
            nodes.append({"Node ID": node_id, "Color": color, "Neighbors": neighbors})
    return nodes

# 解析边数据
def parse_edges(file_content):
    edges = []
    for line in file_content:
        if line.startswith("Edge ID:"):
            parts = line.split(",")
            start_node = int(parts[1].split(":")[1].strip())
            end_node = int(parts[2].split(":")[1].strip())
            edges.append((start_node, end_node))
    return edges

# 读取文件
with open("test.txt", "r") as file:
    content = file.readlines()

# 提取节点和边信息
nodes = parse_nodes(content)
edges = parse_edges(content)

# 创建图
G = nx.Graph()

# 添加节点
for node in nodes:
    color_hex = to_hex([c / 255 for c in node["Color"][:3]])  # 转换颜色到十六进制
    G.add_node(node["Node ID"])#, color=color_hex)

# 添加边
G.add_edges_from(edges)

# 设置节点颜色
# node_colors = [G.nodes[node]["color"] for node in G.nodes]

# 绘制图
plt.figure(figsize=(15, 15))
pos = nx.spring_layout(G, seed=42)  # 使用弹簧布局
nx.draw(G, pos, with_labels=True, node_size=500, edge_color="gray")#node_color=node_colors,
plt.title("Graph Visualization", fontsize=10)
plt.show()
