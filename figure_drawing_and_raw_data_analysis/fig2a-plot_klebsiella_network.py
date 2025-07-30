import pandas as pd
import networkx as nx
from pyvis.network import Network
import random

def random_color():
    return "#" + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)])

def plot_klebsiella_network(csv_file, output_html):
    df = pd.read_csv(csv_file)
    print(f"总记录数: {len(df)}")
    
    # 取proN1和proN2为节点
    node1_list = df['proN1'].astype(str).tolist()
    node2_list = df['proN2'].astype(str).tolist()
    
    G = nx.Graph()
    for n1, n2 in zip(node1_list, node2_list):
        if n1 and n2:
            G.add_edge(n1, n2)
    G.remove_edges_from(nx.selfloop_edges(G))

    # 生成连通子图颜色映射
    components = list(nx.connected_components(G))
    color_map = {}
    for comp in components:
        c = random_color()
        for node in comp:
            color_map[node] = c

    net = Network(notebook=False, width="100%", height="800px", bgcolor="#ffffff", font_color="black")
    # 增加可配置选项
    options = {
        "physics": {
            "enabled": True,
            "barnesHut": {
                "gravitationalConstant": -1500,
                "centralGravity": 0.2,
                "springLength": 120,
                "springConstant": 0.05,
                "damping": 0.15,
                "avoidOverlap": 0.8
            },
            "minVelocity": 0.5,
            "maxVelocity": 30,
            "solver": "barnesHut",
            "stabilization": {
                "enabled": True,
                "iterations": 1500,
                "updateInterval": 25,
                "fit": True
            }
        },
        "interaction": {
            "hover": True,
            "navigationButtons": True,
            "keyboard": True,
            "dragNodes": True,
            "dragView": True,
            "zoomView": True,
            "multiselect": True,
            "selectable": True
        },
        "layout": {
            "improvedLayout": True,
            "hierarchical": {
                "enabled": False
            }
        },
        "configure": {
            "enabled": True,
            "filter": ["physics", "nodes", "edges", "layout", "interaction"]
        }
    }
    net.options = options
    net.from_nx(G)
    for node in net.nodes:
        node_id = node.get("id", node.get("label"))
        node["size"] = 20
        node["label"] = node_id
        node["color"] = color_map.get(node_id, "#4682B4")
        node["title"] = node_id
    for edge in net.edges:
        edge["color"] = "#808080"
        edge["width"] = 2
        edge["title"] = "Protein-Protein Interaction"
    net.show(output_html, notebook=False)
    print(f"网络图已保存至 {output_html}")

if __name__ == "__main__":
    csv_file = 'verifytable/heterodimer_batch123_matched_LIS_LIA_filtered_4_Klebsiella_with_proN.csv'
    output_html = 'klebsiella_protein_network.html'
    plot_klebsiella_network(csv_file, output_html) 