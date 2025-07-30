import pandas as pd
import networkx as nx
from pyvis.network import Network

def plot_combined_network(csv_file, output_html):
    df = pd.read_csv(csv_file, low_memory=False)
    node1_list = df['combined_1'].astype(str).tolist()
    node2_list = df['combined_2'].astype(str).tolist()

    G = nx.Graph()
    for n1, n2 in zip(node1_list, node2_list):
        if n1 and n2:
            G.add_edge(n1, n2)

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

    node1_set = set(node1_list)
    node2_set = set(node2_list)
    for node in net.nodes:
        node_id = node.get("id", node.get("label"))
        if node_id in node1_set:
            node["color"] = "#4682B4"  # 蓝色
        elif node_id in node2_set:
            node["color"] = "#D62728"  # 红色
        else:
            node["color"] = "#888888"
        node["size"] = 20
        node["label"] = node_id
        node["title"] = node_id

    for edge in net.edges:
        edge["color"] = "#808080"
        edge["width"] = 2
        edge["title"] = "Protein-Protein Interaction"

    net.show(output_html, notebook=False)
    print(f"网络图已保存至 {output_html}")

if __name__ == "__main__":
    plot_combined_network(
        "fig444/115w_scores_plus_LIS_0409_human_virus_combined_columns_filtered.csv",
        "fig444/combined_protein_network.html"
    ) 