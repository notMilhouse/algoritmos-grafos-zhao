import networkx as nx
import numpy as np

def create_network(data):
    """Create network graph from STRING data.
    
    Args:
        data (list): List of interaction data dictionaries
        
    Returns:
        nx.Graph: NetworkX graph object
    """
    G = nx.Graph()
    for interaction in data:
        source = interaction['preferredName_A']
        target = interaction['preferredName_B']
        score = float(interaction['score']) / 1000
        G.add_edge(source, target, weight=score)
    return G

def analyze_network_metrics(G):
    """Calculate network metrics for the network.
    
    Args:
        G (nx.Graph): NetworkX graph object
        
    Returns:
        dict: Dictionary of network metrics
    """
    # Basic metrics
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    avg_degree = sum(dict(G.degree()).values()) / num_nodes
    density = nx.density(G)
    
    # Calculate clustering coefficient
    avg_clustering = nx.average_clustering(G)
    
    # Calculate path length metrics if the graph is connected
    try:
        diameter = nx.diameter(G)
        avg_path_length = nx.average_shortest_path_length(G)
    except nx.NetworkXError:
        # If graph is not connected, use largest connected component
        largest_cc = max(nx.connected_components(G), key=len)
        largest_subgraph = G.subgraph(largest_cc)
        diameter = nx.diameter(largest_subgraph)
        avg_path_length = nx.average_shortest_path_length(largest_subgraph)
    
    metrics = {
        "num_nodes": num_nodes,
        "num_edges": num_edges,
        "avg_degree": avg_degree,
        "density": density,
        "avg_clustering": avg_clustering,
        "diameter": diameter,
        "avg_path_length": avg_path_length
    }
    
    return metrics

def get_interaction_type(interaction):
    """Determine interaction type based on evidence scores.
    
    Args:
        interaction (dict): Interaction data dictionary
        
    Returns:
        tuple: (interaction_type, color_code)
    """
    escore = float(interaction['escore'])
    dscore = float(interaction['dscore'])
    tscore = float(interaction['tscore'])
    
    if escore > 0.7:
        return 'physical', '#7E57C2'
    elif dscore > 0.7:
        return 'functional', '#4CAF50'
    elif tscore > 0.7:
        return 'text_mining', '#F44336'
    return 'other', '#90A4AE' 