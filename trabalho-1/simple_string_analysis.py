#!/usr/bin/env python3
# Simple STRING network analysis script
# Fetches protein interaction data and creates a non-radial visualization

import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

def get_string_data(protein, species=9606, required_score=400, limit=100):
    """Fetch protein interaction data from STRING database."""
    string_api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": protein,
        "species": species,
        "network_type": "physical",
        "required_score": required_score,
        "limit": limit
    }
    response = requests.get(string_api_url, params=params)
    return response.json()

def create_network(data):
    """Create network from interaction data."""
    G = nx.Graph()
    for interaction in data:
        source = interaction['preferredName_A']
        target = interaction['preferredName_B']
        # Use combined score as edge weight, normalized to [0,1]
        score = float(interaction['score']) / 1000
        G.add_edge(source, target, weight=score)
    return G

def analyze_network(G):
    """Calculate and print basic network metrics."""
    # Basic metrics
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    density = nx.density(G)
    avg_degree = sum(dict(G.degree()).values()) / num_nodes
    
    # Calculate clustering coefficient
    avg_clustering = nx.average_clustering(G)
    
    # Calculate path length metrics if the graph is connected
    if nx.is_connected(G):
        diameter = nx.diameter(G)
        avg_path_length = nx.average_shortest_path_length(G)
    else:
        # If graph is not connected, use largest connected component
        largest_cc = max(nx.connected_components(G), key=len)
        largest_subgraph = G.subgraph(largest_cc)
        diameter = nx.diameter(largest_subgraph)
        avg_path_length = nx.average_shortest_path_length(largest_subgraph)
    
    # Calculate centrality
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    
    # Print metrics
    print(f"Network Analysis Results for {protein}:")
    print(f"Nodes (proteins): {num_nodes}")
    print(f"Edges (interactions): {num_edges}")
    print(f"Network density: {density:.4f}")
    print(f"Average degree: {avg_degree:.2f}")
    print(f"Average clustering coefficient: {avg_clustering:.4f}")
    print(f"Network diameter: {diameter}")
    print(f"Average path length: {avg_path_length:.4f}")
    
    # Print top proteins by centrality
    print("\nTop 5 proteins by degree centrality:")
    for node, cent in sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"{node}: {cent:.4f}")
    
    print("\nTop 5 proteins by betweenness centrality:")
    for node, cent in sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"{node}: {cent:.4f}")
    
    return degree_centrality, betweenness_centrality

def create_non_radial_visualization(G, data, output_path):
    """Create a non-radial network visualization."""
    plt.figure(figsize=(14, 14))
    
    # Create a better layout that avoids the "hairball" effect
    # 1. Start with a random layout to break any circular patterns
    pos = nx.random_layout(G, seed=42)
    
    # 2. Apply Fruchterman-Reingold layout with tuned parameters
    pos = nx.fruchterman_reingold_layout(G, k=0.15, pos=pos, iterations=50, seed=42)
    
    # 3. Apply spring layout as a refinement with low k to keep nodes close
    pos = nx.spring_layout(G, pos=pos, k=0.1, iterations=50, seed=42)
    
    # 4. Ensure the central protein isn't exactly at center to reduce radial effect
    if protein in G.nodes():
        pos[protein] = np.array([0.1, 0.1])  # Slightly off-center
    
    # Get centrality for node sizes
    betweenness = nx.betweenness_centrality(G)
    degree = dict(G.degree())
    max_degree = max(degree.values())
    
    # Define node colors and sizes
    node_colors = []
    node_sizes = []
    
    for node in G.nodes():
        if node == protein:  # Central protein
            node_colors.append('#FFD700')  # Gold
            node_sizes.append(3000)
        elif degree[node] > max_degree/2:  # High-degree nodes
            node_colors.append('#FF6347')  # Tomato
            node_sizes.append(1000 + 2000 * betweenness[node])
        else:  # Other nodes
            node_colors.append('#1E90FF')  # Dodger blue
            node_sizes.append(500 + 1500 * betweenness[node])
    
    # Edge properties
    edge_colors = []
    edge_widths = []
    
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                # Color edges based on score
                score = float(interaction['score']) / 1000
                if score > 0.9:
                    edge_colors.append('#4CAF50')  # Green for high confidence
                elif score > 0.7:
                    edge_colors.append('#2196F3')  # Blue for medium confidence
                else:
                    edge_colors.append('#9E9E9E')  # Gray for lower confidence
                edge_widths.append(0.5 + 2.5 * score)
                break
    
    # Draw network
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5)
    
    # Draw labels only for important nodes
    important_nodes = []
    # Add central protein
    if protein in G.nodes():
        important_nodes.append(protein)
    
    # Add top 10 nodes by degree
    important_nodes.extend([node for node, d in sorted(degree.items(), 
                                                      key=lambda x: x[1], 
                                                      reverse=True)[:10] 
                           if node != protein])
    
    labels = {node: node for node in important_nodes}
    
    # Draw labels with white background for better visibility
    for node, label in labels.items():
        font_size = 14 if node == protein else 10
        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9)
        plt.text(pos[node][0], pos[node][1], label, 
                horizontalalignment='center', fontsize=font_size, fontweight='bold',
                bbox=bbox_props)
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label=f'{protein} (Central Protein)',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', label='High Connectivity',
                  markerfacecolor='#FF6347', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Low Connectivity',
                  markerfacecolor='#1E90FF', markersize=10),
        plt.Line2D([0], [0], color='#4CAF50', label='High Confidence (>0.9)', linewidth=3),
        plt.Line2D([0], [0], color='#2196F3', label='Medium Confidence (>0.7)', linewidth=3),
        plt.Line2D([0], [0], color='#9E9E9E', label='Lower Confidence', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', 
              title='Network Elements', fontsize=10,
              bbox_to_anchor=(1.1, 1.0))
    
    # Add title with protein name
    plt.title(f"{protein} Protein Interaction Network\nNode size: Importance | Edge width: Confidence score",
              fontsize=16, fontweight='bold', pad=20)
    
    plt.axis('off')
    plt.tight_layout()
    
    # Save and show visualization
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.show()
    print(f"Visualization saved as {output_path}")

if __name__ == "__main__":
    # Set protein of interest
    protein = "INSR"  # Insulin Receptor
    print(f"Analyzing protein: {protein}")
    
    # Fetch data
    data = get_string_data(protein)
    print(f"Fetched {len(data)} interactions")
    
    # Create and analyze network
    G = create_network(data)
    degree_centrality, betweenness_centrality = analyze_network(G)
    
    # Create and save visualization
    create_non_radial_visualization(G, data, f"{protein}_network.png") 