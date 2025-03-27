import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from network_analysis import get_interaction_type

def create_interaction_evidence_network_visualization(G, data, output_path):
    """Create visualization focusing on evidence types and confidence scores.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Layout with INSR centered
    pos = nx.spring_layout(G, k=3, iterations=200, scale=2.0, seed=42)
    pos["INSR"] = (0, 0)
    
    # Node sizes based on betweenness centrality
    betweenness = nx.betweenness_centrality(G)
    node_sizes = [8000 if node == "INSR" else 1000 + 4000 * betweenness[node] for node in G.nodes()]
    node_colors = ['#FFD700' if node == "INSR" else '#1F77B4' for node in G.nodes()]
    
    # Edge properties
    edge_colors = []
    edge_widths = []
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                edge_widths.append(1 + 3 * float(interaction['score']) / 1000)
                break
    
    # Draw network
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.5)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, edgecolors='white', linewidths=2)
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
    
    # Special INSR label
    plt.annotate("INSR (Central Protein)", xy=pos["INSR"], xytext=(0, 0.1),
                ha='center', va='bottom', fontsize=12, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='#FFD700', linewidth=2))
    
    # Legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR (Central Protein)',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3),
        plt.Line2D([0], [0], color='#90A4AE', label='Other Evidence', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
              title='Network Elements', title_fontsize=12, fontsize=10)
    
    plt.title("INSR Protein Interaction Network - Evidence Types\nNode size: Betweenness centrality | Edge width: Confidence score",
              fontsize=16, fontweight='bold', pad=20)
    plt.margins(0.3)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_community_structure_network_visualization(G, data, output_path):
    """Create visualization focusing on community structure, organizing proteins into
    distinct communities displayed in a grid layout.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Detect communities using Louvain method
    communities = nx.community.louvain_communities(G)
    communities = sorted(communities, key=len, reverse=True)
    
    # Calculate grid dimensions
    n_communities = len(communities)
    grid_size = int(np.ceil(np.sqrt(n_communities)))
    
    # Initialize positions
    pos = {}
    box_size = 2.0
    spacing = box_size + 1.0
    
    # Assign positions for each community
    for idx, community in enumerate(communities):
        grid_x = idx % grid_size
        grid_y = idx // grid_size
        subgraph = G.subgraph(community)
        sub_pos = nx.spring_layout(subgraph, k=1.0, iterations=50)
        
        for node in subgraph:
            x, y = sub_pos[node]
            pos[node] = (
                x * box_size + grid_x * spacing - (grid_size * spacing / 2),
                y * box_size + grid_y * spacing - ((grid_size * spacing / 2))
            )
    
    # Assign colors to communities
    community_colors = plt.cm.tab20(np.linspace(0, 1, len(communities)))
    node_colors = []
    for node in G.nodes():
        for idx, community in enumerate(communities):
            if node in community:
                node_colors.append(community_colors[idx])
                break
    
    # Draw network
    nx.draw_networkx_edges(G, pos, width=1, edge_color='gray', alpha=0.5)
    nx.draw_networkx_nodes(G, pos, node_size=1000, node_color=node_colors, edgecolors='white', linewidths=1)
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
    
    # Draw community boxes
    for idx, community in enumerate(communities):
        grid_x = idx % grid_size
        grid_y = idx // grid_size
        
        # Calculate box coordinates
        box_left = grid_x * spacing - (grid_size * spacing / 2) - box_size / 1.5
        box_right = grid_x * spacing - (grid_size * spacing / 2) + box_size / 1.5
        box_bottom = grid_y * spacing - ((grid_size * spacing / 2)) - box_size / 1.5
        box_top = grid_y * spacing - ((grid_size * spacing / 2)) + box_size / 1.5
        
        # Draw box
        rect = plt.Rectangle((box_left, box_bottom), box_right - box_left, box_top - box_bottom,
                            fill=False, edgecolor=community_colors[idx], linewidth=2, alpha=0.6)
        plt.gca().add_patch(rect)
        
        # Add community label
        color = community_colors[idx]
        if "INSR" in community:
            label = f"Community {idx+1} (INSR)"
        else:
            label = f"Community {idx+1}"
        plt.text(box_left + (box_right - box_left) / 2, box_bottom - 0.3,
                label, ha='center', va='top', fontsize=10, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor=color))
    
    plt.title("INSR Protein Interaction Network - Community Structure\nProteins grouped by community detection algorithm",
              fontsize=16, fontweight='bold', pad=20)
    plt.margins(0.2)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_fr_layout_ppi_visualization(G, data, output_path):
    """Create a visualization using Fruchterman-Reingold layout which is
    particularly good for PPI networks with proper node labeling.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Use Fruchterman-Reingold layout which works well for biological networks
    pos = nx.fruchterman_reingold_layout(G, k=0.5, seed=42)
    
    # Make sure INSR is centered
    if "INSR" in G.nodes():
        center = np.array([0, 0])
        offset = center - pos["INSR"]
        # Shift all positions
        for node in pos:
            pos[node] = pos[node] + offset
    
    # Assign node sizes based on degree centrality
    degree_centrality = nx.degree_centrality(G)
    node_sizes = []
    for node in G.nodes():
        if node == "INSR":
            node_sizes.append(5000)  # Larger size for INSR
        else:
            # Make sizes proportional to degree centrality
            node_sizes.append(1000 + 4000 * degree_centrality[node])
    
    # Node colors based on a pleasing color scheme
    # Assign different colors to nodes based on their degree
    degree_values = [d for n, d in G.degree()]
    max_degree = max(degree_values) if degree_values else 1
    
    node_colors = []
    cmap = plt.cm.plasma  # A sequential, perceptually uniform colormap
    
    for node in G.nodes():
        if node == "INSR":
            node_colors.append('#FFD700')  # Gold for INSR
        else:
            # Color based on degree
            degree = G.degree(node)
            node_colors.append(cmap(degree / max_degree))
    
    # Edge properties
    edge_colors = []
    edge_widths = []
    
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                score = float(interaction['score']) / 1000
                # Make edge widths proportional to confidence score
                edge_widths.append(0.5 + 2.5 * score)
                break
    
    # Draw network
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.6)
    
    # Draw nodes with alpha for better visibility
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5, alpha=0.85)
    
    # Draw all labels with varying font sizes based on node importance
    # This creates a more informative visualization
    labels = {}
    font_sizes = {}
    
    for node in G.nodes():
        labels[node] = node
        if node == "INSR":
            font_sizes[node] = 16  # Largest font for INSR
        elif degree_centrality[node] > 0.5:
            font_sizes[node] = 12  # Larger font for important nodes
        elif degree_centrality[node] > 0.2:
            font_sizes[node] = 10  # Medium font for moderately important nodes
        else:
            font_sizes[node] = 8   # Small font for less important nodes
    
    # Draw labels with varying sizes and white background for readability
    for node, label in labels.items():
        plt.text(pos[node][0], pos[node][1], label,
                fontsize=font_sizes[node], ha='center', va='center',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
    
    # Special highlight for INSR
    if "INSR" in G.nodes():
        plt.scatter(pos["INSR"][0], pos["INSR"][1], s=300, color='red', marker='*', alpha=1.0, zorder=10)
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='red', 
                  markersize=15, label='INSR'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(0.2), 
                  markersize=10, label='Low Degree'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(0.5), 
                  markersize=10, label='Medium Degree'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(0.8), 
                  markersize=10, label='High Degree'),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', 
              title='Network Elements', title_fontsize=12, fontsize=10,
              bbox_to_anchor=(1.1, 1.0))
    
    # Add title
    plt.title("INSR Protein Interaction Network (FR Layout)\nNode size: Degree centrality | Node color: Degree | Edge width: Confidence score",
              fontsize=16, fontweight='bold', pad=20)
    
    # Set plot parameters
    plt.axis('off')
    plt.tight_layout()
    
    # Save visualization
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_non_radial_ppi_visualization(G, data, output_path):
    """Create a non-radial PPI network visualization with natural clustering.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Create a better, less radial layout
    # 1. Start with a random layout to break any circular patterns
    pos = nx.random_layout(G, seed=42)
    
    # 2. Apply FR layout with careful parameter tuning to prevent radial formation
    # Lower k values prevent the radial effect
    pos = nx.fruchterman_reingold_layout(G, k=0.15, pos=pos, iterations=50, seed=42)
    
    # 3. Apply spring layout as a final refinement with low k to keep nodes close
    pos = nx.spring_layout(G, pos=pos, k=0.1, iterations=30, seed=42)
    
    # 4. Apply a slight randomization to break any remaining symmetry
    for node in pos:
        pos[node][0] += np.random.normal(0, 0.02)
        pos[node][1] += np.random.normal(0, 0.02)
    
    # 5. Ensure INSR is not exactly at center to break radial pattern
    if "INSR" in G.nodes():
        pos["INSR"] = np.array([0.1, 0.1])  # Slightly off-center
    
    # Create node groupings based on neighbors to enhance natural clusters
    # Group nodes that share many neighbors
    groups = {}
    for node in G.nodes():
        if node == "INSR":
            groups[node] = 0  # Special group for INSR
        else:
            # Use neighborhood similarity to determine groups
            neighbors = set(G.neighbors(node))
            if "INSR" in neighbors:
                # Nodes directly connected to INSR
                groups[node] = 1
            else:
                # Count common neighbors with INSR's neighborhood
                insr_neighbors = set(G.neighbors("INSR")) if "INSR" in G else set()
                common = len(neighbors.intersection(insr_neighbors))
                if common > 0:
                    # Nodes sharing neighbors with INSR
                    groups[node] = 2
                else:
                    # Other nodes
                    groups[node] = 3
    
    # Use a color palette that provides clear distinction
    group_colors = {
        0: '#FFD700',  # Gold for INSR
        1: '#FF6B6B',  # Red for direct connections
        2: '#4ECDC4',  # Teal for secondary connections
        3: '#45B7D1'   # Blue for others
    }
    
    # Assign node colors based on groups
    node_colors = [group_colors[groups[node]] for node in G.nodes()]
    
    # Vary node sizes based on degree and betweenness for visual hierarchy
    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)
    
    node_sizes = []
    for node in G.nodes():
        if node == "INSR":
            node_sizes.append(5000)  # Largest size for INSR
        else:
            size = 800 + 25 * degree[node] + 3000 * betweenness[node]
            node_sizes.append(size)
    
    # Get edge properties
    edge_colors = []
    edge_widths = []
    
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                score = float(interaction['score']) / 1000
                edge_widths.append(0.5 + 2.5 * score)
                break
    
    # Draw network with specific settings to enhance cluster visibility
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, 
                          alpha=0.5, connectionstyle='arc3,rad=0.1')
    
    # Draw nodes with distinct styles
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5, alpha=0.9)
    
    # Draw labels smartly - vary font size and only show important labels
    labels = {}
    important_nodes = []
    
    # Select important nodes to label
    # Always include INSR
    if "INSR" in G.nodes():
        important_nodes.append("INSR")
    
    # Add high-degree nodes (hubs)
    important_nodes.extend([node for node, d in sorted(degree.items(), 
                                                      key=lambda x: x[1], 
                                                      reverse=True)[:15] 
                           if node != "INSR"])
    
    # Add high-betweenness nodes (bridges)
    important_nodes.extend([node for node, b in sorted(betweenness.items(), 
                                                      key=lambda x: x[1], 
                                                      reverse=True)[:10]
                           if node != "INSR" and node not in important_nodes])
    
    # Create labels dictionary for important nodes
    labels = {node: node for node in important_nodes}
    
    # Draw labels with white background for better visibility
    for node, label in labels.items():
        font_size = 16 if node == "INSR" else 10
        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
        plt.text(pos[node][0], pos[node][1] + 0.02, label, 
                horizontalalignment='center', fontsize=font_size, fontweight='bold',
                bbox=bbox_props)
    
    # Add legend for node groups
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FFD700', 
                  markersize=15, label='INSR'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6B6B', 
                  markersize=10, label='Direct INSR Interactors'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#4ECDC4', 
                  markersize=10, label='Secondary Interactors'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#45B7D1', 
                  markersize=10, label='Other Proteins'),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', 
              title='Network Elements', title_fontsize=12, fontsize=10,
              bbox_to_anchor=(1.1, 1.0))
    
    # Add title and adjust plot parameters
    plt.title("INSR Protein Interaction Network - Natural Clusters\nNode color: Interaction proximity | Node size: Importance | Edge width: Confidence",
              fontsize=16, fontweight='bold', pad=20)
    
    plt.axis('off')
    plt.tight_layout()
    
    # Save visualization
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_pathway_enrichment_network_visualization(G, data, enrichment_data, output_path):
    """Create visualization focusing on biological pathways and relationships.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        enrichment_data (list): List of enrichment data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Process enrichment data to get pathway assignments
    pathway_colors = {}
    pathway_assignment = {}
    
    if enrichment_data and isinstance(enrichment_data, list):
        # Sort enrichment terms by p-value
        sorted_terms = sorted(enrichment_data, key=lambda x: float(x.get('pvalue', 1)))
        
        # Get top pathways (excluding generic terms)
        pathways = []
        for term in sorted_terms:
            if isinstance(term, dict) and 'description' in term:
                desc = term['description'].lower()
                # Skip generic terms
                if any(generic in desc for generic in ['regulation', 'process', 'activity', 'function']):
                    continue
                pathways.append(term['description'])
                if len(pathways) >= 6:  # Limit to top 6 pathways
                    break
        
        # Generate colors for pathways
        colors = plt.cm.Set3(np.linspace(0, 1, len(pathways)))
        pathway_colors = dict(zip(pathways, colors))
        
        # Assign proteins to pathways based on enrichment data
        for term in sorted_terms:
            if isinstance(term, dict) and 'description' in term and 'genes' in term:
                pathway = term['description']
                if pathway in pathways:
                    for gene in term['genes']:
                        if gene in G.nodes():
                            pathway_assignment[gene] = pathway
    
    # Assign remaining proteins to 'Other' category
    for node in G.nodes():
        if node not in pathway_assignment:
            pathway_assignment[node] = 'Other'
    
    # Add 'Other' to pathway colors if needed
    if 'Other' not in pathway_colors:
        pathway_colors['Other'] = '#7F7F7F'
    
    # Layout
    pos = nx.spring_layout(G, k=3, iterations=200, scale=2.0, seed=42)
    pos["INSR"] = (0, 0)
    
    # Node sizes based on betweenness centrality
    betweenness = nx.betweenness_centrality(G)
    node_sizes = [8000 if node == "INSR" else 1000 + 4000 * betweenness[node] for node in G.nodes()]
    
    # Node colors based on pathways
    node_colors = []
    for node in G.nodes():
        if node == "INSR":
            node_colors.append('#FFD700')
        else:
            pathway = pathway_assignment[node]
            node_colors.append(pathway_colors[pathway])
    
    # Edge properties
    edge_colors = []
    edge_widths = []
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                edge_widths.append(1 + 3 * float(interaction['score']) / 1000)
                break
    
    # Draw network
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.5)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, edgecolors='white', linewidths=2)
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
    
    # Special INSR label
    plt.annotate("INSR (Central Protein)", xy=pos["INSR"], xytext=(0, 0.1),
                ha='center', va='bottom', fontsize=12, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='#FFD700', linewidth=2))
    
    # Legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR (Central Protein)',
                  markerfacecolor='#FFD700', markersize=15)
    ]
    
    # Add pathway legend elements
    for pathway, color in pathway_colors.items():
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', label=pathway,
                      markerfacecolor=color, markersize=10)
        )
    
    # Add interaction type legend elements
    legend_elements.extend([
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3),
        plt.Line2D([0], [0], color='#90A4AE', label='Other Evidence', linewidth=3)
    ])
    
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
              title='Network Elements', title_fontsize=12, fontsize=10)
    
    # Add enrichment information if available
    if enrichment_data and isinstance(enrichment_data, list) and len(enrichment_data) > 0:
        enrichment_text = "Functional Enrichment:\n"
        for term in enrichment_data[:3]:
            if isinstance(term, dict) and 'description' in term and 'pvalue' in term:
                enrichment_text += f"{term['description']} (p={float(term['pvalue']):.2e})\n"
        
        plt.figtext(0.02, 0.02, enrichment_text,
                   fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.8))
    
    plt.title("INSR Protein Interaction Network - Biological Pathways\nNode size: Betweenness centrality | Edge width: Confidence score",
              fontsize=16, fontweight='bold', pad=20)
    plt.margins(0.3)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_spread_out_centered_visualization(G, data, output_path):
    """Create a spread-out visualization with INSR at center and nodes more spread out.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Modified layout algorithm for better spread
    # 1. Place INSR at center initially
    pos = {}
    pos["INSR"] = np.array([0, 0])
    
    # 2. Initialize other nodes in a uniform circular pattern but far from center
    angle = 0
    angle_step = 2 * np.pi / (G.number_of_nodes() - 1) if G.number_of_nodes() > 1 else 0
    radius = 2.0  # Larger radius for initial spread
    
    for node in G.nodes():
        if node != "INSR":
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pos[node] = np.array([x, y])
            angle += angle_step
    
    # 3. Run Fruchterman-Reingold with increased repulsion
    pos = nx.fruchterman_reingold_layout(G, pos=pos, k=5, iterations=200, weight=None, seed=42)
    
    # 4. Lock INSR at center
    pos["INSR"] = np.array([0, 0])
    
    # 5. Apply scaled spring layout for final adjustments but with INSR fixed
    fixed_nodes = ["INSR"]
    fixed_pos = {node: pos[node] for node in fixed_nodes}
    # Use a compatible way to fix positions
    for i in range(100):  # Manual iterations replacing spring_layout with fixed nodes
        pos_before = pos.copy()
        # Run one iteration of spring layout
        pos = nx.spring_layout(G, pos=pos, k=1.5, iterations=1, scale=3.0)
        # Force INSR to stay at original position
        pos["INSR"] = np.array([0, 0])
        
        # Check if positions have converged
        max_diff = max(np.linalg.norm(pos[n] - pos_before[n]) for n in G.nodes())
        if max_diff < 0.01:
            break
    
    # Node sizes based on betweenness centrality for visual hierarchy
    betweenness = nx.betweenness_centrality(G)
    node_sizes = []
    for node in G.nodes():
        if node == "INSR":
            node_sizes.append(5000)  # Largest size for INSR
        else:
            size = 1000 + 4000 * betweenness[node]
            node_sizes.append(size)
    
    # Node colors using a vibrant colormap for better contrast
    degree = dict(G.degree())
    max_degree = max(degree.values()) if degree.values() else 1
    
    node_colors = []
    cmap = plt.cm.viridis  # Using viridis colormap for better contrast
    
    for node in G.nodes():
        if node == "INSR":
            node_colors.append('#FFD700')  # Gold for INSR
        else:
            # Color based on degree centrality
            node_colors.append(cmap(degree[node] / max_degree))
    
    # Edge properties with adjusted widths for better visibility
    edge_colors = []
    edge_widths = []
    
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                # Increase edge width multiplier
                edge_widths.append(0.8 + 3.5 * float(interaction['score']) / 1000)
                break
    
    # Draw network with straight edges for cleaner look
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=2, alpha=0.9)
    
    # Draw labels for important nodes only to reduce clutter
    # Sort nodes by importance
    important_nodes = sorted([(node, betweenness[node] * degree[node]) 
                              for node in G.nodes() if node != "INSR"], 
                             key=lambda x: x[1], reverse=True)[:15]
    important_nodes = [node[0] for node in important_nodes]
    important_nodes.append("INSR")  # Always include INSR
    
    # Create labels dictionary
    labels = {node: node for node in important_nodes}
    
    # Draw labels with adjusted positions to reduce overlap
    label_pos = {}
    for node in labels:
        # Slightly adjust label position based on quadrant
        x, y = pos[node]
        offset_x = 0
        offset_y = 0
        
        if x >= 0 and y >= 0:  # first quadrant
            offset_x, offset_y = 0.03, 0.03
        elif x < 0 and y >= 0:  # second quadrant
            offset_x, offset_y = -0.03, 0.03
        elif x < 0 and y < 0:  # third quadrant
            offset_x, offset_y = -0.03, -0.03
        else:  # fourth quadrant
            offset_x, offset_y = 0.03, -0.03
        
        label_pos[node] = (x + offset_x, y + offset_y)
    
    # Draw the labels with white background for better visibility
    for node, label in labels.items():
        font_size = 12 if node == "INSR" else 10
        x, y = label_pos[node]
        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
        plt.text(x, y, label, 
                horizontalalignment='center', fontsize=font_size, fontweight='bold',
                bbox=bbox_props)
    
    # Add a nice legend with thicker lines
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', label='High Degree',
                  markerfacecolor=cmap(0.8), markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Medium Degree',
                  markerfacecolor=cmap(0.5), markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Low Degree',
                  markerfacecolor=cmap(0.2), markersize=10),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3),
        plt.Line2D([0], [0], color='#90A4AE', label='Other Evidence', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', 
              title='Network Elements', title_fontsize=12, fontsize=10,
              bbox_to_anchor=(1.1, 1.0))
    
    # Add informative title
    plt.title("INSR Protein Interaction Network - Spread-out View\nNode size: Importance | Node color: Degree | Edge width: Confidence",
              fontsize=16, fontweight='bold', pad=20)
    
    # Increase margins for better spacing
    plt.margins(0.2)
    plt.axis('off')
    plt.tight_layout()
    
    # Save visualization
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_functional_group_visualization(G, data, output_path):
    """Create a visualization that groups proteins by their function into a grid layout.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(22, 22))
    
    # Group nodes by function - using a simplified approach based on STRING interactions
    # 1. Get all interaction types for each protein
    protein_interactions = {}
    for node in G.nodes():
        protein_interactions[node] = []
    
    for interaction in data:
        p1 = interaction['preferredName_A']
        p2 = interaction['preferredName_B']
        int_type, _ = get_interaction_type(interaction)
        
        if p1 in protein_interactions:
            protein_interactions[p1].append(int_type)
        if p2 in protein_interactions:
            protein_interactions[p2].append(int_type)
    
    # 2. Determine the dominant interaction type for each protein
    functional_groups = {}
    for protein, int_types in protein_interactions.items():
        if not int_types:
            dominant_type = 'unknown'
        else:
            type_counts = {}
            for int_type in int_types:
                type_counts[int_type] = type_counts.get(int_type, 0) + 1
            
            # Find the most common type
            dominant_type = max(type_counts.items(), key=lambda x: x[1])[0]
        
        # Special case for INSR
        if protein == "INSR":
            dominant_type = "central"
        
        # Create or append to group
        if dominant_type not in functional_groups:
            functional_groups[dominant_type] = []
        functional_groups[dominant_type].append(protein)
    
    # 3. Define group colors
    group_colors = {
        'central': '#FFD700',       # Gold for INSR
        'physical': '#7E57C2',      # Purple
        'functional': '#4CAF50',    # Green
        'text_mining': '#F44336',   # Red
        'other': '#90A4AE',         # Gray
        'unknown': '#607D8B'        # Dark gray
    }
    
    # 4. Create a grid layout with groups
    # Calculate grid dimensions
    n_groups = len(functional_groups)
    grid_size = int(np.ceil(np.sqrt(n_groups)))
    
    # Initialize positions
    pos = {}
    box_centers = {}
    box_sizes = {}
    
    # First, assign a grid position to each group
    for idx, (group_type, members) in enumerate(functional_groups.items()):
        grid_x = idx % grid_size
        grid_y = idx // grid_size
        
        # Calculate center of the group's box
        center_x = grid_x * 5 - (grid_size * 2.5) + 2.5
        center_y = grid_y * 5 - (grid_size * 2.5) + 2.5
        
        box_centers[group_type] = (center_x, center_y)
        
        # Use a circular layout within each group box
        n_members = len(members)
        box_size = 1.8  # Adjust size based on members
        box_sizes[group_type] = box_size
        
        # Special case for INSR group
        if group_type == 'central':
            # Place INSR at center of its box
            pos["INSR"] = (center_x, center_y)
            continue
            
        # Arrange other group members in a circle
        for i, protein in enumerate(members):
            if n_members == 1:
                angle = 0
            else:
                angle = 2 * np.pi * i / n_members
                
            x = center_x + box_size * 0.7 * np.cos(angle)
            y = center_y + box_size * 0.7 * np.sin(angle)
            pos[protein] = (x, y)
    
    # 5. Apply a mild force-directed layout to refine positions within groups
    # First, determine which nodes to keep fixed
    fixed_nodes = ["INSR"]  # Keep INSR fixed
    fixed_pos = {node: pos[node] for node in fixed_nodes}
    
    # Now refine the layout for each group separately
    for group_type, members in functional_groups.items():
        if group_type == 'central':
            continue  # Skip INSR group
            
        # Create a subgraph for this group
        subgraph = G.subgraph(members)
        
        # Get the current positions for this group
        subgraph_pos = {node: pos[node] for node in members}
        
        # Apply a mild force-directed layout within the constraint of the group box
        center_x, center_y = box_centers[group_type]
        box_size = box_sizes[group_type]
        
        # Refine with a constrained spring layout
        subgraph_pos = nx.spring_layout(subgraph, pos=subgraph_pos, fixed=None, 
                                       iterations=50, k=1.0, seed=42)
        
        # Scale positions to fit within the group box
        for node, (x, y) in subgraph_pos.items():
            # Scale to box size
            x = center_x + x * box_size * 0.7
            y = center_y + y * box_size * 0.7
            pos[node] = (x, y)
    
    # 6. Node properties
    node_colors = []
    node_sizes = []
    
    # Calculate centrality for node sizing
    betweenness = nx.betweenness_centrality(G)
    degree = dict(G.degree())
    
    for node in G.nodes():
        # Determine node color based on its functional group
        for group_type, members in functional_groups.items():
            if node in members:
                node_colors.append(group_colors[group_type])
                break
        
        # Node sizes based on importance
        if node == "INSR":
            node_sizes.append(6000)  # Largest for INSR
        else:
            # Combine degree and betweenness for size
            importance = 800 + 30 * degree[node] + 3000 * betweenness[node]
            node_sizes.append(importance)
    
    # 7. Edge properties
    edge_colors = []
    edge_widths = []
    edge_styles = []
    
    for u, v in G.edges():
        u_group = None
        v_group = None
        
        # Find the groups of the nodes
        for group_type, members in functional_groups.items():
            if u in members:
                u_group = group_type
            if v in members:
                v_group = group_type
            if u_group and v_group:
                break
        
        # Determine edge style based on whether it's intra- or inter-group
        if u_group == v_group:
            # Intra-group edges are solid
            edge_styles.append('solid')
        else:
            # Inter-group edges are dashed
            edge_styles.append('dashed')
        
        # Color and width from interaction data
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                score = float(interaction['score']) / 1000
                edge_widths.append(0.7 + 2.5 * score)
                break
    
    # 8. Draw network
    # Draw edges with different styles for intra- vs inter-group
    for i, (u, v) in enumerate(G.edges()):
        # Draw edges with appropriate style
        style = edge_styles[i]
        if style == 'dashed':
            plt.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], 
                    linestyle='--', color=edge_colors[i], 
                    linewidth=edge_widths[i], alpha=0.6)
        else:
            plt.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], 
                    linestyle='-', color=edge_colors[i], 
                    linewidth=edge_widths[i], alpha=0.7)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5, alpha=0.9)
    
    # 9. Draw group boxes and labels
    for group_type, (center_x, center_y) in box_centers.items():
        box_size = box_sizes.get(group_type, 2.0)
        color = group_colors[group_type]
        
        # Draw a rectangle for each group
        rect = plt.Rectangle((center_x - box_size, center_y - box_size), 
                            2 * box_size, 2 * box_size,
                            fill=False, edgecolor=color, linewidth=3, alpha=0.7)
        plt.gca().add_patch(rect)
        
        # Add group label
        group_label = group_type.replace('_', ' ').title()
        if group_type == 'central':
            group_label = 'INSR (Central Protein)'
            
        plt.text(center_x, center_y - box_size - 0.3, group_label,
                ha='center', va='top', fontsize=12, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor=color, boxstyle='round,pad=0.3'))
    
    # 10. Draw important node labels
    important_nodes = []
    
    # Add INSR
    if "INSR" in G.nodes():
        important_nodes.append("INSR")
    
    # Add top nodes by importance
    node_importance = [(node, degree[node] * betweenness[node]) for node in G.nodes() if node != "INSR"]
    important_nodes.extend([node for node, _ in sorted(node_importance, key=lambda x: x[1], reverse=True)[:20]])
    
    # Draw labels for important nodes
    labels = {node: node for node in important_nodes}
    for node, label in labels.items():
        font_size = 14 if node == "INSR" else 9
        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
        plt.text(pos[node][0], pos[node][1] + 0.1, label, 
                horizontalalignment='center', fontsize=font_size, fontweight='bold',
                bbox=bbox_props)
    
    # 11. Add legend
    legend_elements = []
    
    # Add group colors to legend
    for group_type, color in group_colors.items():
        # Skip unknown if not used
        if group_type == 'unknown' and 'unknown' not in functional_groups:
            continue
            
        label = group_type.replace('_', ' ').title()
        if group_type == 'central':
            label = 'INSR (Central Protein)'
            
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor=color, markersize=10, label=label)
        )
    
    # Add edge styles
    legend_elements.extend([
        plt.Line2D([0], [0], color='gray', label='Intra-group Connection', linewidth=3),
        plt.Line2D([0], [0], color='gray', linestyle='--', label='Inter-group Connection', linewidth=3)
    ])
    
    plt.legend(handles=legend_elements, loc='upper right', 
              title='Functional Groups', title_fontsize=12, fontsize=10,
              bbox_to_anchor=(1.1, 1.0))
    
    # 12. Add title
    plt.title("INSR Protein Interaction Network - Functional Groups\nProteins grouped by dominant interaction type",
              fontsize=16, fontweight='bold', pad=20)
    
    # 13. Final plot settings
    plt.margins(0.2)
    plt.axis('off')
    plt.tight_layout()
    
    # Save visualization
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_circular_layout_visualization(G, data, output_path):
    """Create a circular layout visualization with nodes arranged in a circle.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Create circular layout
    pos = nx.circular_layout(G, scale=8)
    
    # Find and place INSR at the top of the circle
    insr_idx = list(G.nodes()).index("INSR") if "INSR" in G.nodes() else 0
    nodes = list(G.nodes())
    n = len(nodes)
    
    # Rearrange nodes to put INSR at the top
    new_pos = {}
    angle_shift = -np.pi/2  # Start at the top (270 degrees)
    for i, node in enumerate(nodes):
        angle = 2 * np.pi * ((i - insr_idx) % n) / n + angle_shift
        new_pos[node] = np.array([8 * np.cos(angle), 8 * np.sin(angle)])
    
    pos = new_pos
    
    # Node sizes based on degree centrality
    degree_centrality = nx.degree_centrality(G)
    node_sizes = []
    for node in G.nodes():
        if node == "INSR":
            node_sizes.append(3000)  # Larger size for INSR
        else:
            node_sizes.append(1000 + 2000 * degree_centrality[node])
    
    # Node colors based on clustering coefficient
    clustering = nx.clustering(G)
    node_colors = []
    for node in G.nodes():
        if node == "INSR":
            node_colors.append('#FFD700')  # Gold for INSR
        else:
            node_colors.append(plt.cm.viridis(clustering[node]))
    
    # Edge properties
    edge_colors = []
    edge_widths = []
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                edge_widths.append(1 + 2 * float(interaction['score']) / 1000)
                break
    
    # Draw network with curved edges for better visibility in circular layout
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, 
                          alpha=0.7, connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5)
    
    # Add labels with adjusted font sizes
    labels = {}
    for node in G.nodes():
        labels[node] = node
    
    for node, label in labels.items():
        size = 12 if node == "INSR" else 8
        plt.text(pos[node][0], pos[node][1], label,
                ha='center', va='center', fontsize=size, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR (Central Protein)',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3),
        plt.Line2D([0], [0], color='#90A4AE', label='Other Evidence', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='center',
              title='Network Elements', title_fontsize=12, fontsize=10)
    
    plt.title("INSR Protein Interaction Network - Circular Layout\nNode size: Degree centrality | Node color: Clustering coefficient",
              fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_hierarchical_layout_visualization(G, data, output_path):
    """Create a hierarchical layout visualization.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Create a directed version of the graph for hierarchical layout
    DG = nx.DiGraph()
    
    # Add INSR as root node
    for edge in G.edges():
        if "INSR" in edge:
            # Make all edges from INSR outgoing
            if edge[0] == "INSR":
                DG.add_edge(edge[0], edge[1])
            else:
                DG.add_edge(edge[1], edge[0])
        else:
            # For all other edges, maintain the original direction
            DG.add_edge(edge[0], edge[1])
    
    # Add any isolated nodes
    for node in G.nodes():
        if node not in DG:
            DG.add_node(node)
    
    # Compute node levels (distance from INSR)
    levels = {}
    for node in DG.nodes():
        if node == "INSR":
            levels[node] = 0
        elif nx.has_path(DG, "INSR", node):
            levels[node] = nx.shortest_path_length(DG, "INSR", node)
        else:
            # If no path from INSR, place at level 1 but to the side
            levels[node] = 1
    
    # Get max level for vertical spacing
    max_level = max(levels.values()) + 1
    
    # Create a hierarchical layout manually
    pos = {}
    level_counts = {}
    level_current = {}
    
    # Count nodes at each level
    for node, level in levels.items():
        if level not in level_counts:
            level_counts[level] = 0
            level_current[level] = 0
        level_counts[level] += 1
    
    # Compute positions
    for node, level in sorted(levels.items(), key=lambda x: (x[1], x[0])):
        if node == "INSR":
            # Place INSR at the top center
            pos[node] = np.array([0, 10])
        else:
            # Horizontal position based on count at this level
            width = max(15, level_counts[level] * 1.5)
            if level_counts[level] > 1:
                x = width * (level_current[level] / (level_counts[level] - 1) - 0.5)
            else:
                x = 0
            
            # Vertical position based on level
            y = 10 - 2 * level
            
            pos[node] = np.array([x, y])
            level_current[level] += 1
    
    # Node sizes based on betweenness centrality for visual hierarchy
    betweenness = nx.betweenness_centrality(G)
    node_sizes = []
    for node in DG.nodes():
        if node == "INSR":
            node_sizes.append(5000)  # Largest size for INSR
        else:
            size = 1000 + 3000 * betweenness[node]
            node_sizes.append(size)
    
    # Node colors based on hierarchical level
    node_colors = []
    for node in DG.nodes():
        if node == "INSR":
            node_colors.append('#FFD700')  # Gold for INSR
        else:
            # Color based on distance from INSR
            level = levels[node]
            node_colors.append(plt.cm.cool(level / 4))  # Max level is likely 3 or 4
    
    # Edge properties
    edge_colors = []
    edge_widths = []
    for u, v in DG.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                edge_widths.append(1 + 2 * float(interaction['score']) / 1000)
                break
    
    # Draw network
    nx.draw_networkx_edges(DG, pos, edge_color=edge_colors, width=edge_widths, 
                          alpha=0.7, arrows=True)
    nx.draw_networkx_nodes(DG, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5)
    
    # Add labels with optimal font sizes
    for node in DG.nodes():
        size = 16 if node == "INSR" else 10
        plt.text(pos[node][0], pos[node][1], node, ha='center', va='center', 
                fontsize=size, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR (Root Protein)',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', label='Level 1 Proteins',
                  markerfacecolor=plt.cm.cool(1/4), markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Level 2 Proteins',
                  markerfacecolor=plt.cm.cool(2/4), markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Level 3+ Proteins',
                  markerfacecolor=plt.cm.cool(3/4), markersize=10),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right',
              title='Network Elements', title_fontsize=12, fontsize=10)
    
    plt.title("INSR Protein Interaction Network - Hierarchical Layout\nNode size: Betweenness centrality | Node color: Hierarchical level",
              fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_radial_layout_visualization(G, data, output_path):
    """Create a radial layout visualization with INSR at the center.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Create a shell (concentric circles) layout with INSR at center
    # First, sort other nodes by their importance (centrality)
    betweenness = nx.betweenness_centrality(G)
    sorted_nodes = sorted([(n, betweenness[n]) for n in G.nodes() if n != "INSR"], 
                        key=lambda x: x[1], reverse=True)
    
    # Create shells (concentric circles)
    shells = [["INSR"]]  # INSR in the innermost shell
    
    # Divide other nodes into 2 additional shells based on centrality
    n = len(sorted_nodes)
    cutoff1 = n // 3
    cutoff2 = 2 * n // 3
    
    shells.append([node for node, _ in sorted_nodes[:cutoff1]])
    shells.append([node for node, _ in sorted_nodes[cutoff1:cutoff2]])
    shells.append([node for node, _ in sorted_nodes[cutoff2:]])
    
    # Generate shell layout
    pos = nx.shell_layout(G, shells, scale=8)
    
    # Node sizes and colors
    degree = dict(G.degree())
    max_degree = max(degree.values())
    
    # Adjust node sizes based on degree
    node_sizes = []
    for node in G.nodes():
        if node == "INSR":
            node_sizes.append(5000)  # Largest for central node
        else:
            node_sizes.append(1000 + 1500 * degree[node] / max_degree)
    
    # Color nodes by their shell assignment for visual grouping
    node_colors = []
    for node in G.nodes():
        if node == "INSR":
            node_colors.append('#FFD700')  # Gold for INSR
        elif node in shells[1]:
            node_colors.append('#FF5722')  # First shell (high centrality)
        elif node in shells[2]:
            node_colors.append('#2196F3')  # Second shell (medium centrality)
        else:
            node_colors.append('#4CAF50')  # Third shell (low centrality)
    
    # Edge properties with gradient coloring based on distance from center
    edge_colors = []
    edge_widths = []
    
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                
                # Make edges wider based on score and if connected to INSR
                if "INSR" in (u, v):
                    edge_widths.append(2 + 2 * float(interaction['score']) / 1000)
                else:
                    edge_widths.append(0.5 + 2 * float(interaction['score']) / 1000)
                break
    
    # Draw network
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5)
    
    # Add labels with varying font sizes
    for node in G.nodes():
        # Size based on importance
        if node == "INSR":
            size = 16
        elif G.degree(node) > G.number_of_nodes() / 4:
            size = 12
        else:
            size = 8
        
        # Enhanced label with background
        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.7)
        plt.text(pos[node][0], pos[node][1], node, ha='center', va='center', 
                fontsize=size, fontweight='bold', bbox=bbox_props)
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR (Central Protein)',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', label='First Shell (High Centrality)',
                  markerfacecolor='#FF5722', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Second Shell (Medium Centrality)',
                  markerfacecolor='#2196F3', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Third Shell (Low Centrality)',
                  markerfacecolor='#4CAF50', markersize=10),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right',
              title='Network Elements', title_fontsize=12, fontsize=10)
    
    plt.title("INSR Protein Interaction Network - Radial Layout\nNode size: Degree | Node color: Centrality shell",
              fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_spectral_layout_visualization(G, data, output_path):
    """Create a spectral layout visualization based on graph Laplacian.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Create spectral layout
    pos = nx.spectral_layout(G, scale=10)
    
    # Ensure INSR is in a prominent position
    if "INSR" in G.nodes():
        offset = np.array([0, 0]) - pos["INSR"]
        for node in pos:
            pos[node] = pos[node] + offset
    
    # Node sizes based on eigenvector centrality
    try:
        eigenvector_centrality = nx.eigenvector_centrality(G)
    except:
        # Fallback to degree centrality if eigenvector centrality fails
        eigenvector_centrality = nx.degree_centrality(G)
    
    node_sizes = []
    for node in G.nodes():
        if node == "INSR":
            node_sizes.append(5000)  # Largest size for INSR
        else:
            node_sizes.append(1000 + 3000 * eigenvector_centrality[node])
    
    # Node colors based on community detection
    try:
        communities = nx.community.greedy_modularity_communities(G)
        node_to_community = {}
        for i, community in enumerate(communities):
            for node in community:
                node_to_community[node] = i
        
        # Create a colormap
        community_colors = plt.cm.tab10(np.linspace(0, 1, len(communities)))
        
        node_colors = []
        for node in G.nodes():
            if node == "INSR":
                node_colors.append('#FFD700')  # Gold for INSR
            else:
                comm_idx = node_to_community.get(node, 0)
                node_colors.append(community_colors[comm_idx])
    except:
        # Fallback coloring if community detection fails
        node_colors = ['#FFD700' if node == "INSR" else '#1F77B4' for node in G.nodes()]
    
    # Edge properties
    edge_colors = []
    edge_widths = []
    
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                edge_widths.append(1 + 2 * float(interaction['score']) / 1000)
                break
    
    # Draw network
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5)
    
    # Add labels with varying font sizes
    for node in G.nodes():
        size = 16 if node == "INSR" else 10
        plt.text(pos[node][0], pos[node][1], node, ha='center', va='center', 
                fontsize=size, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add legend for interaction types
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR (Central Protein)',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3),
        plt.Line2D([0], [0], color='#90A4AE', label='Other Evidence', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right',
              title='Network Elements', title_fontsize=12, fontsize=10)
    
    plt.title("INSR Protein Interaction Network - Spectral Layout\nNode size: Eigenvector centrality | Node color: Community",
              fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_force_atlas_visualization(G, data, output_path):
    """Create a custom force-atlas like layout, which is popular for PPI networks.
    
    Args:
        G (nx.Graph): NetworkX graph object
        data (list): List of interaction data
        output_path (str): Path to save the visualization
    """
    plt.figure(figsize=(20, 20))
    
    # Start with a random layout
    pos = nx.random_layout(G, seed=42)
    
    # Custom force-directed layout simulation
    repulsion = 0.2
    attraction = 0.05
    max_displacement = 0.5
    temperature = 0.5
    cooling_factor = 0.95
    
    # Force INSR to center
    pos["INSR"] = np.array([0, 0])
    
    for _ in range(100):  # Number of iterations
        # Calculate repulsive forces
        displacement = {node: np.array([0.0, 0.0]) for node in G.nodes()}
        
        # Repulsive forces between nodes
        for i, node1 in enumerate(G.nodes()):
            if node1 == "INSR":
                continue  # Skip INSR
                
            for j, node2 in enumerate(G.nodes()):
                if i != j and node2 != "INSR":
                    delta = pos[node1] - pos[node2]
                    distance = max(0.01, np.linalg.norm(delta))
                    
                    # Repulsive force inversely proportional to distance
                    force = repulsion / distance**2
                    displacement[node1] += delta * force / distance
        
        # Attractive forces along edges
        for edge in G.edges():
            u, v = edge
            if "INSR" in edge:
                continue  # Skip edges with INSR to keep it centered
                
            delta = pos[u] - pos[v]
            distance = max(0.01, np.linalg.norm(delta))
            
            # Attractive force proportional to distance
            force = attraction * distance
            if u != "INSR":
                displacement[u] -= delta * force / distance
            if v != "INSR":
                displacement[v] += delta * force / distance
        
        # Limit displacement by temperature and update positions
        for node in G.nodes():
            if node == "INSR":
                continue  # Keep INSR at center
                
            disp_length = max(0.01, np.linalg.norm(displacement[node]))
            pos[node] += displacement[node] * min(disp_length, temperature) / disp_length
        
        # Cool the system
        temperature *= cooling_factor
    
    # Scale the layout
    min_x = min(x for x, _ in pos.values())
    max_x = max(x for x, _ in pos.values())
    min_y = min(y for _, y in pos.values())
    max_y = max(y for _, y in pos.values())
    
    scale = 5.0 / max(max_x - min_x, max_y - min_y)
    for node in pos:
        pos[node] = pos[node] * scale
    
    # Node sizes based on degree centrality
    degree_centrality = nx.degree_centrality(G)
    node_sizes = []
    for node in G.nodes():
        if node == "INSR":
            node_sizes.append(6000)  # Larger size for INSR
        else:
            node_sizes.append(1000 + 4000 * degree_centrality[node])
    
    # Node colors based on PageRank centrality
    pagerank = nx.pagerank(G)
    node_colors = []
    for node in G.nodes():
        if node == "INSR":
            node_colors.append('#FFD700')  # Gold for INSR
        else:
            node_colors.append(plt.cm.plasma(pagerank[node] / max(pagerank.values())))
    
    # Edge properties with bundling effect
    edge_colors = []
    edge_widths = []
    
    for u, v in G.edges():
        for interaction in data:
            if (interaction['preferredName_A'] == u and interaction['preferredName_B'] == v) or \
               (interaction['preferredName_A'] == v and interaction['preferredName_B'] == u):
                _, color = get_interaction_type(interaction)
                edge_colors.append(color)
                
                # Make edges wider based on score and if connected to INSR
                if "INSR" in (u, v):
                    edge_widths.append(2 + 2 * float(interaction['score']) / 1000)
                else:
                    edge_widths.append(0.5 + 2 * float(interaction['score']) / 1000)
                break
    
    # Draw network with curved edges for a more organic look
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, 
                          alpha=0.5, connectionstyle='arc3,rad=0.1')
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          edgecolors='white', linewidths=1.5, alpha=0.9)
    
    # Add labels with varying font sizes
    for node in G.nodes():
        # Size based on importance
        if node == "INSR":
            size = 16
        elif G.degree(node) > G.number_of_nodes() / 4:
            size = 12
        else:
            size = 8
        
        # Enhanced label with background
        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.7)
        plt.text(pos[node][0], pos[node][1], node, ha='center', va='center', 
                fontsize=size, fontweight='bold', bbox=bbox_props)
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='INSR (Central Protein)',
                  markerfacecolor='#FFD700', markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', label='High PageRank',
                  markerfacecolor=plt.cm.plasma(0.8), markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Medium PageRank',
                  markerfacecolor=plt.cm.plasma(0.5), markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Low PageRank',
                  markerfacecolor=plt.cm.plasma(0.2), markersize=10),
        plt.Line2D([0], [0], color='#7E57C2', label='Physical Interaction', linewidth=3),
        plt.Line2D([0], [0], color='#4CAF50', label='Functional Association', linewidth=3),
        plt.Line2D([0], [0], color='#F44336', label='Text Mining', linewidth=3)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right',
              title='Network Elements', title_fontsize=12, fontsize=10)
    
    plt.title("INSR Protein Interaction Network - Force Atlas Layout\nNode size: Degree centrality | Node color: PageRank",
              fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close() 