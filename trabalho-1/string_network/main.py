import pandas as pd
import os
from data_fetcher import get_string_data, get_enrichment_data
from network_analysis import create_network, analyze_network_metrics
from visualizations import (
    create_interaction_evidence_network_visualization,
    create_community_structure_network_visualization,
    create_fr_layout_ppi_visualization,
    create_non_radial_ppi_visualization,
    create_pathway_enrichment_network_visualization,
    create_functional_group_visualization,
    create_circular_layout_visualization,
    create_radial_layout_visualization,
    create_force_atlas_visualization
)

def main():
    # Ensure output directory exists
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Fetch protein interaction data
    protein = "INSR"
    species = 9606  # Human
    network_type = "functional"
    required_score = 400  # Medium confidence
    limit = 100
    
    interaction_data = get_string_data(protein, species, network_type, required_score, limit)
    print(f"Fetched {len(interaction_data)} interactions")
    
    # Optionally fetch enrichment data
    use_enrichment = True  # Set to True to include enrichment-based visualization
    enrichment_data = None
    if use_enrichment:
        enrichment_data = get_enrichment_data(protein, species)
        print(f"Fetched enrichment data with {len(enrichment_data)} terms")
    
    # Create network
    G = create_network(interaction_data)
    
    # Analyze network metrics
    metrics = analyze_network_metrics(G)
    print("\nNetwork Metrics:")
    for metric, value in metrics.items():
        print(f"{metric}: {value:.3f}")
    
    # Get top proteins by degree centrality
    degree = dict(G.degree())
    top_proteins = sorted(degree.items(), key=lambda x: x[1], reverse=True)[:10]
    print("\nTop proteins by degree:")
    for protein, degree_value in top_proteins:
        print(f"{protein}: {degree_value}")
    
    # Save network data
    network_data = pd.DataFrame(interaction_data)
    output_file = os.path.join(output_dir, "insr_network_edges.csv")
    network_data.to_csv(output_file, index=False)
    print(f"\nNetwork data saved to {output_file}")
    
    # Generate multiple visualizations with different layouts and approaches
    print("\nGenerating visualizations...")
    
    # 1. Evidence-based visualization
    evidence_viz_path = os.path.join(output_dir, "insr_network_evidence.png")
    create_interaction_evidence_network_visualization(G, interaction_data, evidence_viz_path)
    print(f"Evidence-based visualization saved to {evidence_viz_path}")
    
    # 2. Community structure visualization
    community_viz_path = os.path.join(output_dir, "insr_network_communities.png")
    create_community_structure_network_visualization(G, interaction_data, community_viz_path)
    print(f"Community-based visualization saved to {community_viz_path}")
    
    # 3. Fruchterman-Reingold layout visualization
    fr_viz_path = os.path.join(output_dir, "insr_network_fr_layout.png")
    create_fr_layout_ppi_visualization(G, interaction_data, fr_viz_path)
    print(f"FR layout visualization saved to {fr_viz_path}")
    
    # 4. Non-radial natural clustering visualization
    non_radial_viz_path = os.path.join(output_dir, "insr_network_natural_clusters.png")
    create_non_radial_ppi_visualization(G, interaction_data, non_radial_viz_path)
    print(f"Non-radial visualization saved to {non_radial_viz_path}")
    
    # 5. Functional group visualization
    functional_viz_path = os.path.join(output_dir, "insr_network_functional_groups.png")
    create_functional_group_visualization(G, interaction_data, functional_viz_path)
    print(f"Functional group visualization saved to {functional_viz_path}")
    
    # 6. Circular layout visualization
    circular_viz_path = os.path.join(output_dir, "insr_network_circular.png")
    create_circular_layout_visualization(G, interaction_data, circular_viz_path)
    print(f"Circular layout visualization saved to {circular_viz_path}")
    
    # 7. Radial layout visualization
    radial_viz_path = os.path.join(output_dir, "insr_network_radial.png")
    create_radial_layout_visualization(G, interaction_data, radial_viz_path)
    print(f"Radial layout visualization saved to {radial_viz_path}")
    
    # 8. Force Atlas layout visualization
    force_atlas_viz_path = os.path.join(output_dir, "insr_network_force_atlas.png")
    create_force_atlas_visualization(G, interaction_data, force_atlas_viz_path)
    print(f"Force Atlas layout visualization saved to {force_atlas_viz_path}")
    
    # 9. Pathway enrichment visualization (if enrichment data is available)
    if use_enrichment and enrichment_data:
        enrichment_viz_path = os.path.join(output_dir, "insr_network_pathways.png")
        create_pathway_enrichment_network_visualization(G, interaction_data, enrichment_data, enrichment_viz_path)
        print(f"Pathway enrichment visualization saved to {enrichment_viz_path}")
    
    print("\nVisualization files created successfully!")

if __name__ == "__main__":
    main() 