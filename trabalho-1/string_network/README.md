# STRING Network Visualization Tool

This tool fetches protein-protein interaction (PPI) data from the STRING database and generates various network visualizations to analyze the relationships between proteins.

## Visualizations

The tool creates the following visualizations in the `output` folder, each highlighting different aspects of the protein interaction network:

### 1. Evidence-based Visualization

**File:** `output/insr_network_evidence.png`

**Technical Implementation:**
- Uses a spring layout with INSR (central protein) fixed at the center
- Node sizes are proportional to betweenness centrality
- Edge colors represent different interaction types (physical, functional, text mining)
- Edge widths are proportional to confidence scores

**Data Interpretation:**
This visualization focuses on the evidence types and confidence scores of interactions. It helps identify which interactions are supported by physical evidence versus those derived from text mining or functional associations. The size of nodes indicates their importance as connectors between different parts of the network (betweenness centrality).

### 2. Community Structure Visualization

**File:** `output/insr_network_communities.png`

**Technical Implementation:**
- Uses the Louvain method to detect communities in the network
- Arranges communities in a grid layout
- Each community is color-coded and placed in a separate box
- Spring layout is applied within each community

**Data Interpretation:**
This visualization reveals how proteins naturally cluster into functional communities. Proteins within the same box/color are more densely connected to each other than to the rest of the network, suggesting they may be involved in similar biological processes or pathways.

### 3. Fruchterman-Reingold Layout Visualization

**File:** `output/insr_network_fr_layout.png`

**Technical Implementation:**
- Uses Fruchterman-Reingold algorithm, which is particularly well-suited for biological networks
- Node sizes reflect degree centrality
- Node colors use a sequential plasma colormap based on node degree
- Edge widths are proportional to confidence scores

**Data Interpretation:**
The FR layout creates a balanced distribution of nodes with minimal edge crossing, making it easier to visualize the overall network structure. This layout tends to position highly connected nodes (hubs) centrally and place nodes with similar connections close to each other, revealing the natural organization of the network.

### 4. Non-radial Natural Clustering Visualization

**File:** `output/insr_network_natural_clusters.png`

**Technical Implementation:**
- Starts with a random layout to break circular patterns
- Applies a customized force-directed layout to prevent radial formation
- Nodes are grouped by their relationship to INSR (direct connections, secondary connections, etc.)
- Node colors indicate these relationship groups

**Data Interpretation:**
This layout avoids the common "hairball" problem in network visualization by breaking the tendency for nodes to arrange radially around a central hub. It reveals natural clusters based on connection patterns while allowing different regions of the network to spread out more naturally, making it easier to identify functional subgroups.

### 5. Functional Group Visualization

**File:** `output/insr_network_functional_groups.png`

**Technical Implementation:**
- Groups proteins by their dominant interaction type
- Arranges groups in a grid layout
- Uses a custom force-directed layout within each functional group
- Intra-group edges are solid, inter-group edges are dashed

**Data Interpretation:**
This visualization organizes proteins based on their primary functional relationship types (physical interactions, functional associations, text mining evidence). This grouping helps identify which proteins share similar types of evidence and reveals how different functional groups interact with each other.

### 6. Circular Layout Visualization

**File:** `output/insr_network_circular.png`

**Technical Implementation:**
- Arranges nodes in a perfect circle
- Places INSR at the top position
- Node sizes reflect degree centrality
- Node colors represent clustering coefficient
- Uses curved edges for better visibility

**Data Interpretation:**
The circular layout provides a clean, organized view that's useful for visualizing the density and distribution of connections. Edge patterns become more apparent, and it's easier to compare connection patterns between different proteins. The clustering coefficient coloring reveals which proteins participate in tightly connected neighborhoods.

### 7. Radial Layout Visualization

**File:** `output/insr_network_radial.png`

**Technical Implementation:**
- Uses concentric shells based on betweenness centrality
- INSR is positioned at the center
- Other proteins are arranged in shells based on their importance
- Different shells are assigned different colors

**Data Interpretation:**
This layout organizes proteins in concentric circles with INSR at the center, creating a hierarchical view of the network. Proteins are placed in different shells based on their centrality, with the most important bridge nodes (high betweenness) in the inner shells. This reveals the hierarchy of the network and how information might flow through it.

### 8. Force Atlas Layout Visualization

**File:** `output/insr_network_force_atlas.png`

**Technical Implementation:**
- Implements a custom force-directed algorithm inspired by Force Atlas
- Includes customized repulsion and attraction forces
- Node colors are based on PageRank centrality
- Node sizes reflect degree centrality

**Data Interpretation:**
The Force Atlas layout is specifically designed for biological networks, creating a layout where related nodes cluster naturally while maintaining separation between different functional groups. This visualization highlights how proteins cluster based on their interaction patterns, with PageRank coloring revealing which proteins are most influential in the network structure.

### 9. Pathway Enrichment Visualization

**File:** `output/insr_network_pathways.png`

**Technical Implementation:**
- Uses enrichment data from STRING database
- Assigns nodes to specific pathways based on enrichment analysis
- Color codes nodes by their assigned pathway
- Includes enrichment statistics (p-values) in the visualization

**Data Interpretation:**
This visualization incorporates functional enrichment data to reveal biological pathway membership. Proteins are colored according to their involvement in specific biological pathways, allowing you to see how pathway members are distributed and interconnected in the network. This helps identify which pathways are most relevant to your protein of interest.

## Technical Details

- **Library:** Uses NetworkX for network creation and analysis
- **Visualization:** Matplotlib for rendering
- **Data Source:** STRING database API
- **Metrics Calculated:** Degree centrality, betweenness centrality, clustering coefficient, network diameter, and more

For protein interactions, the system uses the confidence score from STRING to determine the edge width, making stronger interactions visually more prominent.

## Output Files

All output files are stored in the `output` directory to keep the project root clean:

- **CSV files:**
  - `output/insr_network_edges.csv` - Raw interaction data from STRING database

- **Visualization files:**
  - All network visualization PNG files described above

You can modify the output directory in the `main.py` file by changing the `output_dir` variable. 