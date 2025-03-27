# Protein-Protein Interaction Network Analysis: INSR Case Study

## Overview
This project analyzes the protein interaction network of the Insulin Receptor (INSR) using data from the STRING database. The analysis focuses on understanding INSR's role in cellular processes through network visualization and metrics analysis.

## Data Collection and Types

### 1. STRING Database Data
We collect two main types of data from the STRING database:

#### a. Protein Interaction Data
- **Source**: STRING database API (`get_string_data`)
- **Target Protein**: INSR (Insulin Receptor)
- **Species**: Homo sapiens (NCBI taxonomy ID: 9606)
- **Network Type**: All interaction types
- **Confidence Threshold**: 0.4 (400 in STRING's 0-1000 scale)
- **Sample Size**: 100 interactions
- **Data Fields**:
  - `preferredName_A/B`: Protein identifiers
  - `score`: Overall interaction confidence (0-1000)
  - `escore`: Experimental evidence score
  - `dscore`: Database evidence score
  - `tscore`: Text mining evidence score

#### b. Functional Enrichment Data
- **Source**: STRING database API (`get_enrichment_data`)
- **Data Fields**:
  - `description`: Pathway/function description
  - `pvalue`: Statistical significance
  - `genes`: List of genes in the pathway
  - `term`: GO term or pathway identifier

### 2. Network Analysis Metrics
We implemented several network metrics to characterize the INSR interaction network:
- **Degree Distribution**: Variance of node degrees
- **Clustering Coefficient**: Average clustering coefficient
- **Characteristic Path Length**: Average shortest path length
- **Centrality Measures**: 
  - Degree centrality
  - Betweenness centrality

### 3. Visualization Strategies
We developed multiple visualization approaches to highlight different aspects of the network:

#### a. Interaction Evidence Network Visualization
- Force-directed layout with INSR at center
- Node size based on betweenness centrality
- Edge colors representing evidence types:
  - Purple: Physical interactions (experimental evidence)
  - Green: Functional associations (database evidence)
  - Red: Text mining evidence
  - Gray: Other evidence types
- Edge width based on confidence score

#### b. Community Structure Network Visualization
- Louvain community detection
- Group-in-a-box layout
- Color-coded communities
- Enhanced visibility of functional modules

#### c. Pathway Enrichment Network Visualization
- Pathway-based node coloring
- Evidence-based edge coloring
- Confidence-based edge widths
- Special highlighting of INSR
- Top enriched pathways displayed

## Results

### 1. Network Statistics
- Total interactions: 1,996
- Average degree: 39.52
- Degree distribution variance: 502.31
- Clustering coefficient: 0.776
- Characteristic path length: 1.605

### 2. Interaction Types Distribution
- Physical interactions: 126
- Functional associations: 242
- Text mining evidence: 265
- Other evidence types: 1,363

### 3. Evidence Distribution
- Both database and experimental evidence: 443
- Only database evidence: 741
- Only experimental evidence: 239
- Neither database nor experimental evidence: 573

### 4. Key Proteins
Top 5 proteins by degree centrality:
1. INSR (1.000)
2. IGF1R (0.900)
3. EGFR (0.890)
4. ERBB2 (0.790)
5. PIK3R1 (0.770)

Top 5 proteins by betweenness centrality:
1. INSR (0.098)
2. EGFR (0.054)
3. IGF1R (0.053)
4. ERBB2 (0.030)
5. SRC (0.026)

### 5. Biological Pathways
The network reveals several key biological pathways:
- Growth Factor Signaling
- Metabolic Regulation
- Cell Signaling
- Signal Transduction
- Cell Cycle
- Apoptosis

## Technical Implementation

### 1. Tools and Libraries
- Python 3.x
- NetworkX for network analysis
- Matplotlib for visualization
- Pandas for data handling
- Requests for API communication

### 2. Key Functions
- `degree_distribution()`: Calculates variance of node degrees
- `clustering_coefficient()`: Computes average clustering coefficient
- `characteristic_path_length()`: Determines average shortest path length
- `create_enhanced_biological_visualization()`: Generates pathway-based visualization

### 3. Visualization Features
- Force-directed layout with INSR at center
- Node sizes based on betweenness centrality
- Edge colors representing evidence types
- Edge widths indicating confidence scores
- Pathway-based node coloring
- Special highlighting of INSR

## Conclusions

### 1. Network Properties
- High clustering coefficient indicates strong modularity
- Short characteristic path length suggests efficient information flow
- High degree variance indicates presence of hub proteins

### 2. Biological Insights
- INSR's central role in multiple signaling pathways
- Strong connection to growth factor signaling
- Integration of metabolic and cell signaling networks
- High confidence in most interactions

### 3. Visualization Effectiveness
- Multiple visualization approaches provide different perspectives
- Community detection reveals functional modules
- Pathway-based coloring highlights biological relationships
- Special highlighting of INSR improves focus on key protein

## References
1. STRING database: https://string-db.org/
2. NetworkX documentation: https://networkx.org/
3. Matplotlib documentation: https://matplotlib.org/ 