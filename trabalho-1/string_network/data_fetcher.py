import requests

def get_string_data(protein, species=9606, network_type="all", required_score=400, limit=100):
    """Fetch protein interaction data from STRING database API.
    
    Args:
        protein (str): Protein identifier
        species (int): NCBI taxonomy ID (default: 9606 for human)
        network_type (str): Type of network to fetch
        required_score (int): Minimum confidence score (0-1000)
        limit (int): Maximum number of interactions to fetch
        
    Returns:
        list: List of interaction data dictionaries
    """
    string_api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": protein,
        "species": species,
        "network_type": network_type,
        "required_score": required_score,
        "limit": limit
    }
    response = requests.get(string_api_url, params=params)
    return response.json()

def get_enrichment_data(protein, species=9606):
    """Fetch functional enrichment data from STRING.
    
    Args:
        protein (str): Protein identifier
        species (int): NCBI taxonomy ID (default: 9606 for human)
        
    Returns:
        list: List of enrichment data dictionaries
    """
    enrichment_url = "https://string-db.org/api/json/enrichment"
    params = {
        "identifiers": protein,
        "species": species,
        "background_string_identifiers": "9606"
    }
    response = requests.get(enrichment_url, params=params)
    return response.json() 