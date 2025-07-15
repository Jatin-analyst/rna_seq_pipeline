import os
import GEOparse

def load_geo_data(accession_id_or_path, destdir="data", parse_platforms=False):
    """
    Load GEO Series Matrix data from a local .txt.gz file or download using accession ID.
    
    Parameters:
        accession_id_or_path (str): Either a GEO accession ID (e.g., 'GSE183947')
                                    or a path to a local .txt.gz file.
        destdir (str): Directory to save file if downloading.
        parse_platforms (bool): Whether to parse platform annotation.
        
    Returns:
        gse (GEOparse.GEO.GSE): Parsed GEO object
    """
    # Case 1: Local file path provided
    if os.path.exists(accession_id_or_path) and accession_id_or_path.endswith(".txt.gz"):
        print(f"Loading local file: {accession_id_or_path}")
        return GEOparse.get_GEO(filepath=accession_id_or_path, parse_platforms=parse_platforms)

    # Case 2: Assume GEO accession ID
    elif len(accession_id_or_path) >= 4 and accession_id_or_path.lower().startswith("gse"):
        print(f"Downloading GEO dataset: {accession_id_or_path}")
        return GEOparse.get_GEO(geo=accession_id_or_path, destdir=destdir, parse_platforms=parse_platforms)

    # Invalid input
    else:
        raise ValueError("Invalid input. Provide a valid GEO accession ID or a .txt.gz file path.")


gse = load_geo_data("GSE301959")

expr_df = gse.pivot_samples("Value")

