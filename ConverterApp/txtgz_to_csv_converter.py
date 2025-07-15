#!/usr/bin/env python3
"""
Advanced TXT.GZ to TSV Converter for Genomic Data
==================================================

Converts .txt.gz files from GEO and other biological databases to clean
tab-separated files with genes as rows and samples as columns.

Features:
- Automatic metadata and footnote removal
- Smart gene and sample column detection
- Duplicate gene handling
- Tab-separated output format
- Batch processing support

Author: Claude AI Assistant
"""

import gzip
import re
import os
import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import pandas as pd
from collections import Counter


class GeneExpressionConverter:
    """
    Advanced converter for genomic expression data files.
    
    Handles messy .txt.gz files from GEO and similar databases,
    cleaning metadata and formatting as gene expression matrices.
    """
    
    def __init__(self, log_level: int = logging.INFO):
        """
        Initialize the converter.
        
        Args:
            log_level: Logging level (logging.INFO, logging.DEBUG, etc.)
        """
        self.setup_logging(log_level)
        self.reset_stats()
        
    def setup_logging(self, level: int) -> None:
        """Configure logging with file and console output."""
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler('gene_converter.log')
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def reset_stats(self) -> None:
        """Reset processing statistics."""
        self.stats = {
            'total_lines': 0,
            'metadata_lines': 0,
            'footnote_lines': 0,
            'inconsistent_lines': 0,
            'final_genes': 0,
            'final_samples': 0
        }
    
    # ========================
    # DATA DETECTION METHODS
    # ========================
    
    def detect_delimiter(self, lines: List[str]) -> str:
        """
        Detect the most common delimiter in the data.
        
        Args:
            lines: List of data lines to analyze
            
        Returns:
            Most likely delimiter character
        """
        delimiters = ['\t', ',', ';', '|']
        delimiter_counts = Counter()
        
        # Sample first 50 non-metadata lines
        sample_lines = []
        for line in lines[:100]:  # Check more lines for better detection
            if line.strip() and not self._is_metadata_line(line):
                sample_lines.append(line)
                if len(sample_lines) >= 50:
                    break
        
        # Count delimiter occurrences
        for line in sample_lines:
            for delim in delimiters:
                delimiter_counts[delim] += line.count(delim)
        
        if delimiter_counts:
            detected_delim = delimiter_counts.most_common(1)[0][0]
            self.logger.info(f"Detected delimiter: '{detected_delim}'")
            return detected_delim
        
        self.logger.warning("No delimiter detected, defaulting to tab")
        return '\t'
    
    def detect_header_row(self, lines: List[str], delimiter: str) -> Optional[int]:
        """
        Find the header row containing column names.
        
        Args:
            lines: List of cleaned data lines
            delimiter: Delimiter character
            
        Returns:
            Index of header row, or None if not found
        """
        for i, line in enumerate(lines):
            if not line.strip():
                continue
                
            fields = line.split(delimiter)
            
            # Header rows typically have:
            # 1. Multiple fields
            # 2. Mostly non-numeric content
            # 3. Common gene column names
            if len(fields) > 1:
                non_numeric_count = sum(1 for field in fields 
                                      if not self._is_numeric(field.strip()))
                
                # Check for gene column indicators
                has_gene_column = any(
                    re.match(r'^(gene|symbol|id_ref|probe)', field.lower().strip())
                    for field in fields
                )
                
                if (non_numeric_count > len(fields) * 0.5 or has_gene_column):
                    self.logger.info(f"Header row detected at line {i + 1}")
                    return i
        
        self.logger.warning("No header row detected")
        return None
    
    def identify_gene_column(self, header: List[str]) -> int:
        """
        Identify which column contains gene names/IDs.
        
        Args:
            header: List of column headers
            
        Returns:
            Index of gene column
        """
        gene_patterns = [
            r'^gene_?symbol',
            r'^gene_?name',
            r'^gene_?id',
            r'^id_ref',
            r'^probe_?id',
            r'^probeset_?id',
            r'^symbol',
            r'^ensembl',
            r'^entrez',
            r'^refseq',
            r'^transcript',
            r'^feature_?id',
            r'^name$',
            r'^id$'
        ]
        
        # Check each column header against patterns
        for i, col_name in enumerate(header):
            col_clean = col_name.lower().strip().replace('"', '')
            
            for pattern in gene_patterns:
                if re.match(pattern, col_clean):
                    self.logger.info(f"Gene column: '{col_name}' at index {i}")
                    return i
        
        # Default to first column if no pattern matches
        self.logger.info(f"Using first column as gene column: '{header[0]}'")
        return 0
    
    def identify_sample_columns(self, header: List[str], gene_col_idx: int) -> List[int]:
        """
        Identify columns containing sample expression data.
        
        Args:
            header: List of column headers
            gene_col_idx: Index of gene column to exclude
            
        Returns:
            List of sample column indices
        """
        sample_columns = []
        
        # Patterns for annotation columns to skip
        annotation_patterns = [
            r'^description',
            r'^annotation',
            r'^definition',
            r'^entrez_?gene',
            r'^go_',
            r'^pathway',
            r'^chromosome',
            r'^chr',
            r'^location',
            r'^position',
            r'^strand',
            r'^length',
            r'^gc_?content',
            r'^synonyms?',
            r'^unigene',
            r'^refseq',
            r'^ensembl',
            r'^omim',
            r'^type',
            r'^class',
            r'^category',
            r'^accession',
            r'^platform',
            r'^sequence'
        ]
        
        for i, col_name in enumerate(header):
            if i == gene_col_idx:
                continue
                
            col_clean = col_name.lower().strip().replace('"', '')
            
            # Check if this is an annotation column
            is_annotation = any(re.match(pattern, col_clean) 
                              for pattern in annotation_patterns)
            
            if not is_annotation:
                sample_columns.append(i)
        
        self.logger.info(f"Found {len(sample_columns)} sample columns")
        return sample_columns
    
    # ========================
    # DATA VALIDATION METHODS
    # ========================
    
    def _is_metadata_line(self, line: str) -> bool:
        """Check if a line contains metadata that should be removed."""
        line = line.strip()
        if not line:
            return True
            
        # GEO and common metadata patterns
        metadata_patterns = [
            r'^!',              # GEO metadata
            r'^#',              # Comments
            r'^@',              # Attributes
            r'^\^',             # Series/sample info
            r'^"!',             # Quoted metadata
            r'^Series_',        # Series information
            r'^Sample_',        # Sample information
            r'^Platform_',      # Platform information
            r'^Database_',      # Database information
            r'^Annotation_',    # Annotation information
            r'^\s*$',           # Empty lines
            r'^"?ID_REF"?\s+"?VALUE"?',  # Some GEO headers
            r'^\d+/\d+/\d+',    # Date lines
            r'^Contact:',       # Contact info
            r'^Contributor:',   # Contributor info
            r'^Data processing:', # Processing info
        ]
        
        return any(re.match(pattern, line, re.IGNORECASE) 
                  for pattern in metadata_patterns)
    
    def _is_footnote_line(self, line: str) -> bool:
        """Check if a line is a footnote or separator."""
        line = line.strip()
        if not line:
            return True
            
        footnote_patterns = [
            r'^\*',             # Asterisk footnotes
            r'^Note:',          # Note lines
            r'^Footer:',        # Footer lines
            r'^Footnote:',      # Footnote lines
            r'^\d+\.',          # Numbered footnotes
            r'^[a-zA-Z]\.',     # Lettered footnotes
            r'^-{3,}',          # Dash separators
            r'^={3,}',          # Equals separators
            r'^_{3,}',          # Underscore separators
            r'^\s*---\s*$',     # Line separators
        ]
        
        return any(re.match(pattern, line, re.IGNORECASE) 
                  for pattern in footnote_patterns)
    
    def _is_numeric(self, value: str) -> bool:
        """Check if a value can be converted to a number."""
        if not value or value.isspace():
            return False
            
        try:
            float(value)
            return True
        except ValueError:
            return False
    
    def _is_consistent_row(self, row: List[str], expected_cols: int) -> bool:
        """
        Check if a data row has consistent structure.
        
        Args:
            row: List of field values
            expected_cols: Expected number of columns
            
        Returns:
            True if row is consistent with expected format
        """
        # Must have at least 70% of expected columns
        if len(row) < expected_cols * 0.7:
            return False
        
        # Must have reasonable amount of non-empty data
        non_empty_count = sum(1 for field in row if field.strip())
        if non_empty_count < expected_cols * 0.3:
            return False
        
        return True
    
    # ========================
    # DATA PROCESSING METHODS
    # ========================
    
    def _clean_field(self, field: str) -> str:
        """
        Clean individual field values.
        
        Args:
            field: Raw field value
            
        Returns:
            Cleaned field value
        """
        # Remove quotes and whitespace
        field = field.strip().strip('"\'')
        
        # Handle missing/null values
        if field.lower() in ['null', 'na', 'n/a', 'nan', '', 'none', 'nd']:
            return 'NA'
        
        # Remove excessive whitespace
        field = ' '.join(field.split())
        
        return field
    
    def _format_gene_expression_matrix(self, header: List[str], 
                                     rows: List[List[str]]) -> Tuple[List[str], List[List[str]]]:
        """
        Format data as gene expression matrix.
        
        Args:
            header: Column headers
            rows: Data rows
            
        Returns:
            Tuple of (formatted_header, formatted_rows)
        """
        # Identify gene and sample columns
        gene_col_idx = self.identify_gene_column(header)
        sample_col_indices = self.identify_sample_columns(header, gene_col_idx)
        
        if not sample_col_indices:
            self.logger.warning("No sample columns found, using all non-gene columns")
            sample_col_indices = [i for i in range(len(header)) if i != gene_col_idx]
        
        # Create new header
        sample_names = [header[i] for i in sample_col_indices]
        new_header = ['Gene'] + sample_names
        
        # Process data rows
        formatted_rows = []
        seen_genes = set()
        
        for row in rows:
            if len(row) <= gene_col_idx:
                continue
                
            # Get gene name
            gene_name = self._clean_field(row[gene_col_idx])
            if not gene_name or gene_name == 'NA':
                continue
            
            # Handle duplicate gene names
            original_gene = gene_name
            counter = 1
            while gene_name in seen_genes:
                gene_name = f"{original_gene}_{counter}"
                counter += 1
            
            if counter > 1:
                self.logger.warning(f"Duplicate gene: {original_gene} -> {gene_name}")
            
            seen_genes.add(gene_name)
            
            # Extract sample values
            sample_values = []
            for col_idx in sample_col_indices:
                if col_idx < len(row):
                    value = self._clean_field(row[col_idx])
                    sample_values.append(value)
                else:
                    sample_values.append('NA')
            
            formatted_rows.append([gene_name] + sample_values)
        
        return new_header, formatted_rows
    
    # ========================
    # MAIN PROCESSING METHODS
    # ========================
    
    def process_file(self, input_file: str, output_file: str = None) -> Optional[str]:
        """
        Process a single .txt.gz file.
        
        Args:
            input_file: Path to input .txt.gz file
            output_file: Path to output .tsv file (optional)
            
        Returns:
            Path to output file if successful, None otherwise
        """
        # Set default output filename
        if output_file is None:
            output_file = str(Path(input_file).with_suffix('.tsv'))
        
        self.logger.info(f"Processing: {input_file} -> {output_file}")
        self.reset_stats()
        
        try:
            # Read compressed file
            with gzip.open(input_file, 'rt', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
        except Exception as e:
            self.logger.error(f"Error reading {input_file}: {e}")
            return None
        
        self.stats['total_lines'] = len(lines)
        self.logger.info(f"Read {len(lines)} lines")
        
        # Clean lines (remove metadata and footnotes)
        clean_lines = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if self._is_metadata_line(line):
                self.stats['metadata_lines'] += 1
                continue
                
            if self._is_footnote_line(line):
                self.stats['footnote_lines'] += 1
                continue
            
            clean_lines.append(line)
        
        if not clean_lines:
            self.logger.error("No valid data lines found")
            return None
        
        # Detect delimiter and header
        delimiter = self.detect_delimiter(clean_lines)
        header_idx = self.detect_header_row(clean_lines, delimiter)
        
        if header_idx is not None:
            header_line = clean_lines[header_idx]
            data_lines = clean_lines[header_idx + 1:]
        else:
            # Generate generic header
            first_line = clean_lines[0]
            num_cols = len(first_line.split(delimiter))
            header_line = delimiter.join([f'Column_{i+1}' for i in range(num_cols)])
            data_lines = clean_lines
            self.logger.info("Generated generic header")
        
        # Parse header
        header = [self._clean_field(field) for field in header_line.split(delimiter)]
        expected_columns = len(header)
        
        # Process data rows
        processed_rows = []
        for line in data_lines:
            row = [self._clean_field(field) for field in line.split(delimiter)]
            
            # Check row consistency
            if not self._is_consistent_row(row, expected_columns):
                self.stats['inconsistent_lines'] += 1
                continue
            
            # Pad or truncate to match header length
            while len(row) < expected_columns:
                row.append('NA')
            row = row[:expected_columns]
            
            processed_rows.append(row)
        
        # Format as gene expression matrix
        formatted_header, formatted_rows = self._format_gene_expression_matrix(
            header, processed_rows
        )
        
        # Update statistics
        self.stats['final_genes'] = len(formatted_rows)
        self.stats['final_samples'] = len(formatted_header) - 1  # Exclude gene column
        
        # Write output file
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                # Write header
                f.write('\t'.join(formatted_header) + '\n')
                
                # Write data rows
                for row in formatted_rows:
                    f.write('\t'.join(str(cell) for cell in row) + '\n')
            
            self.logger.info(f"Successfully created: {output_file}")
            self._print_stats()
            
            return output_file
            
        except Exception as e:
            self.logger.error(f"Error writing {output_file}: {e}")
            return None
    
    def batch_process(self, input_dir: str, output_dir: str = None) -> List[str]:
        """
        Process multiple .txt.gz files in a directory.
        
        Args:
            input_dir: Directory containing .txt.gz files
            output_dir: Output directory (defaults to input_dir)
            
        Returns:
            List of successfully processed output files
        """
        input_path = Path(input_dir)
        
        if output_dir is None:
            output_dir = input_dir
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Find all .txt.gz files
        txtgz_files = list(input_path.glob('*.txt.gz'))
        
        if not txtgz_files:
            self.logger.warning(f"No .txt.gz files found in {input_dir}")
            return []
        
        self.logger.info(f"Processing {len(txtgz_files)} files")
        
        # Process each file
        successful_files = []
        for txtgz_file in txtgz_files:
            output_file = output_path / f"{txtgz_file.stem}.tsv"
            result = self.process_file(str(txtgz_file), str(output_file))
            
            if result:
                successful_files.append(result)
        
        self.logger.info(f"Successfully processed {len(successful_files)}/{len(txtgz_files)} files")
        return successful_files
    
    def validate_output(self, tsv_file: str) -> bool:
        """
        Validate the generated TSV file.
        
        Args:
            tsv_file: Path to TSV file to validate
            
        Returns:
            True if validation passes, False otherwise
        """
        try:
            df = pd.read_csv(tsv_file, sep='\t', index_col=0)
            
            self.logger.info(f"Validation successful:")
            self.logger.info(f"  - Genes: {len(df)}")
            self.logger.info(f"  - Samples: {len(df.columns)}")
            self.logger.info(f"  - Sample names: {list(df.columns[:5])}...")
            self.logger.info(f"  - Gene examples: {list(df.index[:5])}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Validation failed: {e}")
            return False
    
    def _print_stats(self) -> None:
        """Print processing statistics."""
        print("\n" + "="*60)
        print("GENE EXPRESSION CONVERSION STATISTICS")
        print("="*60)
        print(f"Total input lines:        {self.stats['total_lines']:,}")
        print(f"Metadata lines removed:   {self.stats['metadata_lines']:,}")
        print(f"Footnote lines removed:   {self.stats['footnote_lines']:,}")
        print(f"Inconsistent rows removed: {self.stats['inconsistent_lines']:,}")
        print(f"Final genes:              {self.stats['final_genes']:,}")
        print(f"Final samples:            {self.stats['final_samples']:,}")
        print("="*60)


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Convert .txt.gz files to gene expression matrices (TSV format)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s GSE12345.txt.gz                    # Convert single file
  %(prog)s GSE12345.txt.gz -o matrix.tsv      # Specify output name
  %(prog)s /path/to/data/ -b                  # Batch process directory
  %(prog)s data.txt.gz -v --validate          # Verbose with validation
        """
    )
    
    parser.add_argument('input', 
                       help='Input .txt.gz file or directory')
    parser.add_argument('-o', '--output', 
                       help='Output TSV file or directory')
    parser.add_argument('-b', '--batch', 
                       action='store_true',
                       help='Batch process directory')
    parser.add_argument('-v', '--verbose', 
                       action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--validate', 
                       action='store_true',
                       help='Validate output file(s)')
    
    args = parser.parse_args()
    
    # Configure logging level
    log_level = logging.DEBUG if args.verbose else logging.INFO
    
    # Create converter
    converter = GeneExpressionConverter(log_level=log_level)
    
    try:
        if args.batch:
            # Batch processing
            output_files = converter.batch_process(args.input, args.output)
            
            if args.validate:
                print(f"\nValidating {len(output_files)} output files...")
                for output_file in output_files:
                    print(f"Validating: {output_file}")
                    converter.validate_output(output_file)
        else:
            # Single file processing
            output_file = converter.process_file(args.input, args.output)
            
            if output_file and args.validate:
                print(f"\nValidating output file...")
                converter.validate_output(output_file)
    
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

"""
Python Script Usage:
--------------------

# Basic conversion
converter = GeneExpressionConverter()
output_file = converter.process_file('GSE12345.txt.gz')

# Batch processing with validation
converter = GeneExpressionConverter(log_level=logging.DEBUG)
output_files = converter.batch_process('/path/to/geo_data/', '/path/to/output/')

for output_file in output_files:
    converter.validate_output(output_file)

# Custom output filename
converter.process_file('expression_data.txt.gz', 'my_matrix.tsv')

Command Line Usage:
-------------------

# Convert single file
python gene_converter.py GSE12345.txt.gz

# Convert with custom output name
python gene_converter.py GSE12345.txt.gz -o expression_matrix.tsv

# Batch process directory
python gene_converter.py /path/to/data/ -b -o /path/to/output/

# Verbose mode with validation
python gene_converter.py data.txt.gz -v --validate

Output Format:
--------------
Gene      Sample1    Sample2    Sample3    Sample4
ACTB      12.5       11.8       13.2       12.1
GAPDH     15.2       14.9       15.8       15.1  
TP53      8.3        9.1        7.9        8.7
BRCA1     6.7        6.2        6.9        6.4
"""