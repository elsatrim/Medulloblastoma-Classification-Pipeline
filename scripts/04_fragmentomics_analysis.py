#!/usr/bin/env python3
"""
Fragmentomics Analysis Module for MB Classification Pipeline
Extracts fragment patterns from lcWGS BAM files for medulloblastoma subgroup classification

Based on: Markowitz et al. 2025 - npj Precision Oncology
Author: Elsa Sanchez Fernandez
Date: 2024
"""

import pysam
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FragmentomicsAnalyzer:
    """
    Extract and analyze cfDNA fragmentation patterns from BAM files.
    
    Features extracted:
    1. Fragment length distribution
    2. Short-to-long (S/L) ratio
    3. 4bp fragment end motifs
    4. Genome-wide S/L ratios (5 Mb bins)
    """
    
    def __init__(self, bam_file, reference_fasta=None, output_dir=None):
        """
        Initialize fragmentomics analyzer.
        
        Args:
            bam_file: Path to aligned BAM file
            reference_fasta: Path to reference genome (optional, for motif extraction)
            output_dir: Directory for output files
        """
        self.bam_file = Path(bam_file)
        self.reference_fasta = reference_fasta
        self.output_dir = Path(output_dir) if output_dir else self.bam_file.parent
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.sample_name = self.bam_file.stem.replace('.bam', '')
        
        # Storage for extracted features
        self.fragment_lengths = []
        self.end_motifs = []
        self.binned_sl_ratios = {}
        self.metrics = {}
        
        logger.info(f"Initialized fragmentomics analyzer for {self.sample_name}")
    
    
    def extract_fragments(self, min_length=50, max_length=600):
        """
        Extract fragment information from BAM file.
        
        Args:
            min_length: Minimum fragment length to consider (default: 50)
            max_length: Maximum fragment length to consider (default: 600)
        """
        logger.info(f"Opening BAM file: {self.bam_file}")
        
        try:
            bamfile = pysam.AlignmentFile(str(self.bam_file), "rb")
        except Exception as e:
            logger.error(f"Failed to open BAM file: {e}")
            raise
        
        # Open reference if provided
        ref = None
        if self.reference_fasta:
            try:
                ref = pysam.FastaFile(str(self.reference_fasta))
                logger.info(f"Loaded reference genome: {self.reference_fasta}")
            except Exception as e:
                logger.warning(f"Could not load reference genome: {e}")
                logger.warning("End motif extraction will be skipped")
        
        logger.info("Extracting fragments...")
        fragment_count = 0
        processed_reads = 0
        
        for read in bamfile:
            processed_reads += 1
            
            # Progress indicator
            if processed_reads % 1000000 == 0:
                logger.info(f"Processed {processed_reads:,} reads, extracted {fragment_count:,} fragments")
            
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            
            # For paired-end: only process read1 to avoid duplicates
            # For single-end: process all reads
            if read.is_paired and not read.is_read1:
                continue
            
            # Extract fragment length
            if read.is_paired and read.template_length != 0:
                # Paired-end: use insert size
                frag_len = abs(read.template_length)
            else:
                # Single-end or unpaired: use read length
                frag_len = read.query_length if read.query_length else 0
            
            if frag_len > 0 and min_length < frag_len < max_length:
                self.fragment_lengths.append(frag_len)
                fragment_count += 1
            
            # Extract end motif if reference available
            if ref and read.query_sequence and len(read.query_sequence) >= 4:
                try:
                    # Get 5' end motif (first 4bp of fragment)
                    chrom = read.reference_name
                    pos = read.reference_start
                    
                    # Only extract from main chromosomes
                    if chrom in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
                        motif = ref.fetch(chrom, pos, pos + 4).upper()
                        if len(motif) == 4 and all(base in "ACGT" for base in motif):
                            self.end_motifs.append(motif)
                except:
                    pass  # Skip problematic regions
        
        bamfile.close()
        if ref:
            ref.close()
        
        logger.info(f"Extraction complete!")
        logger.info(f"  Total reads processed: {processed_reads:,}")
        logger.info(f"  Fragments extracted: {len(self.fragment_lengths):,}")
        if self.end_motifs:
            logger.info(f"  End motifs extracted: {len(self.end_motifs):,}")
    
    
    def calculate_metrics(self):
        """Calculate fragmentomics metrics."""
        logger.info("Calculating fragmentomics metrics...")
        
        if not self.fragment_lengths:
            logger.error("No fragments extracted! Run extract_fragments() first.")
            return None
        
        # Fragment length statistics
        lengths = np.array(self.fragment_lengths)
        
        self.metrics['n_fragments'] = len(lengths)
        self.metrics['mean_length'] = float(np.mean(lengths))
        self.metrics['median_length'] = float(np.median(lengths))
        self.metrics['std_length'] = float(np.std(lengths))
        self.metrics['min_length'] = float(np.min(lengths))
        self.metrics['max_length'] = float(np.max(lengths))
        
        # Short-to-long ratio (key metric from Markowitz 2025)
        short = np.sum((lengths >= 90) & (lengths <= 150))
        long = np.sum((lengths >= 151) & (lengths <= 220))
        self.metrics['short_count'] = int(short)
        self.metrics['long_count'] = int(long)
        self.metrics['short_long_ratio'] = float(short / long) if long > 0 else 0.0
        
        # Fragment length percentiles
        self.metrics['percentile_25'] = float(np.percentile(lengths, 25))
        self.metrics['percentile_50'] = float(np.percentile(lengths, 50))
        self.metrics['percentile_75'] = float(np.percentile(lengths, 75))
        
        # End motif frequencies
        if self.end_motifs:
            motif_counts = Counter(self.end_motifs)
            total_motifs = len(self.end_motifs)
            
            # Top 20 motifs (matching Markowitz 2025)
            top_motifs = dict(motif_counts.most_common(20))
            self.metrics['top_motifs'] = {
                motif: count / total_motifs 
                for motif, count in top_motifs.items()
            }
            self.metrics['n_unique_motifs'] = len(motif_counts)
            
            # Known Group 4 enriched motifs (from Markowitz 2025)
            group4_motifs = ['AAAA', 'TTTT', 'CCCA', 'TAAA']
            group4_freq = sum(motif_counts.get(m, 0) for m in group4_motifs) / total_motifs
            self.metrics['group4_motif_frequency'] = float(group4_freq)
            
            # Known Group 3 enriched motifs
            group3_motifs = ['GCTG', 'ACAC', 'ACTG']
            group3_freq = sum(motif_counts.get(m, 0) for m in group3_motifs) / total_motifs
            self.metrics['group3_motif_frequency'] = float(group3_freq)
        
        logger.info("Metrics calculation complete!")
        logger.info(f"  Mean fragment length: {self.metrics['mean_length']:.1f} bp")
        logger.info(f"  Short/Long ratio: {self.metrics['short_long_ratio']:.3f}")
        
        return self.metrics
    
    
    def plot_results(self, figsize=(15, 10)):
        """
        Create comprehensive fragmentomics plots.
        
        Args:
            figsize: Figure size (width, height)
        """
        logger.info("Creating fragmentomics plots...")
        
        if not self.fragment_lengths:
            logger.error("No fragments to plot! Run extract_fragments() first.")
            return
        
        fig = plt.figure(figsize=figsize)
        
        # 1. Fragment length distribution
        ax1 = plt.subplot(2, 2, 1)
        plt.hist(self.fragment_lengths, bins=100, edgecolor='black', alpha=0.7, color='steelblue')
        plt.axvline(150, color='red', linestyle='--', linewidth=2, label='Short/Long cutoff (150bp)')
        plt.axvline(np.mean(self.fragment_lengths), color='green', linestyle='--', 
                   linewidth=2, label=f'Mean ({np.mean(self.fragment_lengths):.1f}bp)')
        plt.xlabel('Fragment Length (bp)', fontsize=12)
        plt.ylabel('Count', fontsize=12)
        plt.title(f'Fragment Size Distribution - {self.sample_name}', fontsize=14, fontweight='bold')
        plt.legend()
        plt.grid(alpha=0.3)
        
        # 2. Cumulative distribution
        ax2 = plt.subplot(2, 2, 2)
        sorted_lengths = np.sort(self.fragment_lengths)
        cumulative = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
        plt.plot(sorted_lengths, cumulative, linewidth=2, color='steelblue')
        plt.axvline(150, color='red', linestyle='--', linewidth=2, label='Short/Long cutoff')
        plt.xlabel('Fragment Length (bp)', fontsize=12)
        plt.ylabel('Cumulative Proportion', fontsize=12)
        plt.title('Cumulative Fragment Length Distribution', fontsize=14, fontweight='bold')
        plt.legend()
        plt.grid(alpha=0.3)
        
        # 3. Short vs Long fragments comparison
        ax3 = plt.subplot(2, 2, 3)
        short = sum(1 for l in self.fragment_lengths if 90 <= l <= 150)
        long = sum(1 for l in self.fragment_lengths if 151 <= l <= 220)
        
        categories = ['Short\n(90-150bp)', 'Long\n(151-220bp)']
        counts = [short, long]
        colors = ['#FF6B6B', '#4ECDC4']
        bars = plt.bar(categories, counts, color=colors, edgecolor='black', linewidth=2)
        
        # Add count labels on bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count:,}\n({count/sum(counts)*100:.1f}%)',
                    ha='center', va='bottom', fontsize=11, fontweight='bold')
        
        plt.ylabel('Fragment Count', fontsize=12)
        plt.title(f'Short/Long Ratio = {short/long:.3f}', fontsize=14, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        
        # 4. Top 20 End Motifs (if available)
        ax4 = plt.subplot(2, 2, 4)
        if self.end_motifs and 'top_motifs' in self.metrics:
            motifs = list(self.metrics['top_motifs'].keys())[:20]
            freqs = [self.metrics['top_motifs'][m] for m in motifs]
            
            # Color Group 4 and Group 3 enriched motifs differently
            group4_motifs = ['AAAA', 'TTTT', 'CCCA', 'TAAA']
            group3_motifs = ['GCTG', 'ACAC', 'ACTG']
            
            colors = []
            for m in motifs:
                if m in group4_motifs:
                    colors.append('#4CAF50')  # Green for Group 4
                elif m in group3_motifs:
                    colors.append('#FF9800')  # Orange for Group 3
                else:
                    colors.append('#2196F3')  # Blue for others
            
            bars = plt.bar(range(len(motifs)), freqs, color=colors, edgecolor='black')
            plt.xticks(range(len(motifs)), motifs, rotation=90, fontsize=9)
            plt.ylabel('Frequency', fontsize=12)
            plt.title('Top 20 Fragment End Motifs (4bp)', fontsize=14, fontweight='bold')
            plt.grid(axis='y', alpha=0.3)
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='#4CAF50', edgecolor='black', label='Group 4 enriched'),
                Patch(facecolor='#FF9800', edgecolor='black', label='Group 3 enriched'),
                Patch(facecolor='#2196F3', edgecolor='black', label='Other')
            ]
            plt.legend(handles=legend_elements, loc='upper right', fontsize=9)
        else:
            plt.text(0.5, 0.5, 'End motifs not extracted\n(Reference genome required)', 
                    ha='center', va='center', fontsize=12, transform=ax4.transAxes)
            plt.axis('off')
        
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / f"{self.sample_name}_fragmentomics.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Plot saved: {output_file}")
        
        plt.close()
        
        return output_file
    
    
    def compare_to_literature(self):
        """
        Compare extracted metrics to published Group 4 patterns from Markowitz 2025.
        """
        logger.info("Comparing to literature (Markowitz et al. 2025)...")
        
        if not self.metrics:
            logger.error("No metrics calculated! Run calculate_metrics() first.")
            return None
        
        # Expected Group 4 patterns (from Markowitz 2025)
        # Note: Exact values would need to be extracted from the paper's supplementary data
        # These are approximate based on the paper's figures
        
        literature_values = {
            'group4': {
                'short_long_ratio_range': (0.8, 1.2),  # Approximate from paper
                'mean_length_range': (160, 180),  # Typical for Group 4
                'enriched_motifs': ['AAAA', 'TTTT', 'CCCA', 'TAAA'],
                'depleted_motifs': ['CCAC', 'CCTG']
            },
            'group3': {
                'short_long_ratio_range': (0.9, 1.3),
                'enriched_motifs': ['GCTG', 'ACAC', 'ACTG'],
                'depleted_motifs': ['GAAA']
            }
        }
        
        comparison = {
            'sample_name': self.sample_name,
            'short_long_ratio': self.metrics['short_long_ratio'],
            'mean_length': self.metrics['mean_length']
        }
        
        # Check Group 4 pattern match
        sl_ratio = self.metrics['short_long_ratio']
        g4_sl_range = literature_values['group4']['short_long_ratio_range']
        
        if g4_sl_range[0] <= sl_ratio <= g4_sl_range[1]:
            comparison['sl_ratio_match'] = 'Group 4 range'
            logger.info(f"  S/L ratio ({sl_ratio:.3f}) matches Group 4 range")
        else:
            comparison['sl_ratio_match'] = 'Outside Group 4 range'
            logger.info(f"  S/L ratio ({sl_ratio:.3f}) outside Group 4 range")
        
        # Check motif patterns
        if 'top_motifs' in self.metrics:
            top_10_motifs = list(self.metrics['top_motifs'].keys())[:10]
            
            g4_enriched = literature_values['group4']['enriched_motifs']
            g4_found = [m for m in g4_enriched if m in top_10_motifs]
            
            comparison['group4_enriched_motifs_found'] = g4_found
            comparison['group4_motif_match'] = len(g4_found) >= 2  # At least 2 out of 4
            
            if len(g4_found) >= 2:
                logger.info(f"  Found {len(g4_found)}/4 Group 4 enriched motifs: {g4_found}")
            else:
                logger.info(f"  Only found {len(g4_found)}/4 Group 4 enriched motifs")
        
        # Overall assessment
        if comparison.get('sl_ratio_match') == 'Group 4 range' and \
           comparison.get('group4_motif_match', False):
            comparison['overall_assessment'] = 'CONCORDANT with Group 4 signature'
            logger.info("\n✅ RESULT: Fragment patterns CONCORDANT with Group 4 signature")
        else:
            comparison['overall_assessment'] = 'Pattern match inconclusive'
            logger.info("\n⚠️  RESULT: Fragment pattern match inconclusive")
        
        return comparison
    
    
    def save_metrics(self):
        """Save metrics to JSON file."""
        if not self.metrics:
            logger.error("No metrics to save! Run calculate_metrics() first.")
            return None
        
        output_file = self.output_dir / f"{self.sample_name}_fragmentomics_metrics.json"
        
        with open(output_file, 'w') as f:
            json.dump(self.metrics, f, indent=2)
        
        logger.info(f"Metrics saved: {output_file}")
        return output_file
    
    
    def run_complete_analysis(self):
        """
        Run complete fragmentomics analysis pipeline.
        
        Returns:
            dict: Dictionary with all results and file paths
        """
        logger.info("=" * 70)
        logger.info("STARTING FRAGMENTOMICS ANALYSIS")
        logger.info("=" * 70)
        
        # Extract fragments
        self.extract_fragments()
        
        # Calculate metrics
        metrics = self.calculate_metrics()
        
        # Create plots
        plot_file = self.plot_results()
        
        # Compare to literature
        comparison = self.compare_to_literature()
        
        # Save metrics
        metrics_file = self.save_metrics()
        
        logger.info("=" * 70)
        logger.info("FRAGMENTOMICS ANALYSIS COMPLETE!")
        logger.info("=" * 70)
        
        results = {
            'metrics': metrics,
            'comparison': comparison,
            'plot_file': str(plot_file),
            'metrics_file': str(metrics_file)
        }
        
        return results


def main():
    """Command-line interface for fragmentomics analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Extract fragmentomics patterns from lcWGS BAM files'
    )
    parser.add_argument(
        'bam_file',
        help='Path to input BAM file'
    )
    parser.add_argument(
        '--reference',
        help='Path to reference genome FASTA (optional, for end motif extraction)'
    )
    parser.add_argument(
        '--output-dir',
        help='Output directory (default: same as BAM file)'
    )
    
    args = parser.parse_args()
    
    # Run analysis
    analyzer = FragmentomicsAnalyzer(
        bam_file=args.bam_file,
        reference_fasta=args.reference,
        output_dir=args.output_dir
    )
    
    results = analyzer.run_complete_analysis()
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Fragment count: {results['metrics']['n_fragments']:,}")
    print(f"Mean length: {results['metrics']['mean_length']:.1f} bp")
    print(f"Short/Long ratio: {results['metrics']['short_long_ratio']:.3f}")
    
    if 'comparison' in results and results['comparison']:
        print(f"\nLiterature comparison: {results['comparison']['overall_assessment']}")
    
    print(f"\nPlot: {results['plot_file']}")
    print(f"Metrics: {results['metrics_file']}")
    print("=" * 70)


if __name__ == '__main__':
    main()