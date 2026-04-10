#!/usr/bin/env python3
"""
Feature Generation Module for ARMOR Pipeline

This module provides utilities for generating genomic features from WGS data:
- Gene presence/absence matrices
- SNP matrices
- K-mer frequency matrices
"""

import os
import subprocess
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path

import pandas as pd
import numpy as np
from tqdm import tqdm

logger = logging.getLogger(__name__)


class FeatureGenerator:
    """Generate genomic features from WGS data."""

    def __init__(
        self,
        reference_fasta: Optional[str] = None,
        reference_gff: Optional[str] = None,
        threads: int = 8,
        kmer_size: int = 11,
    ):
        """
        Initialize feature generator.

        Args:
            reference_fasta: Path to reference genome FASTA
            reference_gff: Path to reference genome annotation
            threads: Number of threads for parallel processing
            kmer_size: K-mer size for k-mer counting
        """
        self.reference_fasta = reference_fasta
        self.reference_gff = reference_gff
        self.threads = threads
        self.kmer_size = kmer_size

    def check_dependencies(self) -> Dict[str, bool]:
        """Check if required external tools are installed."""
        tools = {
            "bakta": self._check_command("bakta"),
            "panaroo": self._check_command("panaroo"),
            "snippy": self._check_command("snippy"),
            "gubbins": self._check_command("run_gubbins.py"),
            "jellyfish": self._check_command("jellyfish"),
            "kmc": self._check_command("kmc"),
        }
        return tools

    def _check_command(self, cmd: str) -> bool:
        """Check if a command is available in PATH."""
        try:
            result = subprocess.run(
                ["which", cmd],
                capture_output=True,
                text=True,
            )
            return result.returncode == 0
        except Exception:
            return False

    def generate_gene_matrix(
        self,
        genome_dir: str,
        output_dir: str,
        skip_annotation: bool = False,
    ) -> str:
        """
        Generate gene presence/absence matrix using Bakta and Panaroo.

        Args:
            genome_dir: Directory containing genome FASTA files
            output_dir: Output directory for gene matrix
            skip_annotation: Skip annotation if GFF files exist

        Returns:
            Path to generated gene matrix file
        """
        logger.info("Generating gene presence/absence matrix...")

        # Step 1: Annotation with Bakta
        if not skip_annotation:
            gff_files = self._run_bakta(genome_dir, output_dir)
        else:
            gff_files = list(Path(output_dir).glob("*.gff*"))

        # Step 2: Pangenome analysis with Panaroo
        matrix_file = self._run_panaroo(gff_files, output_dir)

        return matrix_file

    def _run_bakta(self, genome_dir: str, output_dir: str) -> List[str]:
        """Run Bakta annotation on genome files."""
        logger.info("Running Bakta annotation...")

        gff_files = []
        fasta_files = list(Path(genome_dir).glob("*.fna"))
        fasta_files.extend(Path(genome_dir).glob("*.fasta"))

        for fasta in tqdm(fasta_files, desc="Annotating genomes"):
            sample_id = fasta.stem
            sample_out = Path(output_dir) / sample_id

            cmd = [
                "bakta",
                "--db", "db",
                "--output", str(sample_out),
                "--threads", str(self.threads),
                "--prefix", sample_id,
                str(fasta),
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)
                gff_files.append(sample_out / f"{sample_id}.gff3")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to annotate {sample_id}: {e}")

        return gff_files

    def _run_panaroo(self, gff_files: List[Path], output_dir: str) -> str:
        """Run Panaroo pangenome analysis."""
        logger.info("Running Panaroo pangenome analysis...")

        gff_list = ",".join(str(f) for f in gff_files)
        out_prefix = Path(output_dir) / "panaroo"

        cmd = [
            "panaroo",
            "-i", gff_list,
            "-o", str(out_prefix),
            "--clean-mode", "strict",
            "-t", str(self.threads),
            "--algos", "gene",
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Panaroo failed: {e}")
            raise

        # Convert Panaroo output to binary matrix
        matrix_file = self._convert_panaroo_output(
            Path(output_dir) / "panaroo_gene_presence_absence.CSV",
            Path(output_dir) / "gene_matrix.txt",
        )

        return matrix_file

    def _convert_panaroo_output(self, panaroo_file: Path, out_file: Path) -> str:
        """Convert Panaroo CSV to binary presence/absence matrix."""
        logger.info("Converting Panaroo output...")

        # Read Panaroo output
        pa_data = pd.read_csv(panaroo_file, sep="\t", header=0)

        # Extract gene names (first column)
        gene_names = pa_data.iloc[:, 0]

        # Sample columns (columns 9 onwards)
        sample_cols = pa_data.columns[8:]

        # Create binary matrix
        binary_matrix = pd.DataFrame(
            index=gene_names,
            columns=sample_cols,
        )

        for col in sample_cols:
            # Any non-empty value indicates presence
            binary_matrix[col] = pa_data[col].apply(lambda x: 0 if pd.isna(x) or x == "" or x == "-" else 1)

        # Transpose to have samples as rows
        binary_matrix = binary_matrix.T
        binary_matrix.index.name = "ID"
        binary_matrix = binary_matrix.reset_index()

        # Save
        binary_matrix.to_csv(out_file, sep="\t", index=False)
        logger.info(f"Gene matrix saved to {out_file}")

        return str(out_file)

    def generate_snp_matrix(
        self,
        fastq_dir: str,
        output_dir: str,
        skip_gubbins: bool = False,
    ) -> str:
        """
        Generate SNP matrix using Snippy and Gubbins.

        Args:
            fastq_dir: Directory containing FASTQ files
            output_dir: Output directory for SNP matrix
            skip_gubbins: Skip Gubbins recombination filtering

        Returns:
            Path to generated SNP matrix file
        """
        logger.info("Generating SNP matrix...")

        # Step 1: Variant calling with Snippy
        snippy_dirs = self._run_snippy(fastq_dir, output_dir)

        # Step 2: Core SNP identification
        core_vcf, core_aln = self._run_snippy_core(snippy_dirs, output_dir)

        # Step 3: Recombination filtering with Gubbins
        if not skip_gubbins:
            filtered_aln = self._run_gubbins(core_aln, output_dir)
        else:
            filtered_aln = core_aln

        # Step 4: Convert to binary matrix
        matrix_file = self._convert_alignment_to_matrix(
            filtered_aln,
            Path(output_dir) / "SNP_matrix.txt",
        )

        return matrix_file

    def _run_snippy(self, fastq_dir: str, output_dir: str) -> List[Path]:
        """Run Snippy variant calling."""
        logger.info("Running Snippy variant calling...")

        snippy_dirs = []
        fastq_files = list(Path(fastq_dir).glob("*.fastq"))
        fastq_files.extend(Path(fastq_dir).glob("*.fq"))

        for fastq in tqdm(fastq_files, desc="Variant calling"):
            if "_R2" in fastq.name or "_r2" in fastq.name:
                continue  # Skip R2 files

            sample_id = fastq.stem.replace("_R1", "").replace("_r1", "")
            sample_out = Path(output_dir) / "snippy" / sample_id

            # Check if paired
            r2_file = fastq.parent / fastq.name.replace("_R1", "_R2").replace("_r1", "_r2")
            if r2_file.exists():
                cmd = [
                    "snippy",
                    "--cpus", str(self.threads),
                    "--ref", self.reference_fasta,
                    "--R1", str(fastq),
                    "--R2", str(r2_file),
                    "--outdir", str(sample_out),
                    "--prefix", sample_id,
                ]
            else:
                cmd = [
                    "snippy",
                    "--cpus", str(self.threads),
                    "--ref", self.reference_fasta,
                    "--se", str(fastq),
                    "--outdir", str(sample_out),
                    "--prefix", sample_id,
                ]

            try:
                subprocess.run(cmd, check=True)
                snippy_dirs.append(sample_out)
            except subprocess.CalledProcessError as e:
                logger.error(f"Snippy failed for {sample_id}: {e}")

        return snippy_dirs

    def _run_snippy_core(self, snippy_dirs: List[Path], output_dir: str) -> Tuple[str, str]:
        """Run Snippy-core to identify core SNPs."""
        logger.info("Running Snippy-core...")

        core_out = Path(output_dir) / "core"
        core_out.mkdir(exist_ok=True)

        ref_fasta = snippy_dirs[0] / f"{snippy_dirs[0].name}.fasta"

        dirs_list = ",".join(str(d) for d in snippy_dirs)
        cmd = [
            "snippy-core",
            "--ref", str(ref_fasta),
            "--dir", str(core_out),
            dirs_list,
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Snippy-core failed: {e}")
            raise

        return str(core_out / "core.vcf"), str(core_out / "core.full.aln")

    def _run_gubbins(self, alignment_file: str, output_dir: str) -> str:
        """Run Gubbins to filter recombinant regions."""
        logger.info("Running Gubbins...")

        gubbins_out = Path(output_dir) / "gubbins"
        gubbins_out.mkdir(exist_ok=True)

        cmd = [
            "run_gubbins.py",
            "--threads", str(self.threads),
            "--output_dir", str(gubbins_out),
            alignment_file,
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Gubbins failed: {e}")
            raise

        return str(gubbins_out / "gubbins.final_tree.filtered_polymorphic_sites.fasta")

    def _convert_alignment_to_matrix(self, alignment_file: str, out_file: Path) -> str:
        """Convert alignment to binary SNP matrix."""
        logger.info("Converting alignment to SNP matrix...")

        # Read alignment (FASTA format)
        from Bio import SeqIO

        sequences = {}
        for record in SeqIO.parse(alignment_file, "fasta"):
            sequences[record.id] = str(record.seq)

        # Find variant positions
        seq_ids = list(sequences.keys())
        seq_length = len(sequences[seq_ids[0]])

        # Identify polymorphic sites
        variant_positions = []
        for pos in range(seq_length):
            bases = set(seq[pos] for seq in sequences.values())
            if len(bases) > 1 and "-" not in bases:
                variant_positions.append(pos)

        logger.info(f"Found {len(variant_positions)} polymorphic sites")

        # Create binary matrix
        matrix_data = {"ID": seq_ids}
        for pos in variant_positions:
            ref_base = sequences[seq_ids[0]][pos]
            col_name = f"SNP_{pos+1}_{ref_base}"
            matrix_data[col_name] = [
                0 if sequences[sid][pos] == ref_base else 1
                for sid in seq_ids
            ]

        matrix_df = pd.DataFrame(matrix_data)
        matrix_df.to_csv(out_file, sep="\t", index=False)
        logger.info(f"SNP matrix saved to {out_file}")

        return str(out_file)

    def generate_kmer_matrix(
        self,
        fastq_dir: str,
        output_dir: str,
        min_presence: float = 0.01,
        max_presence: float = 0.99,
    ) -> str:
        """
        Generate k-mer frequency matrix using Jellyfish or KMC.

        Args:
            fastq_dir: Directory containing FASTQ files
            output_dir: Output directory for k-mer matrix
            min_presence: Minimum presence rate across samples
            max_presence: Maximum presence rate across samples

        Returns:
            Path to generated k-mer matrix file
        """
        logger.info("Generating k-mer matrix...")

        # Check which tool is available
        use_kmc = self._check_command("kmc")

        if use_kmc:
            logger.info("Using KMC for k-mer counting")
            kmer_counts = self._run_kmc(fastq_dir, output_dir)
        else:
            logger.info("Using Jellyfish for k-mer counting")
            kmer_counts = self._run_jellyfish(fastq_dir, output_dir)

        # Filter and merge k-mers
        matrix_file = self._merge_kmer_counts(
            kmer_counts,
            output_dir,
            min_presence,
            max_presence,
        )

        return matrix_file

    def _run_jellyfish(self, fastq_dir: str, output_dir: str) -> Dict[str, Path]:
        """Run Jellyfish k-mer counting."""
        temp_dir = Path(output_dir) / "jf_temp"
        temp_dir.mkdir(exist_ok=True)

        kmer_counts = {}
        fastq_files = list(Path(fastq_dir).glob("*.fastq"))

        for fastq in tqdm(fastq_files, desc="Counting k-mers"):
            if "_R2" in fastq.name:
                continue

            sample_id = fastq.stem.replace("_R1", "")
            out_prefix = temp_dir / sample_id

            cmd = [
                "jellyfish",
                "count",
                "-C",  # Canonical k-mers
                "-m", str(self.kmer_size),
                "-s", "100M",
                "-t", str(self.threads),
                "-o", str(out_prefix) + ".jf",
                str(fastq),
            ]

            try:
                subprocess.run(cmd, check=True)

                # Dump to text
                dump_file = out_prefix.with_suffix(".txt")
                subprocess.run(
                    ["jellyfish", "dump", "-o", str(dump_file), str(out_prefix) + ".jf"],
                    check=True,
                )
                kmer_counts[sample_id] = dump_file

            except subprocess.CalledProcessError as e:
                logger.error(f"Jellyfish failed for {sample_id}: {e}")

        return kmer_counts

    def _run_kmc(self, fastq_dir: str, output_dir: str) -> Dict[str, Path]:
        """Run KMC k-mer counting."""
        temp_dir = Path(output_dir) / "kmc_temp"
        temp_dir.mkdir(exist_ok=True)

        kmer_counts = {}
        fastq_files = list(Path(fastq_dir).glob("*.fastq"))

        for fastq in tqdm(fastq_files, desc="Counting k-mers (KMC)"):
            if "_R2" in fastq.name:
                continue

            sample_id = fastq.stem.replace("_R1", "")
            out_prefix = temp_dir / sample_id

            cmd = [
                "kmc",
                f"-k{self.kmer_size}",
                f"-t{self.threads}",
                "-m64",
                "-ci2",
                "-cs10000",
                str(fastq),
                str(out_prefix),
                str(temp_dir),
            ]

            try:
                subprocess.run(cmd, check=True)

                # Dump results
                dump_file = out_prefix.with_suffix(".txt")
                subprocess.run(
                    ["kmc_tools", "transform", str(out_prefix), "dump", str(dump_file)],
                    check=True,
                )
                kmer_counts[sample_id] = dump_file

            except subprocess.CalledProcessError as e:
                logger.error(f"KMC failed for {sample_id}: {e}")

        return kmer_counts

    def _merge_kmer_counts(
        self,
        kmer_counts: Dict[str, Path],
        output_dir: str,
        min_presence: float,
        max_presence: float,
    ) -> str:
        """Merge k-mer counts from multiple samples."""
        logger.info("Merging k-mer counts...")

        # Read all k-mer counts
        all_kmers = {}
        for sample_id, count_file in kmer_counts.items():
            kmers = {}
            with open(count_file) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) == 2:
                        kmer, count = parts
                        kmers[kmer] = int(count)
            all_kmers[sample_id] = kmers

        # Get all unique k-mers
        unique_kmers = set()
        for kmers in all_kmers.values():
            unique_kmers.update(kmers.keys())

        logger.info(f"Total unique k-mers: {len(unique_kmers)}")

        # Create matrix
        n_samples = len(all_kmers)
        matrix = pd.DataFrame(
            index=list(all_kmers.keys()),
            columns=list(unique_kmers),
            dtype=float,
        )

        for sample_id, kmers in all_kmers.items():
            for kmer in unique_kmers:
                matrix.loc[sample_id, kmer] = kmers.get(kmer, 0)

        # Filter by presence rate
        presence_rate = (matrix > 0).sum(axis=0) / n_samples
        keep_kmers = (presence_rate >= min_presence) & (presence_rate <= max_presence)

        logger.info(f"K-mers after filtering: {keep_kmers.sum()}")

        matrix = matrix.loc[:, keep_kmers]
        matrix = matrix.reset_index().rename(columns={"index": "ID"})

        # Save
        out_file = Path(output_dir) / "kmer_matrix.txt"
        matrix.to_csv(out_file, sep="\t", index=False)
        logger.info(f"K-mer matrix saved to {out_file}")

        return str(out_file)


def main():
    """Main entry point for feature generation."""
    import click

    @click.group()
    def cli():
        """ARMOR Feature Generation Tool"""
        pass

    @cli.command()
    @click.option("--input", "-i", required=True, help="Genome directory")
    @click.option("--output", "-o", required=True, help="Output directory")
    @click.option("--threads", "-t", default=8, help="Number of threads")
    def gene(input, output, threads):
        """Generate gene presence/absence matrix."""
        generator = FeatureGenerator(threads=threads)
        generator.generate_gene_matrix(input, output)

    @cli.command()
    @click.option("--input", "-i", required=True, help="FASTQ directory")
    @click.option("--reference", "-r", required=True, help="Reference FASTA")
    @click.option("--output", "-o", required=True, help="Output directory")
    @click.option("--threads", "-t", default=8, help="Number of threads")
    @click.option("--skip-gubbins", is_flag=True, help="Skip Gubbins filtering")
    def snp(input, reference, output, threads, skip_gubbins):
        """Generate SNP matrix."""
        generator = FeatureGenerator(reference_fasta=reference, threads=threads)
        generator.generate_snp_matrix(input, output, skip_gubbins)

    @cli.command()
    @click.option("--input", "-i", required=True, help="FASTQ directory")
    @click.option("--output", "-o", required=True, help="Output directory")
    @click.option("--kmer-size", "-k", default=11, help="K-mer size")
    @click.option("--threads", "-t", default=8, help="Number of threads")
    def kmer(input, output, kmer_size, threads):
        """Generate k-mer frequency matrix."""
        generator = FeatureGenerator(threads=threads, kmer_size=kmer_size)
        generator.generate_kmer_matrix(input, output)

    cli()


if __name__ == "__main__":
    main()
