"""
AlphaGenome Atlas Tests
=======================
Tests Atlas data transformation (tidy_atlas_scores), file readers (VCF/TSV/CSV),
Atlas caching logic, and CLI parser integration.

Usage:
    python test_atlas.py
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import anndata

from alphagenome.data import genome
from alphagenome_cli import (
    tidy_atlas_scores,
    read_variants_from_tsv_or_csv,
    read_variants_from_vcf,
    CachedAtlasClient,
    build_parser,
    parse_variant,
    parse_interval,
    get_variant_interval,
    resolve_tissue_term,
)


class TestAtlasIntegration(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.test_path = Path(self.test_dir.name)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_tidy_atlas_scores(self):
        """AnnData yapısının tidy pandas DataFrame'e dönüştürülmesini test et."""
        # Mock AnnData
        var_obj = genome.Variant(
            chromosome="chr22",
            position=36201698,
            reference_bases="A",
            alternate_bases="C",
        )
        obs_df = pd.DataFrame({
            "variant": [var_obj],
            "gene_name": ["APOL4"],
            "gene_id": ["ENSG00000100435.10"],
            "strand": ["+"],
        })
        var_df = pd.DataFrame({
            "name": ["colon_rna_seq"],
            "ontology_curie": ["UBERON:0001157"],
            "strand": ["+"],
        })
        X = np.array([[1.85]], dtype=np.float32)
        layers = {"quantiles": np.array([[0.992]], dtype=np.float32)}

        adata = anndata.AnnData(X=X, obs=obs_df, var=var_df, layers=layers)
        scores_map = {"GeneMaskLFCScorer": adata}

        df = tidy_atlas_scores(scores_map)
        self.assertFalse(df.empty)
        self.assertEqual(len(df), 1)
        self.assertIn("variant", df.columns)
        self.assertIn("variant_scorer", df.columns)
        self.assertIn("raw_score", df.columns)
        self.assertIn("quantile_score", df.columns)
        self.assertIn("gene_name", df.columns)
        self.assertIn("gene_id", df.columns)
        self.assertIn("ontology_curie", df.columns)

        row = df.iloc[0]
        self.assertEqual(row["variant"], "chr22:36201698:A>C")
        self.assertEqual(row["variant_scorer"], "GeneMaskLFCScorer")
        self.assertAlmostEqual(row["raw_score"], 1.85, places=3)
        self.assertAlmostEqual(row["quantile_score"], 0.992, places=3)
        self.assertEqual(row["gene_name"], "APOL4")
        self.assertEqual(row["gene_id"], "ENSG00000100435")
        self.assertEqual(row["ontology_curie"], "UBERON:0001157")

    def test_read_variants_from_tsv(self):
        """TSV formatından varyant okuma testi."""
        tsv_file = self.test_path / "test_variants.tsv"
        tsv_file.write_text("chr\tpos\tref\talt\nchr22\t36201698\tA\tC\nchr1\t1000\tG\tT\n")

        variants = read_variants_from_tsv_or_csv(tsv_file)
        self.assertEqual(len(variants), 2)
        self.assertEqual(variants[0].chromosome, "chr22")
        self.assertEqual(variants[0].position, 36201698)
        self.assertEqual(variants[0].reference_bases, "A")
        self.assertEqual(variants[0].alternate_bases, "C")

    def test_read_variants_from_vcf(self):
        """VCF formatından varyant okuma testi."""
        vcf_file = self.test_path / "test.vcf"
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr22\t36201698\t.\tA\tC,G\t.\tPASS\t.\n"
            "chr1\t5000\t.\tT\t.\t.\tPASS\t.\n"  # ALT '.' atlanmalı
        )
        vcf_file.write_text(vcf_content)

        variants = read_variants_from_vcf(vcf_file)
        self.assertEqual(len(variants), 2)  # C and G
        self.assertEqual(variants[0].chromosome, "chr22")
        self.assertEqual(variants[0].alternate_bases, "C")
        self.assertEqual(variants[1].alternate_bases, "G")

    def test_cached_atlas_client(self):
        """CachedAtlasClient yerel disk önbelleği testi."""
        mock_raw_client = MagicMock()
        mock_raw_client.scorer_metadata.return_value = {"GeneMaskLFCScorer": MagicMock()}
        mock_raw_client.query_variant.return_value = {"dummy": None}

        cache_dir = self.test_path / "atlas_cache"
        client = CachedAtlasClient(mock_raw_client, cache_dir=cache_dir)

        var = genome.Variant(chromosome="chr22", position=36201698, reference_bases="A", alternate_bases="C")

        # 1. çağrı: mock client'a gitmeli
        res1 = client.query_variant(var)
        self.assertEqual(mock_raw_client.query_variant.call_count, 1)

        # 2. çağrı: önbellekten gelmeli, mock çağrılmamalı
        res2 = client.query_variant(var)
        self.assertEqual(mock_raw_client.query_variant.call_count, 1)

    def test_parser_atlas_commands(self):
        """CLI argüman ayrıştırıcısının atlas ve ui alt komutlarını doğrulama."""
        parser = build_parser()

        # atlas variant
        args = parser.parse_args(["atlas", "variant", "chr22:36201698:A>C", "--tissue", "colon"])
        self.assertEqual(args.command, "atlas")
        self.assertEqual(args.atlas_command, "variant")
        self.assertEqual(args.variant, "chr22:36201698:A>C")
        self.assertEqual(args.tissue, "colon")

        # variant --atlas
        args_var = parser.parse_args(["variant", "chr22:36201698:A>C", "--atlas"])
        self.assertTrue(args_var.atlas)

        # ui
        args_ui = parser.parse_args(["ui", "--port", "8080"])
        self.assertEqual(args_ui.command, "ui")
        self.assertEqual(args_ui.port, 8080)

    def test_web_app_endpoints(self):
        """Web UI sunucu uç noktalarının (status, scorers, static files) testi."""
        from http.server import ThreadingHTTPServer
        import urllib.request
        import threading
        from web_app import AlphaGenomeWebHandler

        server = ThreadingHTTPServer(("127.0.0.1", 0), AlphaGenomeWebHandler)
        port = server.server_port
        thread = threading.Thread(target=server.serve_forever)
        thread.daemon = True
        thread.start()

        try:
            # 1. /api/status
            url_status = f"http://127.0.0.1:{port}/api/status"
            with urllib.request.urlopen(url_status) as resp:
                self.assertEqual(resp.status, 200)
                data = json.loads(resp.read().decode("utf-8"))
                self.assertIn("has_key", data)
                self.assertIn("version", data)

            # 2. /api/scorers
            url_scorers = f"http://127.0.0.1:{port}/api/scorers"
            with urllib.request.urlopen(url_scorers) as resp:
                self.assertEqual(resp.status, 200)
                data = json.loads(resp.read().decode("utf-8"))
                self.assertIn("scorers", data)
                self.assertTrue(len(data["scorers"]) > 0)

            # 3. Static file / (index.html)
            url_home = f"http://127.0.0.1:{port}/"
            with urllib.request.urlopen(url_home) as resp:
                self.assertEqual(resp.status, 200)
                html = resp.read().decode("utf-8")
                self.assertIn("AlphaGenome", html)
                self.assertIn("dropZone", html)

        finally:
            server.shutdown()
            server.server_close()

    def test_parse_variant_and_interval(self):
        """Varyant ve genomik aralık string ayrıştırma fonksiyonlarının doğrulanması."""
        # Standart varyant
        v1 = parse_variant("chr10:89010000:G>A")
        self.assertEqual(v1.chromosome, "chr10")
        self.assertEqual(v1.position, 89010000)
        self.assertEqual(v1.reference_bases, "G")
        self.assertEqual(v1.alternate_bases, "A")

        # İki noktalı ve önreksiz varyant
        v2 = parse_variant("22:36201698:A:C")
        self.assertEqual(v2.chromosome, "chr22")
        self.assertEqual(v2.position, 36201698)
        self.assertEqual(v2.reference_bases, "A")
        self.assertEqual(v2.alternate_bases, "C")

        # Aralık ayrıştırma
        iv = parse_interval("chr22:35677410-36725986")
        self.assertEqual(iv.chromosome, "chr22")
        self.assertEqual(iv.start, 35677410)
        self.assertEqual(iv.end, 36725986)

        # Varyant pencere aralığı
        v_iv = get_variant_interval(v1, window_size=1000)
        self.assertEqual(v_iv.chromosome, "chr10")
        self.assertEqual(v_iv.start, 89010000 - 500)
        self.assertEqual(v_iv.end, 89010000 + 500)

        # Doku çözme
        t1 = resolve_tissue_term("UBERON:0001157")
        self.assertEqual(t1, ["UBERON:0001157"])

        # Hatalı varyant ValueError kontrolü
        with self.assertRaises(ValueError):
            parse_variant("invalid_variant_format")


if __name__ == "__main__":
    import json
    unittest.main()

