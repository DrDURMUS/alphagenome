#!/usr/bin/env python3
"""
AlphaGenome CLI — Kullanıcı dostu komut satırı aracı
=====================================================

Kullanım:
    conda activate alphagenome
    python alphagenome_cli.py --help

Alt Komutlar:
    variant    Tek bir varyant analizi (model tahmini veya --atlas)
    vcf        VCF dosyasından toplu varyant analizi
    region     Genomik bölge için ileri tahmin
    list       TSV/CSV dosyasından toplu varyant analizi
    tissues    Mevcut doku/ontoloji terimlerini ara ve listele
    visualize  Sonuçları görselleştir (PNG/PDF/SVG)
    atlas      AlphaGenome Atlas önceden hesaplanmış 9 milyar SNV veritabanı sorguları

Örnekler:
    python alphagenome_cli.py variant chr22:36201698:A>C --tissue colon
    python alphagenome_cli.py variant chr22:36201698:A>C --atlas
    python alphagenome_cli.py atlas variant chr22:36201698:A>C --tissue colon
    python alphagenome_cli.py atlas interval chr22:35677410-35678410 --tissue liver
    python alphagenome_cli.py atlas list variants.tsv --output atlas_scores.tsv
    python alphagenome_cli.py atlas scorers
    python alphagenome_cli.py region chr22:35677410-36725986 --tissue liver
    python alphagenome_cli.py visualize region chr22:35677410-36725986 --tissue colon
    python alphagenome_cli.py visualize variant chr22:36201698:A>C --tissue colon
    python alphagenome_cli.py tissues --search liver
"""

import argparse
import csv
import hashlib
import os
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from alphagenome.atlas import atlas as atlas_api
from alphagenome.data import genome, gene_annotation, transcript
from alphagenome.models import dna_client
from alphagenome.models.dna_model import Organism
from alphagenome.models.variant_scorers import (
    get_recommended_scorers,
    tidy_scores,
)
from alphagenome.visualization import plot_components
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for saving
import matplotlib.pyplot as plt

# ─── Constants ──────────────────────────────────────────────────────────────
DEFAULT_API_KEY = None  # Anahtar artık koda gömülmez — --api-key veya ALPHAGENOME_API_KEY kullanın
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "alphagenome"
GENCODE_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)
REF_ALT_COLORS = {"REF": "dimgrey", "ALT": "red"}

OUTPUT_TYPE_MAP = {
    "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
    "CAGE": dna_client.OutputType.CAGE,
    "ATAC": dna_client.OutputType.ATAC,
    "DNASE": dna_client.OutputType.DNASE,
    "PROCAP": dna_client.OutputType.PROCAP,
    "CHIP_HISTONE": dna_client.OutputType.CHIP_HISTONE,
    "CHIP_TF": dna_client.OutputType.CHIP_TF,
    "SPLICE_SITES": dna_client.OutputType.SPLICE_SITES,
    "SPLICE_JUNCTIONS": dna_client.OutputType.SPLICE_JUNCTIONS,
    "SPLICE_SITE_USAGE": dna_client.OutputType.SPLICE_SITE_USAGE,
    "CONTACT_MAPS": dna_client.OutputType.CONTACT_MAPS,
}

# ANSI renk kodları
GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
CYAN = "\033[96m"
BOLD = "\033[1m"
RESET = "\033[0m"
DIM = "\033[2m"


def ok(msg):
    print(f"{GREEN}✅ {msg}{RESET}")


def warn(msg):
    print(f"{YELLOW}⚠️  {msg}{RESET}")


def err(msg):
    print(f"{RED}❌ {msg}{RESET}")


def info(msg):
    print(f"{CYAN}ℹ️  {msg}{RESET}")


def header(msg):
    width = 60
    print(f"\n{BOLD}{'═' * width}")
    print(f"  {msg}")
    print(f"{'═' * width}{RESET}")


# ─── Cache ──────────────────────────────────────────────────────────────────
def _cache_key(*parts):
    """Argümanlardan SHA-256 cache anahtarı üretir."""
    key_str = "|".join(str(p) for p in parts)
    return hashlib.sha256(key_str.encode()).hexdigest()


class CachedModel:
    """
    AlphaGenome model istemcisini disk önbelleğiyle sarmalar.
    Aynı interval + tissue + modality kombinasyonu tekrar istendiğinde
    API çağrısı yapmadan önbellekten döner.
    """

    def __init__(self, model, cache_dir=DEFAULT_CACHE_DIR, enabled=True):
        self._model = model
        self._cache_dir = Path(cache_dir)
        self._enabled = enabled
        if enabled:
            self._cache_dir.mkdir(parents=True, exist_ok=True)

    def _load(self, key):
        path = self._cache_dir / f"{key}.pkl"
        if path.exists():
            try:
                with open(path, "rb") as f:
                    return True, pickle.load(f)
            except Exception:
                path.unlink(missing_ok=True)  # bozuk dosyayı sil
        return False, None

    def _save(self, key, value):
        path = self._cache_dir / f"{key}.pkl"
        try:
            with open(path, "wb") as f:
                pickle.dump(value, f)
        except Exception as e:
            warn(f"Önbellek yazılamadı: {e}")

    def _terms_key(self, ontology_terms):
        return tuple(sorted(ontology_terms or []))

    def _outputs_key(self, requested_outputs):
        return tuple(sorted(str(o) for o in (requested_outputs or [])))

    def predict_interval(self, interval, ontology_terms=None, requested_outputs=None):
        if not self._enabled:
            return self._model.predict_interval(
                interval, ontology_terms=ontology_terms, requested_outputs=requested_outputs
            )
        key = _cache_key(
            "predict_interval",
            interval.chromosome, interval.start, interval.end,
            self._terms_key(ontology_terms),
            self._outputs_key(requested_outputs),
        )
        hit, result = self._load(key)
        if hit:
            info(f"Önbellekten yüklendi {DIM}(predict_interval · {key[:8]}…){RESET}")
            return result
        result = self._model.predict_interval(
            interval, ontology_terms=ontology_terms, requested_outputs=requested_outputs
        )
        self._save(key, result)
        return result

    def predict_variant(self, interval, variant, ontology_terms=None, requested_outputs=None):
        if not self._enabled:
            return self._model.predict_variant(
                interval, variant, ontology_terms=ontology_terms, requested_outputs=requested_outputs
            )
        key = _cache_key(
            "predict_variant",
            interval.chromosome, interval.start, interval.end,
            variant.chromosome, variant.position,
            variant.reference_bases, variant.alternate_bases,
            self._terms_key(ontology_terms),
            self._outputs_key(requested_outputs),
        )
        hit, result = self._load(key)
        if hit:
            info(f"Önbellekten yüklendi {DIM}(predict_variant · {key[:8]}…){RESET}")
            return result
        result = self._model.predict_variant(
            interval, variant, ontology_terms=ontology_terms, requested_outputs=requested_outputs
        )
        self._save(key, result)
        return result

    def score_variant(self, interval, variant, variant_scorers):
        if not self._enabled:
            return self._model.score_variant(
                interval, variant, variant_scorers=variant_scorers
            )
        key = _cache_key(
            "score_variant",
            interval.chromosome, interval.start, interval.end,
            variant.chromosome, variant.position,
            variant.reference_bases, variant.alternate_bases,
            tuple(sorted(type(s).__name__ for s in variant_scorers)),
        )
        hit, result = self._load(key)
        if hit:
            info(f"Önbellekten yüklendi {DIM}(score_variant · {key[:8]}…){RESET}")
            return result
        result = self._model.score_variant(
            interval, variant, variant_scorers=variant_scorers
        )
        self._save(key, result)
        return result

    def score_variants(self, intervals, variants, variant_scorers, **kwargs):
        # Toplu işlemde her varyantı ayrı ayrı önbellekli score_variant üzerinden çalıştırır.
        scorer_names = tuple(sorted(type(s).__name__ for s in variant_scorers))
        results = []
        cache_hits = 0
        for interval, variant in zip(intervals, variants):
            key = _cache_key(
                "score_variant",
                interval.chromosome, interval.start, interval.end,
                variant.chromosome, variant.position,
                variant.reference_bases, variant.alternate_bases,
                scorer_names,
            )
            hit, cached = self._load(key)
            if hit:
                results.append(cached)
                cache_hits += 1
            else:
                results.append(None)  # placeholder

        missing_intervals = [intervals[i] for i, r in enumerate(results) if r is None]
        missing_variants = [variants[i] for i, r in enumerate(results) if r is None]

        if missing_intervals:
            if cache_hits:
                info(f"{cache_hits}/{len(variants)} varyant önbellekten yüklendi, "
                     f"{len(missing_intervals)} tanesi API'den alınıyor...")
            api_results = self._model.score_variants(
                missing_intervals, missing_variants,
                variant_scorers=variant_scorers, **kwargs
            )
            api_iter = iter(api_results)
            mi_idx = 0
            for i, r in enumerate(results):
                if r is None:
                    api_val = next(api_iter)
                    results[i] = api_val
                    key = _cache_key(
                        "score_variant",
                        missing_intervals[mi_idx].chromosome,
                        missing_intervals[mi_idx].start,
                        missing_intervals[mi_idx].end,
                        missing_variants[mi_idx].chromosome,
                        missing_variants[mi_idx].position,
                        missing_variants[mi_idx].reference_bases,
                        missing_variants[mi_idx].alternate_bases,
                        scorer_names,
                    )
                    self._save(key, api_val)
                    mi_idx += 1
        else:
            info(f"Tüm {len(variants)} varyant önbellekten yüklendi — API çağrısı yapılmadı!")

        return results

    def output_metadata(self):
        if not self._enabled:
            return self._model.output_metadata()
        key = "output_metadata"
        hit, result = self._load(key)
        if hit:
            return result
        result = self._model.output_metadata()
        self._save(key, result)
        return result


class CachedAtlasClient:
    """
    AlphaGenome Atlas istemcisini yerel disk önbelleğiyle sarmalar.
    Önceden hesaplanmış 9 milyar SNV veritabanından çekilen varyant ve aralık
    skorlarını önbelleğe alarak tekrarlanan çağrıları hızlandırır.
    """

    def __init__(self, client, cache_dir=DEFAULT_CACHE_DIR, enabled=True):
        self._client = client
        self._cache_dir = Path(cache_dir)
        self._enabled = enabled
        if enabled:
            self._cache_dir.mkdir(parents=True, exist_ok=True)

    def _load(self, key):
        path = self._cache_dir / f"atlas_{key}.pkl"
        if path.exists():
            try:
                with open(path, "rb") as f:
                    return True, pickle.load(f)
            except Exception:
                path.unlink(missing_ok=True)
        return False, None

    def _save(self, key, value):
        path = self._cache_dir / f"atlas_{key}.pkl"
        try:
            with open(path, "wb") as f:
                pickle.dump(value, f)
        except Exception as e:
            warn(f"Atlas önbellek yazılamadı: {e}")

    def scorer_metadata(self):
        if not self._enabled:
            return self._client.scorer_metadata()
        key = "scorer_metadata"
        hit, result = self._load(key)
        if hit:
            return result
        result = self._client.scorer_metadata()
        self._save(key, result)
        return result

    def query_variant(self, variant, *, requested_scorers=None, ontology_terms=None, gene_ids=None, gene_names=None):
        if not self._enabled:
            return self._client.query_variant(
                variant, requested_scorers=requested_scorers,
                ontology_terms=ontology_terms, gene_ids=gene_ids, gene_names=gene_names
            )
        key = _cache_key(
            "variant",
            variant.chromosome, variant.position, variant.reference_bases, variant.alternate_bases,
            tuple(sorted(requested_scorers or [])),
            tuple(sorted(str(o) for o in (ontology_terms or []))),
            tuple(sorted(gene_ids or [])),
            tuple(sorted(gene_names or [])),
        )
        hit, result = self._load(key)
        if hit:
            info(f"Önbellekten yüklendi {DIM}(Atlas variant · {key[:8]}…){RESET}")
            return result
        result = self._client.query_variant(
            variant, requested_scorers=requested_scorers,
            ontology_terms=ontology_terms, gene_ids=gene_ids, gene_names=gene_names
        )
        self._save(key, result)
        return result

    def query_variants(self, variants, *, requested_scorers=None, ontology_terms=None, gene_ids=None, gene_names=None, progress_bar=True, max_workers=10):
        if not self._enabled:
            return self._client.query_variants(
                variants, requested_scorers=requested_scorers,
                ontology_terms=ontology_terms, gene_ids=gene_ids, gene_names=gene_names,
                progress_bar=progress_bar, max_workers=max_workers
            )
        batch_key = _cache_key(
            "batch",
            len(variants),
            tuple((v.chromosome, v.position, v.reference_bases, v.alternate_bases) for v in variants),
            tuple(sorted(requested_scorers or [])),
            tuple(sorted(str(o) for o in (ontology_terms or []))),
            tuple(sorted(gene_ids or [])),
            tuple(sorted(gene_names or [])),
        )
        hit, result = self._load(batch_key)
        if hit:
            info(f"Önbellekten yüklendi {DIM}(Atlas batch · {len(variants)} varyant · {batch_key[:8]}…){RESET}")
            return result
        result = self._client.query_variants(
            variants, requested_scorers=requested_scorers,
            ontology_terms=ontology_terms, gene_ids=gene_ids, gene_names=gene_names,
            progress_bar=progress_bar, max_workers=max_workers
        )
        self._save(batch_key, result)
        return result

    def query_interval(self, interval, *, requested_scorers=None, ontology_terms=None, gene_ids=None, gene_names=None, progress_bar=True, max_workers=10):
        if not self._enabled:
            return self._client.query_interval(
                interval, requested_scorers=requested_scorers,
                ontology_terms=ontology_terms, gene_ids=gene_ids, gene_names=gene_names,
                progress_bar=progress_bar, max_workers=max_workers
            )
        key = _cache_key(
            "interval",
            interval.chromosome, interval.start, interval.end,
            tuple(sorted(requested_scorers or [])),
            tuple(sorted(str(o) for o in (ontology_terms or []))),
            tuple(sorted(gene_ids or [])),
            tuple(sorted(gene_names or [])),
        )
        hit, result = self._load(key)
        if hit:
            info(f"Önbellekten yüklendi {DIM}(Atlas interval · {key[:8]}…){RESET}")
            return result
        result = self._client.query_interval(
            interval, requested_scorers=requested_scorers,
            ontology_terms=ontology_terms, gene_ids=gene_ids, gene_names=gene_names,
            progress_bar=progress_bar, max_workers=max_workers
        )
        self._save(key, result)
        return result


# ─── Helpers ────────────────────────────────────────────────────────────────
def get_api_key(args, required=True):
    """API anahtarını argümanlardan veya env değişkeninden al."""
    key = getattr(args, "api_key", None)
    if key:
        return key
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if key:
        return key
    if DEFAULT_API_KEY:
        return DEFAULT_API_KEY
    if not required:
        return None
    err("API anahtarı bulunamadı.")
    info("Şunlardan birini kullanın:")
    info("  export ALPHAGENOME_API_KEY=<anahtarınız>   (önerilen)")
    info("  veya komuta  --api-key <anahtarınız>  ekleyin")
    raise RuntimeError("AlphaGenome API anahtarı bulunamadı.")


def create_model(args):
    """Model istemcisini oluştur ve önbelleğe sar."""
    api_key = get_api_key(args)
    info(f"API anahtarı ile bağlanıyor... (***{api_key[-4:]})")
    model = dna_client.create(api_key)
    ok("Model bağlantısı başarılı!")

    no_cache = getattr(args, "no_cache", False)
    cache_dir = getattr(args, "cache_dir", None) or DEFAULT_CACHE_DIR
    if no_cache:
        info("Önbellekleme devre dışı (--no-cache)")
        return CachedModel(model, enabled=False)
    info(f"Önbellek dizini: {DIM}{cache_dir}{RESET}")
    return CachedModel(model, cache_dir=cache_dir, enabled=True)


def create_atlas(args):
    """AlphaGenome Atlas istemcisini oluştur ve önbelleğe sar."""
    api_key = get_api_key(args)
    info(f"AlphaGenome Atlas'a bağlanılıyor... (***{api_key[-4:]})")
    client = atlas_api.create(api_key)
    ok("AlphaGenome Atlas bağlantısı başarılı!")

    no_cache = getattr(args, "no_cache", False)
    cache_dir = getattr(args, "cache_dir", None) or DEFAULT_CACHE_DIR
    if no_cache:
        info("Atlas önbellekleme devre dışı (--no-cache)")
        return CachedAtlasClient(client, enabled=False)
    info(f"Atlas önbellek dizini: {DIM}{cache_dir}{RESET}")
    return CachedAtlasClient(client, cache_dir=cache_dir, enabled=True)


def resolve_modalities(modality_names):
    """Modality isimlerini OutputType'a çevir."""
    output_types = []
    for name in modality_names:
        name_upper = name.upper()
        if name_upper not in OUTPUT_TYPE_MAP:
            err(f"Bilinmeyen modality: {name}")
            info(f"Geçerli seçenekler: {', '.join(OUTPUT_TYPE_MAP.keys())}")
            sys.exit(1)
        output_types.append(OUTPUT_TYPE_MAP[name_upper])
    return output_types


def parse_variant(variant_str):
    """
    Varyant string'ini parse et ve genome.Variant döndür.
    Desteklenen formatlar:
        chr10:89010000:G>A
        chr10:89010000:G:A
        chr10:89010000:G/A
        chr10-89010000-G-A
        chr10 89010000 G A
        10:89010000:G>A
    Hata durumunda ValueError fırlatır.
    """
    if not variant_str:
        raise ValueError("Varyant belirtilmedi")
    s = str(variant_str).strip()
    match = re.match(
        r"(?:chr)?(\w+)[:\-_ \t]+(\d+)[:\-_ \t]+([ACGTN]+)[>:/ \t\-_]+([ACGTN]+)",
        s,
        re.IGNORECASE,
    )
    if not match:
        raise ValueError(
            f"Geçersiz varyant formatı: '{variant_str}'. Beklenen format: chr10:89010000:G>A veya 10:89010000:G:A"
        )

    chrom_raw = match.group(1)
    chrom = f"chr{chrom_raw}" if not chrom_raw.lower().startswith("chr") else chrom_raw
    pos = int(match.group(2))
    ref = match.group(3).upper()
    alt = match.group(4).upper()

    return genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt,
    )


def parse_variant_str(variant_str):
    """CLI komut satırı için varyant string'ini parse et; hata durumunda sys.exit(1)."""
    try:
        return parse_variant(variant_str)
    except ValueError as e:
        err(str(e))
        info("Beklenen format: chr22:36201698:A>C  veya  22:36201698:A:C")
        sys.exit(1)


def parse_interval(interval_str):
    """
    Genomik bölge/aralık string'ini parse et ve genome.Interval döndür.
    Desteklenen formatlar:
        chr22:35677410-36725986
        22:35677410-36725986
        chr22 35677410 36725986
    Hata durumunda ValueError fırlatır.
    """
    if not interval_str:
        raise ValueError("Genomik aralık belirtilmedi")
    s = str(interval_str).strip()
    match = re.match(r"(?:chr)?(\w+)[:\-_ \t]+(\d+)[\-_ \t:]+(\d+)", s, re.IGNORECASE)
    if not match:
        raise ValueError(
            f"Geçersiz aralık formatı: '{interval_str}'. Beklenen format: chr22:35677410-36725986"
        )

    chrom_raw = match.group(1)
    chrom = f"chr{chrom_raw}" if not chrom_raw.lower().startswith("chr") else chrom_raw
    start = int(match.group(2))
    end = int(match.group(3))

    return genome.Interval(chromosome=chrom, start=start, end=end)


def parse_region_str(region_str):
    """CLI komut satırı için bölge parse et; hata durumunda sys.exit(1)."""
    try:
        return parse_interval(region_str)
    except ValueError as e:
        err(str(e))
        info("Beklenen format: chr22:35677410-36725986")
        sys.exit(1)


def get_variant_interval(variant, window_size=1_048_576):
    """Varyant etrafında pencere aralığı döndürür (auto_interval_for_variant alias)."""
    return auto_interval_for_variant(variant, window_size=window_size)


# AlphaGenome metadata'sında geçen ontoloji CURIE önekleri.
# (UBERON: anatomik yapı, CL: hücre tipi, EFO/CLO: hücre hattı, NTR: yeni terim)
ONTOLOGY_PREFIXES = ("UBERON:", "CL:", "EFO:", "CLO:", "NTR:")
# Doku-bağımsız (ontoloji boyutu olmayan) modaliteler — doku filtresi uygulanmaz.
TISSUE_AGNOSTIC_MODALITIES = {"SPLICE_SITES", "SPLICE_JUNCTIONS", "CONTACT_MAPS"}


def _looks_like_curie(s):
    """'EFO:0002067' gibi bir ontoloji terimi mi? (PREFIX:rakamlar)"""
    return bool(re.match(r"^[A-Z]+:\d+$", s.strip()))


def warn_missing_modality_terms(model, ontology_terms, modality_names):
    """Çözülen terimlerin istenen her modalitede mevcut olup olmadığını bildir.

    Bir terim yalnızca RNA_SEQ'te bulunup başka bir modalitede yoksa, o
    modalite için sessizce 0 track döner — burada bunu kullanıcıya gösteririz.
    """
    if not ontology_terms or not modality_names:
        return
    meta = model.output_metadata()
    requested = set()
    for term in ontology_terms:
        for name in modality_names:
            name_up = name.upper()
            if name_up in TISSUE_AGNOSTIC_MODALITIES:
                continue  # doku filtresi anlamsız
            md = getattr(meta, name_up.lower(), None)
            cols = getattr(md, "columns", [])
            if md is None or "ontology_curie" not in cols:
                continue
            requested.add((name_up, term))
            avail = set(md["ontology_curie"].dropna().unique())
            if term not in avail:
                warn(f"{name_up}: '{term}' bu modalitede yok → 0 track döner")


def find_ontology_terms(model, tissue_query):
    """
    Doku adına göre eşleşen TÜM ontoloji terimlerini bul.
    Büyük/küçük harf duyarsız arama yapar.
    Tüm eşleşen terimler API'ye gönderilir → farklı track'ler olarak görünür.
    """
    # Virgülle ayrılmış birden çok CURIE doğrudan verilebilir
    # (ör. "EFO:0001187,EFO:0002067"). Her parçayı ayrı terim olarak doğrula.
    parts = [p.strip() for p in tissue_query.split(",") if p.strip()]
    if parts and all(p.startswith(ONTOLOGY_PREFIXES) or _looks_like_curie(p)
                     for p in parts):
        invalid = [p for p in parts if not _looks_like_curie(p)]
        if invalid:
            err(f"Geçersiz ontoloji terimi: {', '.join(invalid)}")
            info("CURIE biçimi: <tip>:<id>  (ör. EFO:0002067)")
            sys.exit(1)
        ok(f"{len(parts)} doğrudan ontoloji terimi: {', '.join(parts)}")
        return parts

    meta = model.output_metadata()
    rna_meta = meta.rna_seq
    query = tissue_query.lower()

    # Biosample name'de arama
    matches = rna_meta[
        rna_meta["biosample_name"].str.lower().str.contains(query, na=False)
        | rna_meta["name"].str.lower().str.contains(query, na=False)
    ]

    if matches.empty:
        # gtex_tissue'da dene
        matches = rna_meta[
            rna_meta["gtex_tissue"].str.lower().str.contains(query, na=False)
        ]

    if matches.empty:
        warn(f"'{tissue_query}' ile eşleşen doku bulunamadı.")
        info("Mevcut dokuları listelemek için: python alphagenome_cli.py tissues --search <arama>")
        sys.exit(1)

    # Benzersiz ontoloji terimlerini al — hepsini kullan
    unique_terms = list(matches["ontology_curie"].unique())
    term_names = matches.drop_duplicates("ontology_curie")[["ontology_curie", "biosample_name"]]
    ok(f"{len(unique_terms)} doku eşleşmesi bulundu — tümü kullanılacak:")
    for _, row in term_names.iterrows():
        info(f"  → {row['biosample_name']} ({row['ontology_curie']})")

    return unique_terms


def auto_interval_for_variant(variant, window_size=1_048_576):
    """Varyant etrafında otomatik 1Mb pencere oluştur."""
    half = window_size // 2
    start = max(0, variant.position - half)
    end = start + window_size
    return genome.Interval(
        chromosome=variant.chromosome,
        start=start,
        end=end,
    )


def save_dataframe(df, output_path):
    """DataFrame'i dosyaya kaydet."""
    path = Path(output_path)
    if path.suffix == ".csv":
        df.to_csv(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False)
    ok(f"Sonuçlar kaydedildi: {path}")


def load_transcripts(interval, transcript_set="mane"):
    """GENCODE transkriptlerini yükle.

    transcript_set:
      "mane"    — MANE Select (klinik standart kanonik transkript). Bir gen için
                  MANE Select yoksa o gen bu sette görünmez.
      "longest" — gen başına en uzun izoform.
    """
    info("GENCODE gen anotasyonları yükleniyor...")
    try:
        gtf = pd.read_feather(GENCODE_URL)
        gtf_filtered = gene_annotation.filter_protein_coding(gtf)

        if transcript_set == "mane":
            mane = gene_annotation.filter_to_mane_select_transcript(gtf_filtered)
            extractor = transcript.TranscriptExtractor(mane)
            transcripts = extractor.extract(interval)
            if transcripts:
                ok(f"{len(transcripts)} MANE Select transkript bulundu")
                return transcripts
            # MANE Select bu bölgede yoksa en uzun izoforma düş
            warn("Bu bölgede MANE Select transkript yok — en uzun izoforma geçiliyor")

        gtf_filtered = gene_annotation.filter_to_longest_transcript(gtf_filtered)
        extractor = transcript.TranscriptExtractor(gtf_filtered)
        transcripts = extractor.extract(interval)
        if transcripts:
            ok(f"{len(transcripts)} transkript bulundu (en uzun izoform)")
        else:
            warn("Bu bölgede protein-kodlayan gen bulunamadı")
        return transcripts
    except Exception as e:
        warn(f"Gen anotasyonları yüklenemedi: {e}")
        return []


SUPPORTED_IMG_FORMATS = {".png", ".pdf", ".svg", ".svgz", ".jpg", ".jpeg", ".eps",
                         ".tif", ".tiff", ".webp", ".ps", ".pgf", ".raw", ".rgba"}


def save_figure(fig, output_path):
    """Matplotlib figürünü kaydet."""
    path = Path(output_path)
    if path.suffix.lower() not in SUPPORTED_IMG_FORMATS:
        warn(f"'{path.suffix}' görsel formatı değil — otomatik olarak .png ekleniyor")
        path = Path(str(path) + ".png")
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    ok(f"Görsel kaydedildi: {path}")
    plt.close(fig)


# ─── Sonuç yorumlama ─────────────────────────────────────────────────────────
# Scorer adını biyolojik mekanizmaya eşleyen anahtar kelimeler.
MECHANISM_KEYWORDS = {
    "Splicing (splice bölgesi/kullanımı/junction)": ("splice",),
    "Gen ifadesi (RNA/CAGE/PROCAP)": ("rna", "cage", "expression", "procap", "polyadenyl"),
    "Kromatin/düzenleyici (ATAC/DNase/ChIP/Contact)":
        ("atac", "dnase", "chip", "histone", "contact", "tf"),
}


def _mechanism_of(scorer_name):
    s = str(scorer_name).lower()
    for label, keys in MECHANISM_KEYWORDS.items():
        if any(k in s for k in keys):
            return label
    return "Diğer"


# Band seviyeleri (0..3) ve renkleri.
EFFECT_LEVELS = ["ihmal edilebilir", "zayıf", "orta", "BELİRGİN"]
LEVEL_COLORS = [GREEN, CYAN, YELLOW, RED]

# Anlamlılık eşikleri: |quantile_score| — arka plana göre ne kadar uç.
QUANTILE_THR = (0.90, 0.50, 0.20)
# Etki büyüklüğü eşikleri: |raw_score| — gen-düzeyi scorer'lar için log-fold-change
# veya olasılık farkı (örn. 0.1 ≈ %7 ifade değişimi; 0.5 ≈ 1.4×). HEURİSTİK.
RAW_THR = (0.50, 0.20, 0.10)


def _level(absval, thr):
    """Bir |değeri| 0..3 ordinal seviyeye çevirir."""
    if absval >= thr[0]:
        return 3
    if absval >= thr[1]:
        return 2
    if absval >= thr[2]:
        return 1
    return 0


def _combined_level(abs_quantile, abs_raw, has_raw=True):
    """İstatistiksel olarak doğru band: etki büyüklüğü VE anlamlılığın daha zayıfı.

    Bir etki ancak hem arka plana göre uç (quantile) HEM DE yeterli büyüklükte
    (raw) ise yüksek seviyeye çıkar. raw yoksa yalnız quantile seviyesi.
    """
    q_lvl = _level(abs_quantile, QUANTILE_THR)
    if not has_raw:
        return q_lvl
    return min(q_lvl, _level(abs_raw, RAW_THR))


def _band_from_level(level):
    return EFFECT_LEVELS[level], LEVEL_COLORS[level]


def _raw_interpretable(scorer_name, output_type):
    """raw_score'un birimi büyüklük yorumuna uygun mu?

    raw_score yalnızca AYNI scorer içinde kıyaslanabilir. Evrensel eşik (RAW_THR)
    sadece birimi bilinen scorer'lar için anlamlı:
      • '...LFC...'  → log-fold-change (ifade) → eşik uygulanır
      • SPLICE_SITES / SPLICE_SITE_USAGE → olasılık farkı [0–1] → eşik uygulanır
    GeneMaskActive, CenterMask, Contact vb. farklı/bilinmeyen birimlerde → uygulanmaz.
    """
    if "LFC" in str(scorer_name):
        return True
    # Yalnızca splice bölgesi OLASILIK skorları (SPLICE_SITES / SPLICE_SITE_USAGE)
    # SpliceAI Δ ile aynı [0–1] anlamdadır → 0.2/0.5 eşikleri uygulanabilir.
    # SPLICE_JUNCTIONS raw'ı junction KULLANIM metriğidir, KALİBRE DEĞİL; SpliceAI
    # ile aynı ölçekte değildir (over-call riski) → güvenilir sayılmaz, yalnız
    # quantile ile gösterilir ve 'zayıf' ile sınırlanır.
    if str(output_type) in ("SPLICE_SITES", "SPLICE_SITE_USAGE"):
        return True
    return False


def interpret_results(variant, score_df=None, results_rows=None, gene=None):
    """Skor/etki sonuçlarını insan-okunur bir yorum bloğu olarak yazdırır.

    gene verilirse yorum yalnızca o gene odaklanır (1 Mb penceredeki diğer
    genleri dışlar). Yalnızca KARAR DESTEK aracıdır; ACMG sınıflandırması yapmaz.
    """
    print(f"\n{BOLD}🧭 Yorum (otomatik özet — karar destek):{RESET}")

    # ── Kalibre skor tablosu varsa onu kullan (tercih edilen yol) ──
    if score_df is not None and not score_df.empty and (
            "quantile_score" in score_df.columns or "raw_score" in score_df.columns):
        # quantile_score kalibre ve modaliteler arası kıyaslanabilir → tercih et
        if "quantile_score" in score_df.columns:
            val_col, score_kind = "quantile_score", "quantile"
        else:
            val_col, score_kind = "raw_score", "raw"
        has_raw = "raw_score" in score_df.columns
        mech_col = "output_type" if "output_type" in score_df.columns else (
            "variant_scorer" if "variant_scorer" in score_df.columns else None)
        scorer_col = "variant_scorer" if "variant_scorer" in score_df.columns else mech_col
        gene_col = "gene_name" if "gene_name" in score_df.columns else None

        df = score_df.dropna(subset=[val_col]).copy()

        # Gene odağı (varsa)
        if gene and gene_col:
            mask = df[gene_col].astype(str).str.upper() == gene.upper()
            if not mask.any():  # sembol eşleşmezse 'içeriyor' ile dene
                mask = df[gene_col].astype(str).str.upper().str.contains(
                    gene.upper(), na=False)
            if mask.any():
                df = df[mask]
                ok(f"Yorum '{gene}' genine odaklandı ({len(df)} kayıt)")
            else:
                warn(f"'{gene}' skor tablosunda yok — tüm genler kullanılıyor")

        has_q = "quantile_score" in df.columns
        ot_col = "output_type" if "output_type" in df.columns else mech_col
        if df.empty:
            warn("Skor tablosunda sayısal değer yok.")
            return

        # Her mekanizma için temsilci satırı seç. KURAL:
        #  • raw birimi yorumlanabilir scorer'lar (LFC, splice olasılığı) tercih edilir
        #    → band = min(büyüklük, anlamlılık).
        #  • Yorumlanamayan raw (Active/CenterMask...) yalnız quantile verir ve
        #    büyüklük doğrulanamadığından band 'zayıf' ile sınırlanır.
        #  • Bir mekanizmada güvenilir scorer varsa, güvenilmez olanlar yok sayılır.
        per_mech = {}  # mech -> cand
        for _, row in df.iterrows():
            mech = _mechanism_of(row[mech_col]) if mech_col else "Diğer"
            q_signed = float(row["quantile_score"]) if has_q else float(row[val_col])
            raw_signed = float(row["raw_score"]) if has_raw else q_signed
            scorer = row[scorer_col] if scorer_col else "?"
            ot = row[ot_col] if ot_col else None
            trusted = bool(has_raw and has_q and _raw_interpretable(scorer, ot))
            if trusted:
                level = _combined_level(abs(q_signed), abs(raw_signed), True)
            else:
                level = min(_level(abs(q_signed), QUANTILE_THR), 1)  # büyüklük yok → ≤ zayıf
            cand = {
                "level": level, "q": abs(q_signed), "raw": abs(raw_signed),
                "q_signed": q_signed, "raw_signed": raw_signed,
                "scorer": scorer, "gene": row[gene_col] if gene_col else None,
                "trusted": trusted,
            }
            mag = abs(raw_signed) if trusted else abs(q_signed)
            key = (trusted, level, mag)  # trusted birincil → güvenilir scorer kazanır
            if mech not in per_mech:
                per_mech[mech] = cand
            else:
                ex = per_mech[mech]
                if key > (ex["trusted"], ex["level"],
                          ex["raw"] if ex["trusted"] else ex["q"]):
                    per_mech[mech] = cand

        scope = f"'{gene}' geni" if gene else "1 Mb penceredeki TÜM genler"
        info(f"Skor sütunu: raw_score (etki büyüklüğü) × quantile_score (anlamlılık) "
             f"· kapsam: {scope} · {len(df)} kayıt")
        print(f"  {DIM}Band = büyüklük × anlamlılık; raw yalnızca birimi bilinen "
              f"scorer'larda (LFC/splice) büyüklük sayılır.{RESET}")

        for mech, c in sorted(per_mech.items(),
                              key=lambda kv: (kv[1]["level"], kv[1]["raw"]), reverse=True):
            band, color = _band_from_level(c["level"])
            yon = (" ↑ (artış)" if c["raw_signed"] > 0
                   else (" ↓ (azalış)" if c["raw_signed"] < 0 else ""))
            gtxt = f" · gen: {c['gene']}" if c["gene"] else ""
            qtxt = f"quantile={c['q']:.3f}"
            if c["trusted"]:
                detail = f"raw={c['raw_signed']:+.3f}, {qtxt}"
            else:
                detail = f"{qtxt}, raw birimi bilinmiyor → büyüklük doğrulanamadı"
            print(f"  {color}{band:>16}{RESET}  {mech}{yon}  "
                  f"{DIM}({detail}, {c['scorer']}{gtxt}){RESET}")

        top_mech, tc = max(per_mech.items(),
                           key=lambda kv: (kv[1]["level"], kv[1]["trusted"], kv[1]["q"]))
        band, color = _band_from_level(tc["level"])
        print()
        if tc["level"] == 0:
            print(f"  {GREEN}➜ AlphaGenome bu varyant için kayda değer bir fonksiyonel "
                  f"etki ÖNGÖRMÜYOR{RESET} {DIM}(etki büyüklüğü ve/veya anlamlılık düşük).{RESET}")
            print(f"    {DIM}ACMG çerçevesinde BS3 (fonksiyonel — zararsız) yönünde "
                  f"destekleyici sinyal olabilir.{RESET}")
        else:
            gtxt = f" ({tc['gene']} geninde)" if tc["gene"] else ""
            magtxt = (f"raw={tc['raw_signed']:+.3f}, quantile={tc['q']:.3f}"
                      if tc["trusted"] else
                      f"quantile={tc['q']:.3f} (büyüklük doğrulanamadı)")
            print(f"  {color}➜ Baskın öngörülen etki: {top_mech}{gtxt} — {band} "
                  f"({magtxt}).{RESET}")
            print(f"    {DIM}Güçlü ve zarar yönünde bir etkiyse ACMG'de PS3 "
                  f"(fonksiyonel — zararlı) yönünde destekleyici sinyal olabilir.{RESET}")

        # Splicing belirgin ama gen düzeyi ifade korunuyorsa: ekzon atlanması /
        # kriptik splice → lokal RNA-seq kapsaması yeniden dağılır (gözlemlenen düşüş).
        splice_lvl = max((c["level"] for m, c in per_mech.items() if "Splic" in m),
                         default=0)
        expr_lvl = max((c["level"] for m, c in per_mech.items() if "ifade" in m),
                       default=0)
        rna_local = next((float(r["max_abs_effect"]) for r in (results_rows or [])
                          if str(r.get("modality", "")).upper() == "RNA_SEQ"), None)
        if splice_lvl >= 2 and expr_lvl == 0:
            extra = (f" RNA-seq izinde lokal |fark| ~{rna_local:.1f} görülebilir."
                     if rna_local and rna_local >= 1.0 else "")
            print(f"    {YELLOW}ℹ Not: gen düzeyi ifade ~sabit ama splicing güçlü → "
                  f"ekzon atlanması/kriptik splice beklenir;{RESET}"
                  f"{DIM} spesifik ekzonlarda RNA-seq kapsaması düşüp başka bölgelerde "
                  f"artabilir (toplam transkript korunur).{extra}{RESET}")

        if not gene and gene_col:
            print(f"    {DIM}İpucu: ilgilendiğiniz gene odaklanmak için --gene <SEMBOL> "
                  f"ekleyin (pencerede {df[gene_col].nunique()} gen var).{RESET}")

    # ── --score yoksa: yalnızca ham iz farkları var ──
    elif results_rows:
        warn("(--score verilmedi) — aşağıdaki değerler HAM iz farklarıdır; "
             "ölçekler modaliteler arasında FARKLIDIR ve kıyaslanamaz.")
        info("Kalibre, kıyaslanabilir skorlar için komuta --score ekleyin.")
        for r in sorted(results_rows, key=lambda x: -float(x["max_abs_effect"])):
            mech = _mechanism_of(r["modality"])
            print(f"  {CYAN}{r['modality']:<18}{RESET} "
                  f"{DIM}max|fark|={float(r['max_abs_effect']):.4f}  "
                  f"ort|fark|={float(r['mean_abs_effect']):.6f}  ({mech}){RESET}")
        print(f"\n  {YELLOW}➜ Ham farklar mekanizma kararı için yeterli değil — "
              f"--score ile kalibre skor alın.{RESET}")
    else:
        warn("Yorumlanacak sayısal sonuç yok.")
        return

    # Uyarılar
    print(f"\n  {YELLOW}⚠ Uyarılar:{RESET}")
    print(f"    {DIM}• Skor büyüklüğü GÖRECELİDİR; klinik eşik değildir. Aynı gendeki "
          f"bilinen zararlı/zararsız varyantların skorlarıyla kıyaslayın.{RESET}")
    print(f"    {DIM}• 'Yön' (artış/azalış) ≠ patojenite; gen-hastalık mekanizmasına bağlıdır.{RESET}")
    print(f"    {DIM}• Splice tahminlerini SpliceAI ile DOĞRULAYIN — modeller farklı "
          f"sonuç verebilir; SpliceAI splice için klinik altın standarttır.{RESET}")
    print(f"    {DIM}• Doku seçimi sonucu etkiler; ilgili dokuyu seçtiğinizden emin olun.{RESET}")
    print(f"    {DIM}• Koordinatın hg38/GRCh38 olduğunu doğrulayın.{RESET}")
    print(f"    {DIM}• Bu özet karar desteğidir; tek başına ACMG sınıflandırması yapmaz.{RESET}")


# ─── AlphaGenome Atlas Helpers & Commands ──────────────────────────────────
def read_variants_from_tsv_or_csv(list_path):
    """TSV veya CSV dosyasından varyantları oku."""
    path = Path(list_path)
    if not path.exists():
        err(f"Dosya bulunamadı: {path}")
        sys.exit(1)

    info(f"Dosya okunuyor: {path}")
    sep = "," if path.suffix == ".csv" else "\t"

    variants = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=sep)
        header_row = next(reader, None)

        if header_row:
            first_val = header_row[0].lower().strip()
            if first_val not in ("chr", "chrom", "chromosome", "#chrom", "#chr"):
                # İlk satır doğrudan veri mi?
                if (first_val.startswith("chr") or first_val.isdigit()) and len(header_row) >= 4:
                    chrom = header_row[0] if header_row[0].startswith("chr") else f"chr{header_row[0]}"
                    try:
                        variants.append(
                            genome.Variant(
                                chromosome=chrom,
                                position=int(header_row[1]),
                                reference_bases=header_row[2].upper(),
                                alternate_bases=header_row[3].upper(),
                            )
                        )
                    except ValueError:
                        pass

        for row in reader:
            if len(row) < 4:
                continue
            if row[0].startswith("#"):
                continue
            chrom = row[0] if row[0].startswith("chr") else f"chr{row[0]}"
            try:
                pos = int(row[1])
            except ValueError:
                continue
            ref = row[2].upper()
            alt = row[3].upper()
            variants.append(
                genome.Variant(
                    chromosome=chrom,
                    position=pos,
                    reference_bases=ref,
                    alternate_bases=alt,
                )
            )

    if not variants:
        err("Dosyada varyant bulunamadı.")
        info("Beklenen format (TSV): chr\tpos\tref\talt")
        sys.exit(1)

    ok(f"{len(variants)} varyant okundu")
    return variants


def read_variants_from_vcf(vcf_path):
    """VCF dosyasından varyantları oku."""
    path = Path(vcf_path)
    if not path.exists():
        err(f"VCF dosyası bulunamadı: {path}")
        sys.exit(1)

    variants = []
    info(f"VCF okunuyor: {path}")
    opener = open
    if str(path).endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref = parts[3].upper()
            for alt in parts[4].upper().split(","):
                alt = alt.strip()
                if alt in (".", "*", ""):
                    continue
                variants.append(
                    genome.Variant(
                        chromosome=chrom,
                        position=pos,
                        reference_bases=ref,
                        alternate_bases=alt,
                    )
                )

    if not variants:
        err("VCF dosyasında varyant bulunamadı.")
        sys.exit(1)

    ok(f"{len(variants)} varyant okundu")
    return variants


def resolve_atlas_scorers(atlas_client, scorers_arg):
    """Atlas için istenen scorer'ları belirle."""
    if scorers_arg:
        if isinstance(scorers_arg, list):
            scorers = []
            for item in scorers_arg:
                scorers.extend([s.strip() for s in item.split(",") if s.strip()])
            return scorers
        return [s.strip() for s in scorers_arg.split(",") if s.strip()]

    # Kullanıcı belirtmediyse Atlas metadata'sından tüm mevcut scorer'ları sorgula
    try:
        metadata = atlas_client.scorer_metadata()
        if metadata:
            return list(metadata.keys())
    except Exception as e:
        warn(f"Atlas scorer metadata alınamadı ({e}). Varsayılan scorer listesi deneniyor.")

    return [
        "GeneMaskLFCScorer",
        "GeneMaskActiveScorer",
        "CenterMaskScorer",
        "SpliceJunctionScorer",
    ]


def resolve_atlas_ontology(args, tissue_query):
    """Atlas sorgusu için doku ontoloji terimlerini çözümle."""
    if not tissue_query:
        return None

    # Doğrudan CURIE formatı kontrolü (UBERON:0001157 vb.)
    parts = [p.strip() for p in tissue_query.split(",") if p.strip()]
    if parts and all(p.startswith(ONTOLOGY_PREFIXES) or _looks_like_curie(p) for p in parts):
        return parts

    # Model metadata üzerinden doku araması yap
    try:
        model = create_model(args)
        return find_ontology_terms(model, tissue_query)
    except Exception as e:
        warn(f"Doku metadatası taranamadı: {e}. Terim doğrudan kullanılıyor.")
        return [tissue_query]


def resolve_tissue_term(tissue_query, args=None):
    """
    Doku veya hücre tipi sorgusunu (örn. 'colon', 'liver', 'UBERON:0001157')
    ontoloji CURIE terim listesine çözümler.
    """
    if not tissue_query:
        return None
    if isinstance(tissue_query, (list, tuple)):
        return list(tissue_query)
    tissue_str = str(tissue_query).strip()
    if not tissue_str:
        return None

    # Doğrudan CURIE formatı kontrolü
    parts = [p.strip() for p in tissue_str.split(",") if p.strip()]
    if parts and all(p.startswith(ONTOLOGY_PREFIXES) or _looks_like_curie(p) for p in parts):
        return parts

    try:
        if args is None:
            class DummyArgs:
                api_key = None
                no_cache = False
                cache_dir = None
            args = DummyArgs()
        return resolve_atlas_ontology(args, tissue_str)
    except Exception:
        return [tissue_str]


def tidy_atlas_scores(scores_map):
    """
    AlphaGenome Atlas'tan dönen Mapping[str, AnnData] yapısını
    kullanıcı dostu, filtrelenebilir ve sıralanabilir bir pandas DataFrame'e dönüştürür.
    """
    if not scores_map:
        return pd.DataFrame()

    all_records = []
    for scorer_name, adata in scores_map.items():
        if adata is None:
            continue
        try:
            shape = adata.shape
        except Exception:
            continue
        if shape[0] == 0 or shape[1] == 0:
            continue

        n_obs, n_var = shape
        has_quantiles = hasattr(adata, "layers") and "quantiles" in adata.layers

        for i in range(n_obs):
            obs_row = adata.obs.iloc[i].to_dict() if (hasattr(adata, "obs") and not adata.obs.empty) else {}
            var_obj = obs_row.get("variant")
            if hasattr(var_obj, "chromosome"):
                var_str = f"{var_obj.chromosome}:{var_obj.position}:{var_obj.reference_bases}>{var_obj.alternate_bases}"
            elif var_obj is not None:
                var_str = str(var_obj)
            else:
                var_str = None

            gene_id = obs_row.get("gene_id")
            if gene_id and "." in str(gene_id):
                gene_id_clean = str(gene_id).split(".")[0]
            else:
                gene_id_clean = gene_id

            gene_name = obs_row.get("gene_name")
            gene_strand = obs_row.get("strand")

            for j in range(n_var):
                var_row = adata.var.iloc[j].to_dict() if (hasattr(adata, "var") and not adata.var.empty) else {}
                raw_score = float(adata.X[i, j])
                rec = {
                    "variant": var_str,
                    "variant_scorer": scorer_name,
                    "raw_score": raw_score,
                }
                if has_quantiles:
                    rec["quantile_score"] = float(adata.layers["quantiles"][i, j])
                if gene_name:
                    rec["gene_name"] = str(gene_name)
                if gene_id_clean:
                    rec["gene_id"] = str(gene_id_clean)
                if gene_strand:
                    rec["gene_strand"] = str(gene_strand)
                if "name" in var_row:
                    rec["track_name"] = str(var_row["name"])
                if "ontology_curie" in var_row:
                    rec["ontology_curie"] = str(var_row["ontology_curie"])
                if "strand" in var_row:
                    rec["track_strand"] = str(var_row["strand"])

                all_records.append(rec)

    if not all_records:
        return pd.DataFrame()
    return pd.DataFrame(all_records)


def cmd_atlas_variant(args):
    """AlphaGenome Atlas'tan tek bir varyantın önceden hesaplanmış skorlarını sorgula."""
    header("AlphaGenome Atlas — Varyant Skor Sorgulama")
    atlas = create_atlas(args)
    variant = parse_variant_str(args.variant)
    info(f"Varyant: {variant.chromosome}:{variant.position} "
         f"{variant.reference_bases}>{variant.alternate_bases}")

    ontology_terms = None
    if getattr(args, "tissue", None):
        ontology_terms = resolve_atlas_ontology(args, args.tissue)

    requested_scorers = resolve_atlas_scorers(atlas, getattr(args, "scorers", None))
    info(f"Kullanılan Scorer'lar: {', '.join(requested_scorers)}")

    gene_names = [args.gene] if getattr(args, "gene", None) else None

    print(f"\n{DIM}AlphaGenome Atlas sorgulanıyor (önceden hesaplanmış 9 milyar SNV veritabanı)...{RESET}")
    raw_scores = atlas.query_variant(
        variant,
        requested_scorers=requested_scorers,
        ontology_terms=ontology_terms,
        gene_names=gene_names,
    )
    ok("Atlas sorgusu tamamlandı!")

    df = tidy_atlas_scores(raw_scores)
    if df is not None and not df.empty:
        ok(f"Atlas skor tablosu: {df.shape[0]} satır × {df.shape[1]} sütun")
        sort_col = "quantile_score" if "quantile_score" in df.columns else "raw_score"

        top_df = df
        if getattr(args, "gene", None) and "gene_name" in df.columns:
            g = args.gene.upper()
            gmask = df["gene_name"].astype(str).str.upper() == g
            if not gmask.any():
                gmask = df["gene_name"].astype(str).str.upper().str.contains(g, na=False)
            if gmask.any():
                top_df = df[gmask]

        top = top_df.reindex(
            top_df[sort_col].abs().sort_values(ascending=False).index
        ).head(15)

        gtag = f" · gen: {args.gene}" if getattr(args, "gene", None) else ""
        print(f"\n{BOLD}🏆 Atlas En Yüksek Skorlar (|{sort_col}|{gtag}):{RESET}")
        display_cols = [c for c in
                        ["variant", "gene_name", "variant_scorer", "track_name",
                         "ontology_curie", "raw_score", "quantile_score"]
                        if c in top.columns]
        print(top[display_cols].to_string(index=False))

        if getattr(args, "output", None):
            save_dataframe(df, args.output)
    else:
        warn("Bu varyant ve kriterler için Atlas skoru bulunamadı.")

    # Otomatik yorum bloğu
    if not getattr(args, "no_interpret", False):
        interpret_results(variant, score_df=df, gene=getattr(args, "gene", None))

    ok("Atlas varyant analizi tamamlandı!")


def cmd_atlas_interval(args):
    """Genomik bir aralıktaki tüm varyantların Atlas skorlarını sorgula."""
    header("AlphaGenome Atlas — Bölge/Aralık Skor Sorgulama")
    atlas = create_atlas(args)
    interval = parse_region_str(args.interval)
    width = interval.end - interval.start
    info(f"Bölge: {interval.chromosome}:{interval.start}-{interval.end} ({width:,} bp)")
    if width > 50_000:
        warn(f"Bölge genişliği {width:,} bp (~{width * 3:,} varyant). Atlas sorgusu biraz zaman alabilir.")

    ontology_terms = None
    if getattr(args, "tissue", None):
        ontology_terms = resolve_atlas_ontology(args, args.tissue)

    requested_scorers = resolve_atlas_scorers(atlas, getattr(args, "scorers", None))
    info(f"Kullanılan Scorer'lar: {', '.join(requested_scorers)}")

    gene_names = [args.gene] if getattr(args, "gene", None) else None
    progress_bar = not getattr(args, "no_progress", False)
    workers = getattr(args, "workers", 10) or 10

    print(f"\n{DIM}Bölgedeki tüm tek nükleotid varyantları Atlas'tan çekiliyor...{RESET}")
    raw_scores = atlas.query_interval(
        interval,
        requested_scorers=requested_scorers,
        ontology_terms=ontology_terms,
        gene_names=gene_names,
        progress_bar=progress_bar,
        max_workers=workers,
    )
    ok("Atlas bölge sorgusu tamamlandı!")

    df = tidy_atlas_scores(raw_scores)
    if df is not None and not df.empty:
        ok(f"Toplam Atlas skoru: {df.shape[0]} satır × {df.shape[1]} sütun")
        sort_col = "quantile_score" if "quantile_score" in df.columns else "raw_score"
        top = df.reindex(
            df[sort_col].abs().sort_values(ascending=False).index
        ).head(15)

        print(f"\n{BOLD}🏆 Bölgedeki En Etkili Varyantlar (|{sort_col}|):{RESET}")
        display_cols = [c for c in
                        ["variant", "gene_name", "variant_scorer", "track_name",
                         "raw_score", "quantile_score"]
                        if c in top.columns]
        print(top[display_cols].to_string(index=False))

        out_path = getattr(args, "output", None) or f"atlas_{interval.chromosome}_{interval.start}_{interval.end}.tsv"
        save_dataframe(df, out_path)
    else:
        warn("Bu bölge için Atlas skoru döndürülmedi.")

    ok("Atlas bölge analizi tamamlandı!")


def cmd_atlas_list(args):
    """TSV/CSV dosyasındaki varyantları toplu olarak AlphaGenome Atlas'tan sorgula."""
    header("AlphaGenome Atlas — Dosyadan Toplu Varyant Analizi")
    atlas = create_atlas(args)

    variants = read_variants_from_tsv_or_csv(args.list_file)

    ontology_terms = None
    if getattr(args, "tissue", None):
        ontology_terms = resolve_atlas_ontology(args, args.tissue)

    requested_scorers = resolve_atlas_scorers(atlas, getattr(args, "scorers", None))
    info(f"Kullanılan Scorer'lar: {', '.join(requested_scorers)}")

    gene_names = [args.gene] if getattr(args, "gene", None) else None
    progress_bar = not getattr(args, "no_progress", False)
    workers = getattr(args, "workers", 10) or 10

    print(f"\n{DIM}{len(variants)} varyant AlphaGenome Atlas üzerinden sorgulanıyor...{RESET}")
    raw_scores = atlas.query_variants(
        variants,
        requested_scorers=requested_scorers,
        ontology_terms=ontology_terms,
        gene_names=gene_names,
        progress_bar=progress_bar,
        max_workers=workers,
    )
    ok("Atlas toplu sorgusu tamamlandı!")

    df = tidy_atlas_scores(raw_scores)
    if df is not None and not df.empty:
        ok(f"Toplam Atlas skoru: {df.shape[0]} satır × {df.shape[1]} sütun")
        sort_col = "quantile_score" if "quantile_score" in df.columns else "raw_score"
        print(f"\n{BOLD}📊 Özet:{RESET}")
        if "variant" in df.columns:
            print(f"  Benzersiz varyant sayısı: {df['variant'].nunique()}")
        if "gene_name" in df.columns:
            print(f"  Etkilenen gen sayısı: {df['gene_name'].nunique()}")
        print(f"  Max |skor|: {df[sort_col].abs().max():.6f}")

        top5 = df.reindex(df[sort_col].abs().sort_values(ascending=False).index).head(10)
        display_cols = [c for c in ["variant", "gene_name", "variant_scorer", "raw_score", "quantile_score"] if c in top5.columns]
        if display_cols:
            print(f"\n{BOLD}🏆 En yüksek skorlar (|{sort_col}|):{RESET}")
            print(top5[display_cols].to_string(index=False))

        output_file = getattr(args, "output", None) or (Path(args.list_file).stem + "_atlas_scores.tsv")
        save_dataframe(df, output_file)
    else:
        warn("Varyantlar için Atlas skoru döndürülmedi.")

    ok("Atlas toplu varyant analizi tamamlandı!")


def cmd_atlas_vcf(args):
    """VCF dosyasındaki varyantları AlphaGenome Atlas'tan toplu sorgula."""
    header("AlphaGenome Atlas — VCF Toplu Varyant Analizi")
    atlas = create_atlas(args)

    variants = read_variants_from_vcf(args.vcf_file)

    ontology_terms = None
    if getattr(args, "tissue", None):
        ontology_terms = resolve_atlas_ontology(args, args.tissue)

    requested_scorers = resolve_atlas_scorers(atlas, getattr(args, "scorers", None))
    info(f"Kullanılan Scorer'lar: {', '.join(requested_scorers)}")

    gene_names = [args.gene] if getattr(args, "gene", None) else None
    progress_bar = not getattr(args, "no_progress", False)
    workers = getattr(args, "workers", 10) or 10

    print(f"\n{DIM}{len(variants)} VCF varyantı AlphaGenome Atlas üzerinden sorgulanıyor...{RESET}")
    raw_scores = atlas.query_variants(
        variants,
        requested_scorers=requested_scorers,
        ontology_terms=ontology_terms,
        gene_names=gene_names,
        progress_bar=progress_bar,
        max_workers=workers,
    )
    ok("Atlas VCF sorgusu tamamlandı!")

    df = tidy_atlas_scores(raw_scores)
    if df is not None and not df.empty:
        ok(f"Toplam Atlas skoru: {df.shape[0]} satır × {df.shape[1]} sütun")
        sort_col = "quantile_score" if "quantile_score" in df.columns else "raw_score"
        print(f"\n{BOLD}📊 Özet:{RESET}")
        if "variant" in df.columns:
            print(f"  Benzersiz varyant sayısı: {df['variant'].nunique()}")
        if "gene_name" in df.columns:
            print(f"  Etkilenen gen sayısı: {df['gene_name'].nunique()}")
        print(f"  Max |skor|: {df[sort_col].abs().max():.6f}")

        top5 = df.reindex(df[sort_col].abs().sort_values(ascending=False).index).head(10)
        display_cols = [c for c in ["variant", "gene_name", "variant_scorer", "raw_score", "quantile_score"] if c in top5.columns]
        if display_cols:
            print(f"\n{BOLD}🏆 En yüksek skorlar (|{sort_col}|):{RESET}")
            print(top5[display_cols].to_string(index=False))

        output_file = getattr(args, "output", None) or (Path(args.vcf_file).stem + "_atlas_scores.tsv")
        save_dataframe(df, output_file)
    else:
        warn("VCF varyantları için Atlas skoru döndürülmedi.")

    ok("Atlas VCF analizi tamamlandı!")


def cmd_atlas_scorers(args):
    """AlphaGenome Atlas'ta mevcut tüm scorer'ları ve açıklamalarını listele."""
    header("AlphaGenome Atlas — Mevcut Scorer'lar")
    atlas = create_atlas(args)

    info("Atlas scorer metadatası alınıyor...")
    meta = atlas.scorer_metadata()
    if not meta:
        warn("Atlas scorer metadatası bulunamadı.")
        return

    search = getattr(args, "search", None)
    if search:
        search = search.lower()

    print(f"\n{BOLD}{'Scorer Adı':<35} {'İşaretli':<10} {'Track Sayısı':<15} {'Özet / Dokular'}{RESET}")
    print(f"{'─' * 80}")

    count = 0
    for name, s_meta in sorted(meta.items()):
        tm = s_meta.track_metadata
        n_tracks = len(tm) if (tm is not None and hasattr(tm, "__len__")) else 0
        is_signed_str = "Evet (+/-)" if s_meta.is_signed else "Hayır (pozitif)"

        tracks_desc = f"{n_tracks} track" if n_tracks > 0 else "Gen odaklı / ontolojisiz"

        # Arama filtresi
        if search:
            match_name = search in name.lower()
            match_tracks = False
            if tm is not None and not tm.empty:
                for col in ("name", "biosample_name", "ontology_curie"):
                    if col in tm.columns and tm[col].astype(str).str.lower().str.contains(search).any():
                        match_tracks = True
                        break
            if not (match_name or match_tracks):
                continue

        count += 1
        print(f"  {CYAN}{name:<33}{RESET} {is_signed_str:<10} {tracks_desc:<15}")

    print(f"\n{DIM}Toplam {count} scorer listelendi.{RESET}")
    print(f"{DIM}Kullanım: python alphagenome_cli.py atlas variant <variant> --scorers <scorer1> <scorer2>{RESET}")


def cmd_atlas(args):
    """Atlas üst komut yönlendiricisi."""
    if not getattr(args, "atlas_command", None):
        print(f"\n{YELLOW}Lütfen bir Atlas alt komutu seçin:{RESET}")
        print(f"  {CYAN}variant{RESET}   Tek bir varyantın Atlas skorlarını sorgula")
        print(f"  {CYAN}interval{RESET}  Genomik aralıktaki tüm varyantların skorlarını çek")
        print(f"  {CYAN}list{RESET}      TSV/CSV dosyasındaki varyantları Atlas'tan sorgula")
        print(f"  {CYAN}vcf{RESET}       VCF dosyasındaki varyantları Atlas'tan sorgula")
        print(f"  {CYAN}scorers{RESET}   Atlas'taki tüm scorer'ları listele")
        print(f"\nÖrnek: python alphagenome_cli.py atlas variant chr22:36201698:A>C")
        sys.exit(1)


# ─── Command: variant ──────────────────────────────────────────────────────
def cmd_variant(args):
    """Tek bir varyantın etkisini analiz et."""
    if getattr(args, "atlas", False):
        return cmd_atlas_variant(args)
    header("Tekli Varyant Analizi")
    model = create_model(args)

    variant = parse_variant_str(args.variant)
    info(f"Varyant: {variant.chromosome}:{variant.position} "
         f"{variant.reference_bases}>{variant.alternate_bases}")

    # Doku çözümleme
    ontology_terms = find_ontology_terms(model, args.tissue) if args.tissue else None

    # Modality çözümleme
    output_types = resolve_modalities(args.modality)
    info(f"Modality: {', '.join(args.modality)}")
    warn_missing_modality_terms(model, ontology_terms, args.modality)

    # Interval oluştur
    interval = auto_interval_for_variant(variant)
    info(f"Pencere: {interval.chromosome}:{interval.start}-{interval.end} "
         f"({interval.end - interval.start:,} bp)")

    # ---- Predict Variant ----
    print(f"\n{DIM}Tahmin yapılıyor...{RESET}")
    outputs = model.predict_variant(
        interval=interval,
        variant=variant,
        ontology_terms=ontology_terms,
        requested_outputs=output_types,
    )

    # Sonuçları göster
    print(f"\n{BOLD}📊 Sonuçlar:{RESET}")
    results_rows = []
    for mod_name in args.modality:
        attr = mod_name.lower()
        ref_data = getattr(outputs.reference, attr, None)
        alt_data = getattr(outputs.alternate, attr, None)
        if ref_data is not None and alt_data is not None:
            # Bu doku/ontoloji için track yoksa dizi boş olur → atla
            if ref_data.values.size == 0 or ref_data.values.shape[-1] == 0:
                warn(f"{mod_name}: bu doku için track yok (0 track) — atlanıyor")
                continue
            diff = alt_data.values - ref_data.values
            max_eff = np.abs(diff).max()
            mean_eff = np.abs(diff).mean()
            print(f"  {CYAN}{mod_name}:{RESET}")
            print(f"    REF aralığı : [{ref_data.values.min():.4f}, {ref_data.values.max():.4f}]")
            print(f"    ALT aralığı : [{alt_data.values.min():.4f}, {alt_data.values.max():.4f}]")
            print(f"    Max |etki|  : {max_eff:.6f}")
            print(f"    Ort |etki|  : {mean_eff:.6f}")
            results_rows.append({
                "variant": f"{variant.chromosome}:{variant.position}:{variant.reference_bases}>{variant.alternate_bases}",
                "modality": mod_name,
                "max_abs_effect": max_eff,
                "mean_abs_effect": mean_eff,
            })

    # ---- Score Variant ----
    df = None
    if args.score:
        print(f"\n{DIM}Varyant skorlanıyor...{RESET}")
        recommended = get_recommended_scorers(Organism.HOMO_SAPIENS)
        scores = model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=recommended,
        )
        df = tidy_scores(scores)
        if df is not None and not df.empty:
            ok(f"Skor tablosu: {df.shape[0]} satır × {df.shape[1]} sütun")
            # tidy_scores sütunları: gene_name, output_type, variant_scorer,
            # raw_score, quantile_score (varsa). Sıralama için quantile'ı tercih et.
            sort_col = "quantile_score" if "quantile_score" in df.columns else (
                "raw_score" if "raw_score" in df.columns else None)
            # --gene verildiyse tabloyu da o gene filtrele (yorumla tutarlı)
            top_df = df
            if getattr(args, "gene", None) and "gene_name" in df.columns:
                g = args.gene.upper()
                gmask = df["gene_name"].astype(str).str.upper() == g
                if not gmask.any():
                    gmask = df["gene_name"].astype(str).str.upper().str.contains(g, na=False)
                if gmask.any():
                    top_df = df[gmask]
            if sort_col:
                top = top_df.reindex(
                    top_df[sort_col].abs().sort_values(ascending=False).index
                ).head(10)
                gtag = f" · gen: {args.gene}" if getattr(args, "gene", None) else ""
                print(f"\n{BOLD}🏆 En yüksek skorlar (|{sort_col}|{gtag}):{RESET}")
                display_cols = [c for c in
                                ["gene_name", "output_type", "variant_scorer",
                                 "raw_score", "quantile_score"]
                                if c in top.columns]
                if display_cols:
                    print(top[display_cols].to_string(index=False))
            if args.output:
                save_dataframe(df, args.output)
        else:
            warn("Skor verisi döndürülmedi.")

    # Özet çıktıyı kaydet
    if args.output and not args.score:
        if results_rows:
            save_dataframe(pd.DataFrame(results_rows), args.output)

    # Otomatik yorum bloğu
    if not getattr(args, "no_interpret", False):
        interpret_results(variant, score_df=df, results_rows=results_rows,
                          gene=getattr(args, "gene", None))

    ok("Tekli varyant analizi tamamlandı!")


# ─── Command: vcf ──────────────────────────────────────────────────────────
def cmd_vcf(args):
    """VCF dosyasından varyantları analiz et."""
    if getattr(args, "atlas", False):
        return cmd_atlas_vcf(args)
    header("VCF Dosyasından Toplu Analiz")
    model = create_model(args)

    vcf_path = Path(args.vcf_file)
    variants = read_variants_from_vcf(vcf_path)

    # Doku çözümleme
    ontology_terms = find_ontology_terms(model, args.tissue) if args.tissue else None

    # Toplu skorlama
    info("Varyantlar skorlanıyor...")
    intervals = [auto_interval_for_variant(v) for v in variants]
    recommended = get_recommended_scorers(Organism.HOMO_SAPIENS)

    scores_list = model.score_variants(
        intervals=intervals,
        variants=variants,
        variant_scorers=recommended,
        progress_bar=True,
        max_workers=args.workers,
    )

    df = tidy_scores(scores_list)
    if df is not None and not df.empty:
        ok(f"Toplam skor: {df.shape[0]} satır × {df.shape[1]} sütun")
        print(f"\n{BOLD}📊 Özet:{RESET}")
        if "variant_id" in df.columns:
            df["variant_id"] = df["variant_id"].astype(str)
            print(f"  Benzersiz varyant sayısı: {df['variant_id'].nunique()}")
        if "gene_name" in df.columns:
            print(f"  Etkilenen gen sayısı: {df['gene_name'].nunique()}")
        if "score" in df.columns:
            print(f"  Max |skor|: {df['score'].abs().max():.6f}")

        if args.output:
            save_dataframe(df, args.output)
        else:
            # Sonuçları varsayılan dosyaya kaydet
            default_out = vcf_path.stem + "_scores.tsv"
            save_dataframe(df, default_out)
    else:
        warn("Skor verisi döndürülmedi.")

    ok("VCF analizi tamamlandı!")


# ─── Command: region ───────────────────────────────────────────────────────
def cmd_region(args):
    """Genomik bölge için ileri tahmin yap."""
    header("Bölge Bazlı İleri Tahmin")
    model = create_model(args)

    interval = parse_region_str(args.region)
    info(f"Bölge: {interval.chromosome}:{interval.start}-{interval.end}")
    info(f"Uzunluk: {interval.end - interval.start:,} bp")

    # Doku çözümleme
    ontology_terms = find_ontology_terms(model, args.tissue) if args.tissue else None

    # Modality çözümleme
    output_types = resolve_modalities(args.modality)
    info(f"Modality: {', '.join(args.modality)}")
    warn_missing_modality_terms(model, ontology_terms, args.modality)

    print(f"\n{DIM}Tahmin yapılıyor...{RESET}")
    outputs = model.predict_interval(
        interval=interval,
        ontology_terms=ontology_terms,
        requested_outputs=output_types,
    )

    print(f"\n{BOLD}📊 Sonuçlar:{RESET}")
    results = {}
    for mod_name in args.modality:
        attr = mod_name.lower()
        data = getattr(outputs, attr, None)
        if data is not None:
            if data.values.size == 0 or data.values.shape[-1] == 0:
                warn(f"{mod_name}: bu doku için track yok (0 track) — atlanıyor")
                continue
            print(f"  {CYAN}{mod_name}:{RESET}")
            print(f"    Çıktı aralığı: {data.interval}")
            print(f"    Veri boyutu  : {data.values.shape}")
            print(f"    Değer aralığı: [{data.values.min():.4f}, {data.values.max():.4f}]")
            print(f"    Ortalama     : {data.values.mean():.4f}")
            results[mod_name] = data

    if args.output and results:
        # Her modality için ayrı dosya ya da birleşik tablo kaydet
        summary = []
        for mod_name, data in results.items():
            summary.append({
                "region": f"{interval.chromosome}:{interval.start}-{interval.end}",
                "modality": mod_name,
                "shape": str(data.values.shape),
                "min": data.values.min(),
                "max": data.values.max(),
                "mean": data.values.mean(),
            })
        save_dataframe(pd.DataFrame(summary), args.output)

    ok("Bölge analizi tamamlandı!")


# ─── Command: list ──────────────────────────────────────────────────────────
def cmd_list(args):
    """TSV/CSV dosyasından varyantları analiz et."""
    if getattr(args, "atlas", False):
        return cmd_atlas_list(args)
    header("Dosyadan Toplu Varyant Analizi")
    model = create_model(args)

    list_path = Path(args.list_file)
    if not list_path.exists():
        err(f"Dosya bulunamadı: {list_path}")
        sys.exit(1)

    # Doku çözümleme
    ontology_terms = find_ontology_terms(model, args.tissue) if args.tissue else None

    # Dosyayı oku
    info(f"Dosya okunuyor: {list_path}")
    sep = "," if list_path.suffix == ".csv" else "\t"

    variants = []
    with open(list_path, "r") as f:
        reader = csv.reader(f, delimiter=sep)
        header_row = next(reader, None)

        # Header kontrolü — eğer ilk satır veri gibiyse, onu da işle
        has_header = False
        if header_row:
            first_val = header_row[0].lower().strip()
            if first_val in ("chr", "chrom", "chromosome", "#chrom", "#chr"):
                has_header = True
            elif re.match(r"^chr", first_val, re.IGNORECASE):
                # İlk satır doğrudan veri
                chrom = header_row[0] if header_row[0].startswith("chr") else f"chr{header_row[0]}"
                variants.append(
                    genome.Variant(
                        chromosome=chrom,
                        position=int(header_row[1]),
                        reference_bases=header_row[2].upper(),
                        alternate_bases=header_row[3].upper(),
                    )
                )

        for row in reader:
            if len(row) < 4:
                continue
            if row[0].startswith("#"):
                continue
            chrom = row[0] if row[0].startswith("chr") else f"chr{row[0]}"
            try:
                pos = int(row[1])
            except ValueError:
                continue
            ref = row[2].upper()
            alt = row[3].upper()
            variants.append(
                genome.Variant(
                    chromosome=chrom,
                    position=pos,
                    reference_bases=ref,
                    alternate_bases=alt,
                )
            )

    if not variants:
        err("Dosyada varyant bulunamadı.")
        info("Beklenen format (TSV): chr\\tpos\\tref\\talt")
        sys.exit(1)

    ok(f"{len(variants)} varyant okundu")

    # Toplu skorlama
    info("Varyantlar skorlanıyor...")
    intervals = [auto_interval_for_variant(v) for v in variants]
    recommended = get_recommended_scorers(Organism.HOMO_SAPIENS)

    scores_list = model.score_variants(
        intervals=intervals,
        variants=variants,
        variant_scorers=recommended,
        progress_bar=True,
        max_workers=args.workers,
    )

    df = tidy_scores(scores_list)
    if df is not None and not df.empty:
        ok(f"Toplam skor: {df.shape[0]} satır × {df.shape[1]} sütun")
        print(f"\n{BOLD}📊 Özet:{RESET}")
        if "variant_id" in df.columns:
            df["variant_id"] = df["variant_id"].astype(str)
            print(f"  Benzersiz varyant sayısı: {df['variant_id'].nunique()}")
        if "gene_name" in df.columns:
            print(f"  Etkilenen gen sayısı: {df['gene_name'].nunique()}")
        if "score" in df.columns:
            top5 = df.reindex(df["score"].abs().sort_values(ascending=False).index).head(5)
            display_cols = [c for c in ["variant_id", "gene_name", "scorer_name", "score"] if c in top5.columns]
            if display_cols:
                print(f"\n{BOLD}🏆 En yüksek 5 skor:{RESET}")
                print(top5[display_cols].to_string(index=False))

        output_file = args.output or (list_path.stem + "_scores.tsv")
        save_dataframe(df, output_file)
    else:
        warn("Skor verisi döndürülmedi.")

    ok("Toplu analiz tamamlandı!")


# ─── Command: tissues ───────────────────────────────────────────────────────
def cmd_tissues(args):
    """Mevcut doku/ontoloji terimlerini listele ve ara."""
    header("Doku / Ontoloji Terimi Arama")
    model = create_model(args)

    meta = model.output_metadata()
    rna_meta = meta.rna_seq

    # Benzersiz dokuları çıkar
    tissues = rna_meta[["ontology_curie", "biosample_name", "biosample_type", "biosample_life_stage"]].drop_duplicates()
    tissues = tissues.sort_values("biosample_name")

    if args.search:
        query = args.search.lower()
        mask = (
            tissues["biosample_name"].str.lower().str.contains(query, na=False)
            | tissues["ontology_curie"].str.lower().str.contains(query, na=False)
            | tissues["biosample_type"].str.lower().str.contains(query, na=False)
        )
        # gtex_tissue'da da ara
        if "gtex_tissue" in rna_meta.columns:
            gtex_matches = rna_meta[
                rna_meta["gtex_tissue"].str.lower().str.contains(query, na=False)
            ]["ontology_curie"].unique()
            mask = mask | tissues["ontology_curie"].isin(gtex_matches)

        filtered = tissues[mask]

        if filtered.empty:
            warn(f"'{args.search}' ile eşleşen doku bulunamadı.")
            info("Tüm dokuları listelemek için: python alphagenome_cli.py tissues")
            return

        print(f"\n{BOLD}'{args.search}' araması — {len(filtered)} sonuç:{RESET}\n")
        for _, row in filtered.iterrows():
            print(f"  {CYAN}{row['ontology_curie']:<20}{RESET} "
                  f"{row['biosample_name']:<35} "
                  f"{DIM}({row['biosample_type']}, {row['biosample_life_stage']}){RESET}")
    else:
        print(f"\n{BOLD}Toplam {len(tissues)} benzersiz doku/hücre tipi:{RESET}\n")
        for _, row in tissues.iterrows():
            print(f"  {CYAN}{row['ontology_curie']:<20}{RESET} "
                  f"{row['biosample_name']:<35} "
                  f"{DIM}({row['biosample_type']}){RESET}")

    print(f"\n{DIM}Kullanım: --tissue <doku_adı> veya --tissue UBERON:0001157{RESET}")


# ─── Command: visualize ─────────────────────────────────────────────────────
def cmd_visualize(args):
    """Sonuçları görselleştir."""
    header("AlphaGenome Görselleştirme")
    model = create_model(args)

    # Doku çözümleme
    ontology_terms = find_ontology_terms(model, args.tissue) if args.tissue else None

    # Modality çözümleme
    output_types = resolve_modalities(args.modality)
    warn_missing_modality_terms(model, ontology_terms, args.modality)

    def _safe_ylabel(mod_name, tdata_obj):
        """Metadata sütunlarına göre güvenli ylabel şablonu oluşturur."""
        try:
            cols = set(tdata_obj.metadata.columns)
        except Exception:
            return mod_name
        parts = [mod_name + ":"]
        if "biosample_name" in cols:
            parts.append("{biosample_name}")
        elif "name" in cols:
            parts.append("{name}")
        if "strand" in cols:
            parts.append(":{strand}")
        return " ".join(parts) if len(parts) > 1 else mod_name

    # ── Variant mode ──
    if args.mode == "variant":
        variant = parse_variant_str(args.target)
        interval = auto_interval_for_variant(variant)
        info(f"Varyant: {variant.chromosome}:{variant.position} "
             f"{variant.reference_bases}>{variant.alternate_bases}")
        info(f"Pencere: {interval.chromosome}:{interval.start}-{interval.end}")
        info(f"Modality: {', '.join(args.modality)}")

        # Gen anotasyonları
        transcripts = load_transcripts(interval, transcript_set=args.transcripts)

        # Tahmin
        print(f"\n{DIM}Tahmin yapılıyor (REF ve ALT)...{RESET}")
        output = model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=ontology_terms,
            requested_outputs=output_types,
        )
        ok("Tahmin tamamlandı!")

        # Plot oluştur
        components = []
        if transcripts:
            components.append(plot_components.TranscriptAnnotation(transcripts))

        for mod_name in args.modality:
            attr = mod_name.lower()
            ref_data = getattr(output.reference, attr, None)
            alt_data = getattr(output.alternate, attr, None)
            if ref_data is not None and alt_data is not None:
                # Check if data actually has tracks (columns > 0)
                n_tracks = ref_data.values.shape[1] if hasattr(ref_data, 'values') and ref_data.values.ndim >= 2 else 0
                if n_tracks == 0:
                    warn(f"{mod_name}: Bu doku için veri yok (0 track). "
                         f"Farklı bir doku/ontoloji terimi deneyin.")
                    continue
                ylabel = _safe_ylabel(mod_name, ref_data)
                components.append(
                    plot_components.OverlaidTracks(
                        tdata={"REF": ref_data, "ALT": alt_data},
                        colors=REF_ALT_COLORS,
                        ylabel_template=ylabel,
                    )
                )
            else:
                warn(f"{mod_name}: Bu modalite için tahmin verisi döndürülmedi.")

        if not components:
            err("Görselleştirilecek veri bulunamadı.")
            return

        # Zoom penceresi — varyant ±16kb
        zoom = args.zoom or 2**15  # 32768 bp default
        plot_interval = genome.Interval(
            chromosome=variant.chromosome,
            start=max(interval.start, variant.position - zoom // 2),
            end=min(interval.end, variant.position + zoom // 2),
        )

        title = (f"Variant Effect: {variant.chromosome}:{variant.position} "
                 f"{variant.reference_bases}>{variant.alternate_bases}")

        fig = plot_components.plot(
            components,
            interval=plot_interval,
            title=title,
            annotations=[
                plot_components.VariantAnnotation([variant], alpha=0.8)
            ],
            fig_width=20,
        )

        out_path = args.output or f"variant_{variant.chromosome}_{variant.position}.png"
        save_figure(fig, out_path)

    # ── Region mode ──
    elif args.mode == "region":
        interval = parse_region_str(args.target)
        info(f"Bölge: {interval.chromosome}:{interval.start}-{interval.end}")
        info(f"Modality: {', '.join(args.modality)}")

        # Gen anotasyonları
        transcripts = load_transcripts(interval, transcript_set=args.transcripts)

        # Tahmin
        print(f"\n{DIM}Tahmin yapılıyor...{RESET}")
        output = model.predict_interval(
            interval=interval,
            ontology_terms=ontology_terms,
            requested_outputs=output_types,
        )
        ok("Tahmin tamamlandı!")

        # Plot oluştur
        components = []
        if transcripts:
            components.append(plot_components.TranscriptAnnotation(transcripts))

        for mod_name in args.modality:
            attr = mod_name.lower()
            data = getattr(output, attr, None)
            if data is not None:
                # Check if data actually has tracks (columns > 0)
                n_tracks = data.values.shape[1] if hasattr(data, 'values') and data.values.ndim >= 2 else 0
                if n_tracks == 0:
                    warn(f"{mod_name}: Bu doku için veri yok (0 track). "
                         f"Farklı bir doku/ontoloji terimi deneyin.")
                    continue
                ylabel = _safe_ylabel(mod_name, data)
                components.append(
                    plot_components.Tracks(
                        tdata=data,
                        ylabel_template=ylabel,
                    )
                )
            else:
                warn(f"{mod_name}: Bu modalite için tahmin verisi döndürülmedi.")

        if not components:
            err("Görselleştirilecek veri bulunamadı.")
            return

        # Zoom penceresi (opsiyonel)
        plot_interval = interval
        if args.zoom:
            center = (interval.start + interval.end) // 2
            plot_interval = genome.Interval(
                chromosome=interval.chromosome,
                start=max(interval.start, center - args.zoom // 2),
                end=min(interval.end, center + args.zoom // 2),
            )

        title = f"Region: {interval.chromosome}:{interval.start}-{interval.end}"

        fig = plot_components.plot(
            components,
            interval=plot_interval,
            title=title,
            fig_width=20,
        )

        out_path = args.output or f"region_{interval.chromosome}_{interval.start}_{interval.end}.png"
        save_figure(fig, out_path)

    ok("Görselleştirme tamamlandı!")


# ─── Command: ui ───────────────────────────────────────────────────────────
def cmd_ui(args):
    """Modern Web Arayüzünü Başlat."""
    header("AlphaGenome Web UI — Modern Görsel Analiz Platformu")
    try:
        from web_app import start_web_server
    except ImportError:
        err("Web UI modülü (web_app.py) bulunamadı.")
        sys.exit(1)

    host = getattr(args, "host", "127.0.0.1")
    port = getattr(args, "port", 8000)
    open_browser = not getattr(args, "no_browser", False)

    start_web_server(host=host, port=port, open_browser=open_browser)


# ─── Argument Parsing ───────────────────────────────────────────────────────
def build_parser():
    parser = argparse.ArgumentParser(
        description="AlphaGenome CLI — Kolay genomik analiz aracı",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Örnekler:
  %(prog)s variant chr22:36201698:A>C --tissue colon
  %(prog)s variant chr22:36201698:A>C --atlas
  %(prog)s atlas variant chr22:36201698:A>C --tissue colon
  %(prog)s atlas interval chr22:35677410-35678410 --tissue liver
  %(prog)s atlas list variants.tsv --output atlas_scores.tsv
  %(prog)s atlas scorers
  %(prog)s ui --port 8000
  %(prog)s region chr22:35677410-36725986 --tissue liver --modality RNA_SEQ ATAC
  %(prog)s vcf variants.vcf --tissue brain --output results.tsv
  %(prog)s list variants.tsv --tissue colon
  %(prog)s tissues --search liver
  %(prog)s visualize variant chr22:36201698:A>C --tissue colon
  %(prog)s visualize region chr22:35677410-36725986 --tissue colon --modality RNA_SEQ ATAC
        """,
    )
    parser.add_argument(
        "--api-key",
        help="AlphaGenome API anahtarı (yoksa ALPHAGENOME_API_KEY env değişkeni kullanılır)",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Önbelleklemeyi devre dışı bırak — her zaman API'ye gönder",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help=f"Önbellek dizini (varsayılan: {DEFAULT_CACHE_DIR})",
    )

    subparsers = parser.add_subparsers(dest="command", help="Alt komut seçin")

    # ── variant ──
    p_var = subparsers.add_parser(
        "variant",
        help="Tek varyant analizi",
        description="Tek bir genetik varyantın etkisini analiz eder.",
    )
    p_var.add_argument("variant", help="Varyant (chr22:36201698:A>C formatında)")
    p_var.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_var.add_argument(
        "--modality", "-m", nargs="+", default=["RNA_SEQ"],
        help="Çıktı tipleri (varsayılan: RNA_SEQ). Seçenekler: RNA_SEQ, CAGE, ATAC, DNASE, PROCAP, SPLICE_SITES..."
    )
    p_var.add_argument("--score", "-s", action="store_true", help="Varyantı skorla (detaylı etki analizi)")
    p_var.add_argument("--atlas", action="store_true", help="Model tahmini yerine AlphaGenome Atlas'tan sorgula")
    p_var.add_argument("--scorers", nargs="+", help="Atlas için scorer filtreleri")
    p_var.add_argument("--output", "-o", help="Sonuçları dosyaya kaydet (TSV/CSV)")
    p_var.add_argument("--no-interpret", action="store_true",
                       help="Otomatik yorum/özet bloğunu gösterme")
    p_var.add_argument("--gene", "-g",
                       help="Yorumu belirli bir gene odakla (ör. FAS); "
                            "1 Mb penceredeki diğer genleri dışlar")
    p_var.set_defaults(func=cmd_variant)

    # ── vcf ──
    p_vcf = subparsers.add_parser(
        "vcf",
        help="VCF dosyasından toplu analiz",
        description="VCF dosyasındaki tüm varyantları toplu olarak skorlar.",
    )
    p_vcf.add_argument("vcf_file", help="VCF dosyası (.vcf veya .vcf.gz)")
    p_vcf.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_vcf.add_argument("--atlas", action="store_true", help="Model inferansı yerine AlphaGenome Atlas'tan sorgula")
    p_vcf.add_argument("--scorers", nargs="+", help="Atlas için scorer filtreleri")
    p_vcf.add_argument("--output", "-o", help="Çıktı dosyası (varsayılan: <vcf>_scores.tsv)")
    p_vcf.add_argument("--workers", "-w", type=int, default=5, help="Paralel çalışan sayısı (varsayılan: 5)")
    p_vcf.set_defaults(func=cmd_vcf)

    # ── region ──
    p_reg = subparsers.add_parser(
        "region",
        help="Genomik bölge analizi",
        description="Belirtilen genomik bölge için ileri tahmin yapar.",
    )
    p_reg.add_argument("region", help="Bölge (chr22:35677410-36725986 formatında)")
    p_reg.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_reg.add_argument(
        "--modality", "-m", nargs="+", default=["RNA_SEQ"],
        help="Çıktı tipleri (varsayılan: RNA_SEQ)"
    )
    p_reg.add_argument("--output", "-o", help="Sonuçları dosyaya kaydet")
    p_reg.set_defaults(func=cmd_region)

    # ── list ──
    p_lst = subparsers.add_parser(
        "list",
        help="TSV/CSV dosyasından toplu analiz",
        description="TSV veya CSV dosyasındaki varyantları toplu olarak skorlar.\n"
                    "Dosya formatı: chr  pos  ref  alt (tab veya virgülle ayrılmış)",
    )
    p_lst.add_argument("list_file", help="Varyant listesi dosyası (.tsv veya .csv)")
    p_lst.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_lst.add_argument("--atlas", action="store_true", help="Model inferansı yerine AlphaGenome Atlas'tan sorgula")
    p_lst.add_argument("--scorers", nargs="+", help="Atlas için scorer filtreleri")
    p_lst.add_argument("--output", "-o", help="Çıktı dosyası (varsayılan: <dosya>_scores.tsv)")
    p_lst.add_argument("--workers", "-w", type=int, default=5, help="Paralel çalışan sayısı (varsayılan: 5)")
    p_lst.set_defaults(func=cmd_list)

    # ── tissues ──
    p_tis = subparsers.add_parser(
        "tissues",
        help="Mevcut doku/ontoloji terimlerini listele",
        description="AlphaGenome'da mevcut olan doku ve hücre tiplerini listeler.",
    )
    p_tis.add_argument("--search", "-s", help="Doku adında arama yapılacak kelime")
    p_tis.set_defaults(func=cmd_tissues)

    # ── visualize ──
    p_vis = subparsers.add_parser(
        "visualize",
        help="Sonuçları görselleştir (PNG/PDF/SVG)",
        description="Genomik tahminleri grafik olarak görselleştirir.\n"
                    "Gen anotasyonları, çoklu modality ve varyant etkilerini içerir.",
    )
    p_vis.add_argument("mode", choices=["variant", "region"], help="Görselleştirme modu")
    p_vis.add_argument("target", help="Varyant (chr:pos:ref>alt) veya bölge (chr:start-end)")
    p_vis.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_vis.add_argument(
        "--modality", "-m", nargs="+", default=["RNA_SEQ"],
        help="Çıktı tipleri (varsayılan: RNA_SEQ)"
    )
    p_vis.add_argument("--output", "-o", help="Çıktı dosyası (varsayılan: otomatik .png)")
    p_vis.add_argument(
        "--zoom", "-z", type=int, default=None,
        help="Zoom penceresi boyutu (bp). Variant modu varsayılanı: 32768 bp"
    )
    p_vis.add_argument(
        "--transcripts", choices=["mane", "longest"], default="mane",
        help="Gösterilecek transkript: mane=MANE Select kanonik (varsayılan), "
             "longest=gen başına en uzun izoform"
    )
    p_vis.set_defaults(func=cmd_visualize)

    # ── atlas ──
    p_atlas = subparsers.add_parser(
        "atlas",
        help="AlphaGenome Atlas (9 milyar önceden hesaplanmış SNV veritabanı)",
        description="AlphaGenome Atlas — İnsan genomundaki 9 milyar olası tek nükleotid varyantının\n"
                    "önceden hesaplanmış düzenleyici etki tahminlerini (1 petabayt veri, AVI skorları)\n"
                    "model çıkarımı beklemeden anında sorgulayın.",
    )
    p_atlas.set_defaults(func=cmd_atlas)
    atlas_subparsers = p_atlas.add_subparsers(dest="atlas_command", help="Atlas alt komutu seçin")

    # atlas variant
    p_at_var = atlas_subparsers.add_parser(
        "variant",
        help="Tek varyantın Atlas skorlarını sorgula",
        description="Tek bir varyantın önceden hesaplanmış Atlas skorlarını sorgular.",
    )
    p_at_var.add_argument("variant", help="Varyant (chr22:36201698:A>C formatında)")
    p_at_var.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi (ör. colon veya UBERON:0001157)")
    p_at_var.add_argument("--gene", "-g", help="Belirli bir gene odakla (ör. APOL4)")
    p_at_var.add_argument("--scorers", "-s", nargs="+", help="İstenen scorer'lar")
    p_at_var.add_argument("--output", "-o", help="Sonuçları dosyaya kaydet (TSV/CSV)")
    p_at_var.add_argument("--no-interpret", action="store_true", help="Otomatik yorum bloğunu gösterme")
    p_at_var.set_defaults(func=cmd_atlas_variant)

    # atlas interval
    p_at_int = atlas_subparsers.add_parser(
        "interval",
        help="Genomik aralıktaki tüm varyantların Atlas skorlarını sorgula",
        description="Genomik aralıktaki tüm olası tek nükleotid varyantlarının Atlas skorlarını çeker.",
    )
    p_at_int.add_argument("interval", help="Genomik aralık (chr22:35677410-35678410 formatında)")
    p_at_int.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_at_int.add_argument("--gene", "-g", help="Gen filtresi")
    p_at_int.add_argument("--scorers", "-s", nargs="+", help="İstenen scorer'lar")
    p_at_int.add_argument("--output", "-o", help="Çıktı dosyası (varsayılan: atlas_<aralık>.tsv)")
    p_at_int.add_argument("--workers", "-w", type=int, default=10, help="Paralel çalışan sayısı (varsayılan: 10)")
    p_at_int.add_argument("--no-progress", action="store_true", help="İlerleme çubuğunu gizle")
    p_at_int.set_defaults(func=cmd_atlas_interval)

    # atlas list
    p_at_lst = atlas_subparsers.add_parser(
        "list",
        help="TSV/CSV dosyasındaki varyantları Atlas'tan toplu sorgula",
        description="Dosyadaki varyantları paralel olarak AlphaGenome Atlas'tan sorgular.",
    )
    p_at_lst.add_argument("list_file", help="Varyant listesi dosyası (.tsv veya .csv)")
    p_at_lst.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_at_lst.add_argument("--gene", "-g", help="Gen filtresi")
    p_at_lst.add_argument("--scorers", "-s", nargs="+", help="İstenen scorer'lar")
    p_at_lst.add_argument("--output", "-o", help="Çıktı dosyası (varsayılan: <dosya>_atlas_scores.tsv)")
    p_at_lst.add_argument("--workers", "-w", type=int, default=10, help="Paralel çalışan sayısı (varsayılan: 10)")
    p_at_lst.add_argument("--no-progress", action="store_true", help="İlerleme çubuğunu gizle")
    p_at_lst.set_defaults(func=cmd_atlas_list)

    # atlas vcf
    p_at_vcf = atlas_subparsers.add_parser(
        "vcf",
        help="VCF dosyasındaki varyantları Atlas'tan toplu sorgula",
        description="VCF dosyasındaki varyantları paralel olarak AlphaGenome Atlas'tan sorgular.",
    )
    p_at_vcf.add_argument("vcf_file", help="VCF dosyası (.vcf veya .vcf.gz)")
    p_at_vcf.add_argument("--tissue", "-t", help="Doku adı veya ontoloji terimi")
    p_at_vcf.add_argument("--gene", "-g", help="Gen filtresi")
    p_at_vcf.add_argument("--scorers", "-s", nargs="+", help="İstenen scorer'lar")
    p_at_vcf.add_argument("--output", "-o", help="Çıktı dosyası (varsayılan: <vcf>_atlas_scores.tsv)")
    p_at_vcf.add_argument("--workers", "-w", type=int, default=10, help="Paralel çalışan sayısı (varsayılan: 10)")
    p_at_vcf.add_argument("--no-progress", action="store_true", help="İlerleme çubuğunu gizle")
    p_at_vcf.set_defaults(func=cmd_atlas_vcf)

    # atlas scorers
    p_at_scr = atlas_subparsers.add_parser(
        "scorers",
        help="Mevcut tüm Atlas scorer'larını ve açıklamalarını listele",
        description="Atlas veritabanındaki tüm önceden hesaplanmış scorer'ları, yön işaretlerini ve track sayılarını listeler.",
    )
    p_at_scr.add_argument("--search", "-s", help="Scorer veya doku adında arama terimi")
    p_at_scr.set_defaults(func=cmd_atlas_scorers)

    # ── ui ──
    p_ui = subparsers.add_parser(
        "ui",
        help="Modern Web Arayüzünü Başlat (Drag & Drop, Atlas, Görselleştirme)",
        description="Sürükle-bırak dosya yükleme, varyant arama ve görsel analiz destekli Web UI başlatır.",
    )
    p_ui.add_argument("--host", default="127.0.0.1", help="Sunucu host adresi (varsayılan: 127.0.0.1)")
    p_ui.add_argument("--port", type=int, default=8000, help="Sunucu portu (varsayılan: 8000)")
    p_ui.add_argument("--no-browser", action="store_true", help="Web tarayıcısını otomatik olarak açma")
    p_ui.set_defaults(func=cmd_ui)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        print(f"\n{YELLOW}Bir alt komut seçin: variant, vcf, region, list, tissues, visualize, atlas, ui{RESET}")
        sys.exit(1)

    try:
        args.func(args)
    except KeyboardInterrupt:
        print(f"\n{YELLOW}İptal edildi.{RESET}")
        sys.exit(130)
    except Exception as e:
        err(f"Hata: {e}")
        raise


if __name__ == "__main__":
    main()
