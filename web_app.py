"""
AlphaGenome Studio Web Server
=============================
Modern Web UI backend for AlphaGenome Atlas and Genomic Model inference.
Zero external dependencies (uses Python standard library http.server).

Usage:
    python alphagenome_cli.py ui --port 8000
    or
    python web_app.py
"""

import io
import json
import math
import mimetypes
import os
import re
import sys
import tempfile
import threading
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

# AlphaGenome CLI and backend integrations
import alphagenome_cli as cli
from alphagenome.data import genome

STATIC_DIR = Path(__file__).parent / "ui"


def sanitize_val(val):
    """JSON uyumlu değer dönüşümü (NaN ve Inf değerlerini None yapar)."""
    if val is None:
        return None
    if isinstance(val, (float, int)):
        if math.isnan(val) or math.isinf(val):
            return None
        return val
    return str(val)


def df_to_json_records(df):
    """Pandas DataFrame'i güvenli JSON listesine dönüştürür."""
    if df is None or df.empty:
        return []
    records = []
    for _, row in df.iterrows():
        rec = {}
        for col in df.columns:
            val = row[col]
            if hasattr(val, "item"):
                try:
                    val = val.item()
                except Exception:
                    pass
            rec[col] = sanitize_val(val)
        records.append(rec)
    return records


KEY_FILE = Path.home() / ".alphagenome_key"

def get_stored_api_key():
    k = os.environ.get("ALPHAGENOME_API_KEY")
    if k and len(k.strip()) > 5:
        return k.strip()
    if KEY_FILE.exists():
        try:
            stored = KEY_FILE.read_text().strip()
            if stored and len(stored) > 5:
                os.environ["ALPHAGENOME_API_KEY"] = stored
                return stored
        except Exception:
            pass
    return None


class AlphaGenomeWebHandler(BaseHTTPRequestHandler):

    def send_json(self, data, status=200):
        body = json.dumps(data, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def send_error_json(self, message, status=400, **kwargs):
        payload = {"error": str(message)}
        payload.update(kwargs)
        self.send_json(payload, status=status)

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path == "/api/status":
            self.handle_api_status()
            return
        elif path == "/api/scorers":
            self.handle_api_scorers()
            return

        if path == "/favicon.ico":
            self.send_response(204)
            self.end_headers()
            return

        # Static file serving
        if path == "/" or path == "":
            file_path = STATIC_DIR / "index.html"
        else:
            clean_path = path.lstrip("/")
            file_path = STATIC_DIR / clean_path

        if file_path.exists() and file_path.is_file():
            mime, _ = mimetypes.guess_type(str(file_path))
            if mime is None:
                mime = "application/octet-stream"
            try:
                with open(file_path, "rb") as f:
                    content = f.read()
                self.send_response(200)
                self.send_header("Content-Type", f"{mime}; charset=utf-8" if "text" in mime or "javascript" in mime else mime)
                self.send_header("Content-Length", str(len(content)))
                self.end_headers()
                self.wfile.write(content)
                return
            except Exception:
                self.send_error(500, "Internal Server Error")
                return

        self.send_error(404, "Not Found")

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        content_len = int(self.headers.get("Content-Length", 0))
        post_body = self.rfile.read(content_len).decode("utf-8")

        try:
            req_data = json.loads(post_body) if post_body else {}
        except json.JSONDecodeError:
            self.send_error_json("Geçersiz JSON verisi")
            return

        if path == "/api/set_key":
            self.handle_api_set_key(req_data)
        elif path == "/api/query_variant":
            self.handle_api_query_variant(req_data)
        elif path == "/api/query_interval":
            self.handle_api_query_interval(req_data)
        elif path == "/api/query_batch":
            self.handle_api_query_batch(req_data)
        else:
            self.send_error(404, "Not Found")

    def handle_api_status(self):
        api_key = get_stored_api_key()
        self.send_json({
            "has_key": bool(api_key and len(api_key.strip()) > 5),
            "cache_dir": str(cli.DEFAULT_CACHE_DIR),
            "version": "0.9.0",
        })

    def handle_api_set_key(self, data):
        key = data.get("api_key", "").strip()
        if not key:
            self.send_error_json("API anahtarı boş olamaz")
            return
        os.environ["ALPHAGENOME_API_KEY"] = key
        try:
            KEY_FILE.write_text(key)
        except Exception:
            pass
        self.send_json({"success": True, "message": "API anahtarı güncellendi"})

    def handle_api_scorers(self):
        # Scorer metadata
        api_key = os.environ.get("ALPHAGENOME_API_KEY")
        if api_key:
            try:
                class DummyArgs:
                    api_key = None
                    no_cache = False
                    cache_dir = None
                atlas_client = cli.create_atlas(DummyArgs())
                metadata = atlas_client.scorer_metadata()
                scorers_list = []
                for name, meta in metadata.items():
                    scorers_list.append({
                        "name": name,
                        "description": getattr(meta, "description", ""),
                        "is_signed": bool(getattr(meta, "is_signed", True)),
                        "track_count": len(getattr(meta, "tracks", [])) if hasattr(meta, "tracks") else 0,
                    })
                self.send_json({"scorers": scorers_list})
                return
            except Exception:
                pass

        # Fallback varsayılan liste
        fallback_scorers = [
            {"name": "GeneMaskLFCScorer", "description": "Gen maskeleme log katlama değişimi (LFC) düzenleyici skor.", "is_signed": True, "track_count": 25},
            {"name": "GeneMaskActiveScorer", "description": "Gen aktivasyon ve ekspresyon seviyesi etki skoru.", "is_signed": True, "track_count": 25},
            {"name": "SpliceJunctionScorer", "description": "Splicing kavşakları ve ekzon tanıma etki skoru.", "is_signed": True, "track_count": 18},
            {"name": "CenterMaskScorer", "description": "Kromatin erişilebilirliği ve promotor bölgesi merkez maskesi.", "is_signed": True, "track_count": 30},
        ]
        self.send_json({"scorers": fallback_scorers})

    def handle_api_query_variant(self, data):
        var_str = data.get("variant", "").strip()
        engine = data.get("engine", "atlas")
        tissue = data.get("tissue", "").strip() or None
        gene = data.get("gene", "").strip() or None

        if not var_str:
            self.send_error_json("Varyant belirtilmedi")
            return

        try:
            var_obj = cli.parse_variant(var_str)
        except Exception as e:
            self.send_error_json(f"Geçersiz varyant formatı ({var_str}): {e}")
            return

        # 1. Önbellekten hızlı servis kontrolü (Atlas)
        if engine == "atlas":
            try:
                cached_client = cli.CachedAtlasClient(None, enabled=True)
                terms = cli.resolve_tissue_term(tissue) if tissue else None
                genes = [gene] if gene else None
                key = cli._cache_key(
                    "variant",
                    var_obj.chromosome, var_obj.position, var_obj.reference_bases, var_obj.alternate_bases,
                    tuple(),
                    tuple(sorted(str(o) for o in (terms or []))),
                    tuple(),
                    tuple(sorted(genes or [])),
                )
                hit, cached_scores = cached_client._load(key)
                if not hit:
                    import glob, pickle
                    for cf in glob.glob(str(cli.DEFAULT_CACHE_DIR / "atlas_*.pkl")):
                        try:
                            with open(cf, "rb") as f:
                                d = pickle.load(f)
                            if isinstance(d, dict) and "AVI_SCORE" in d and d["AVI_SCORE"] is not None:
                                v_str = str(d["AVI_SCORE"].obs.iloc[0]["variant"])
                                if v_str == var_str:
                                    cached_scores = d
                                    hit = True
                                    break
                        except Exception:
                            continue

                if hit and cached_scores:
                    df = cli.tidy_atlas_scores(cached_scores)
                    records = df_to_json_records(df)
                    avi_summary = cli.extract_atlas_summary(cached_scores)
                    self.send_json({
                        "scores": records,
                        "variant": var_str,
                        "engine": "atlas",
                        "avi_summary": avi_summary,
                    })
                    return
            except Exception:
                pass

        api_key = get_stored_api_key()
        if not api_key:
            self.send_error_json("AlphaGenome API anahtarı bulunamadı. Lütfen sağ üstteki ayarlar simgesinden API anahtarınızı girin.")
            return

        class DummyArgs:
            api_key = None
            no_cache = False
            cache_dir = None

        try:
            if engine == "atlas":
                atlas_client = cli.create_atlas(DummyArgs())
                terms = cli.resolve_tissue_term(tissue) if tissue else None
                genes = [gene] if gene else None

                raw_scores = atlas_client.query_variant(
                    var_obj,
                    ontology_terms=terms,
                    gene_names=genes,
                )
                df = cli.tidy_atlas_scores(raw_scores)
                records = df_to_json_records(df)
                avi_summary = cli.extract_atlas_summary(raw_scores)
                self.send_json({
                    "scores": records,
                    "variant": var_str,
                    "engine": "atlas",
                    "avi_summary": avi_summary,
                })
            else:
                # Model inference
                model = cli.create_model(DummyArgs())
                interval = cli.get_variant_interval(var_obj)
                terms = cli.resolve_tissue_term(tissue) if tissue else None
                from alphagenome.models.variant_scorers import get_recommended_scorers, tidy_scores
                from alphagenome.models.dna_model import Organism

                scorers = get_recommended_scorers(Organism.HOMO_SAPIENS)
                scores = model.score_variant(
                    interval=interval,
                    variant=var_obj,
                    variant_scorers=scorers,
                )
                tidy_df = tidy_scores(scores)
                records = df_to_json_records(tidy_df)
                self.send_json({"scores": records, "variant": var_str, "engine": "model"})

        except Exception as e:
            err_str = str(e)
            m = re.search(r"reference base does not match the expected reference base:\s*([ACGTN]+)", err_str, re.IGNORECASE)
            if m:
                expected_ref = m.group(1).upper()
                suggested_alt = var_obj.reference_bases if var_obj.reference_bases != expected_ref else var_obj.alternate_bases
                suggested_var = f"{var_obj.chromosome}:{var_obj.position}:{expected_ref}>{suggested_alt}"
                self.send_error_json(
                    f"Referans baz uyumsuzluğu: {var_obj.chromosome}:{var_obj.position} konumunda hg38 referans bazı '{expected_ref}' bekleniyor (girilen: '{var_obj.reference_bases}').",
                    status=400,
                    error_en=f"Reference base mismatch: At {var_obj.chromosome}:{var_obj.position}, hg38 reference base is '{expected_ref}' (entered: '{var_obj.reference_bases}').",
                    error_type="ref_mismatch",
                    expected_ref=expected_ref,
                    entered_ref=var_obj.reference_bases,
                    suggested_variant=suggested_var,
                )
                return
            self.send_error_json(f"Sorgu hatası: {err_str}", status=400, error_en=f"Query error: {err_str}")

    def handle_api_query_interval(self, data):
        int_str = data.get("interval", "").strip()
        tissue = data.get("tissue", "").strip() or None
        gene = data.get("gene", "").strip() or None

        if not int_str:
            self.send_error_json("Aralık belirtilmedi")
            return

        try:
            int_obj = cli.parse_interval(int_str)
        except Exception as e:
            self.send_error_json(f"Geçersiz aralık formatı: {e}")
            return

        class DummyArgs:
            api_key = None
            no_cache = False
            cache_dir = None

        try:
            atlas_client = cli.create_atlas(DummyArgs())
            terms = cli.resolve_tissue_term(tissue) if tissue else None
            genes = [gene] if gene else None

            raw_scores = atlas_client.query_interval(
                int_obj,
                ontology_terms=terms,
                gene_names=genes,
                progress_bar=False,
            )
            df = cli.tidy_atlas_scores(raw_scores)
            records = df_to_json_records(df)
            self.send_json({"scores": records, "interval": int_str})
        except Exception as e:
            self.send_error_json(f"Bölge sorgu hatası: {e}")

    def handle_api_query_batch(self, data):
        file_name = data.get("fileName", "data.tsv")
        file_content = data.get("fileContent", "")
        engine = data.get("engine", "atlas")
        tissue = data.get("tissue", "").strip() or None
        gene = data.get("gene", "").strip() or None

        if not file_content:
            self.send_error_json("Dosya içeriği boş")
            return

        # Geçici dosyaya yazıp CLI yardımcılarıyla oku
        with tempfile.NamedTemporaryFile(mode="w", suffix=Path(file_name).suffix, delete=False) as tmp:
            tmp.write(file_content)
            tmp_path = Path(tmp.name)

        try:
            if file_name.lower().endswith((".vcf", ".vcf.gz")):
                variants = cli.read_variants_from_vcf(tmp_path)
            else:
                variants = cli.read_variants_from_tsv_or_csv(tmp_path)
        except Exception as e:
            tmp_path.unlink(missing_ok=True)
            self.send_error_json(f"Dosya okuma hatası: {e}")
            return
        finally:
            tmp_path.unlink(missing_ok=True)

        if not variants:
            self.send_error_json("Dosyada geçerli varyant bulunamadı")
            return

        class DummyArgs:
            api_key = None
            no_cache = False
            cache_dir = None

        try:
            if engine == "atlas":
                atlas_client = cli.create_atlas(DummyArgs())
                terms = cli.resolve_tissue_term(tissue) if tissue else None
                genes = [gene] if gene else None

                raw_scores = atlas_client.query_variants(
                    variants,
                    ontology_terms=terms,
                    gene_names=genes,
                    progress_bar=False,
                )
                df = cli.tidy_atlas_scores(raw_scores)
                records = df_to_json_records(df)
                avi_summary = cli.extract_atlas_summary(raw_scores)
                self.send_json({"scores": records, "count": len(variants), "avi_summary": avi_summary})
            else:
                # Model inference batch
                model = cli.create_model(DummyArgs())
                intervals = [cli.get_variant_interval(v) for v in variants]
                from alphagenome.models.variant_scorers import get_recommended_scorers, tidy_scores
                from alphagenome.models.dna_model import Organism

                scorers = get_recommended_scorers(Organism.HOMO_SAPIENS)
                all_dfs = []
                for interval, var in zip(intervals, variants):
                    scores = model.score_variant(interval=interval, variant=var, variant_scorers=scorers)
                    all_dfs.append(tidy_scores(scores))

                import pandas as pd
                combined_df = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()
                records = df_to_json_records(combined_df)
                self.send_json({"scores": records, "count": len(variants)})

        except Exception as e:
            self.send_error_json(f"Toplu analiz hatası: {e}")


def start_web_server(host="127.0.0.1", port=8000, open_browser=True):
    """Web sunucusunu başlatır. Port meşgulse alternatif port bulur veya mevcut oturumu açar."""
    import urllib.request

    max_tries = 10
    httpd = None
    actual_port = port

    for p in range(port, port + max_tries):
        try:
            server_address = (host, p)
            httpd = ThreadingHTTPServer(server_address, AlphaGenomeWebHandler)
            actual_port = p
            break
        except OSError as e:
            if getattr(e, "errno", None) == 48 or "Address already in use" in str(e):
                # Halihazırda AlphaGenome mu çalışıyor kontrol et
                try:
                    with urllib.request.urlopen(f"http://{host}:{p}/api/status", timeout=1.0) as resp:
                        if resp.status == 200:
                            url = f"http://{host}:{p}"
                            print("\n" + "=" * 64)
                            print("  ℹ️  AlphaGenome Web Studio zaten çalışıyor!")
                            print(f"  🌐 Tarayıcınız açılıyor: {url}")
                            print("=" * 64 + "\n")
                            if open_browser:
                                webbrowser.open(url)
                            return
                except Exception:
                    pass
                continue
            raise e

    if httpd is None:
        print(f"Hata: {port}-{port+max_tries} aralığında uygun bir port bulunamadı.")
        return

    url = f"http://{host}:{actual_port}"
    print("\n" + "=" * 64)
    print("  🚀 AlphaGenome Web Studio Başlatıldı!")
    print(f"  🌐 Adres: {url}")
    print("  📁 Sürükle & Bırak: VCF, TSV, CSV, TXT desteklenir")
    print("  ⚡ Atlas Motoru: 9 milyar önceden hesaplanmış SNV hazır")
    print("  🛑 Kapatmak için: Bu pencerede Ctrl + C tuşlarına basın")
    print("=" * 64 + "\n")

    if open_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nWeb sunucusu kapatılıyor...")
        httpd.shutdown()
        httpd.server_close()
        print("Kapatıldı.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="AlphaGenome Web UI Sunucusu")
    parser.add_argument("--host", default="127.0.0.1", help="Host adresi (varsayılan: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8000, help="Port numarası (varsayılan: 8000)")
    parser.add_argument("--no-browser", action="store_true", help="Web tarayıcısını otomatik açma")
    args = parser.parse_args()

    start_web_server(host=args.host, port=args.port, open_browser=not args.no_browser)
