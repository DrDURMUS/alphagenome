# AlphaGenome CLI — Kapsamlı Kullanım Kılavuzu

> **AlphaGenome**, Google DeepMind tarafından geliştirilen, tek nükleotid düzeyinde çoklu genomik modaliteleri tahmin edebilen bir yapay zeka modelidir. Bu kılavuz, AlphaGenome API'sine kolay erişim sağlayan komut satırı aracının detaylı kullanımını kapsar.

---

## İçindekiler

1. [Kurulum & Ön Hazırlık](#1-kurulum--ön-hazırlık)
2. [Hızlı Başlangıç](#2-hızlı-başlangıç)
3. [Alt Komutlar](#3-alt-komutlar)
   - [3.1 `variant` — Tek Varyant Analizi](#31-variant--tek-varyant-analizi)
   - [3.2 `vcf` — VCF'den Toplu Analiz](#32-vcf--vcfden-toplu-analiz)
   - [3.3 `region` — Bölge Tahmini](#33-region--bölge-tahmini)
   - [3.4 `list` — Dosyadan Toplu Analiz](#34-list--dosyadan-toplu-analiz)
   - [3.5 `tissues` — Doku/Ontoloji Arama](#35-tissues--dokuontoloji-arama)
   - [3.6 `visualize` — Görselleştirme](#36-visualize--görselleştirme)
4. [Desteklenen Modaliteler](#4-desteklenen-modaliteler)
5. [Doku / Ontoloji Sistemi](#5-doku--ontoloji-sistemi)
6. [Giriş Formatları](#6-giriş-formatları)
7. [Çıktı Formatları](#7-çıktı-formatları)
8. [Gelişmiş Örnek Çalıştırmalar](#8-gelişmiş-örnek-çalıştırmalar)
10. [API Referansı](#10-api-referansı)
11. [AlphaGenome Atlas (9 Milyar SNV Veritabanı)](#11-alphagenome-atlas-9-milyar-snv-veritabanı)
12. [AlphaGenome Studio Web Arayüzü (Drag & Drop)](#12-alphagenome-studio-web-arayüzü-drag--drop)

---

## 1. Kurulum & Ön Hazırlık

### Sistem Gereksinimleri

| Bileşen | Minimum |
|---------|---------|
| Python | 3.11+ |
| RAM | 4 GB |
| İnternet | API erişimi için gerekli |

### Ortam Kurulumu

```bash
# 1. Conda ortamını aktifleştir
conda activate alphagenome

# 2. CLI aracının bulunduğu dizine geç
cd ~/alphagenome

# 3. Doğrulama — yardım menüsü
python alphagenome_cli.py --help
```

### API Anahtarı Yapılandırması

API anahtarı 3 yöntemle sağlanabilir (öncelik sırasına göre):

```bash
# Yöntem 1: Komut satırı argümanı (en yüksek öncelik)
python alphagenome_cli.py variant chr22:36201698:A>C --api-key YOUR_KEY

# Yöntem 2: Ortam değişkeni
export ALPHAGENOME_API_KEY="YOUR_KEY"
python alphagenome_cli.py variant chr22:36201698:A>C

# Yöntem 3: Yerleşik anahtar (varsayılan, kod içinde tanımlı)
python alphagenome_cli.py variant chr22:36201698:A>C
```

> [!TIP]
> Kalıcı yapılandırma için ortam değişkenini `~/.zshrc` dosyanıza ekleyin:
> ```bash
> echo 'export ALPHAGENOME_API_KEY="YOUR_KEY"' >> ~/.zshrc
> source ~/.zshrc
> ```

---

## 2. Hızlı Başlangıç

```bash
conda activate alphagenome
cd ~/alphagenome

# ✅ Tek bir varyantın etkisini hemen analiz et
python alphagenome_cli.py variant chr22:36201698:A>C --tissue colon

# ✅ Mevcut dokuları ara
python alphagenome_cli.py tissues --search liver

# ✅ Sonuçları görselleştir
python alphagenome_cli.py visualize variant chr22:36201698:A>C --tissue colon
```

---

## 3. Alt Komutlar

### 3.1 `variant` — Tek Varyant Analizi

Tek bir genetik varyantın belirtilen doku üzerindeki etkisini analiz eder.

#### Sözdizimi

```
python alphagenome_cli.py variant VARYANT [seçenekler]
```

#### Parametreler

| Parametre | Kısa | Zorunlu | Açıklama |
|-----------|-------|---------|----------|
| `VARYANT` | — | ✅ | Varyant tanımı (aşağıdaki formatlara bakın) |
| `--tissue` | `-t` | — | Doku adı veya UBERON/CL ontoloji terimi |
| `--modality` | `-m` | — | Çıktı tipi (varsayılan: `RNA_SEQ`) |
| `--score` | `-s` | — | Detaylı varyant skorlaması yapılsın mı |
| `--output` | `-o` | — | Sonuçları dosyaya kaydet (.tsv veya .csv) |

#### Varyant Formatları

```
chr22:36201698:A>C      ← Tercih edilen format
chr22:36201698:A:C      ← İki nokta ayırıcı
22:36201698:A>C         ← chr ön eki opsiyonel
```

#### Örnekler

```bash
# Temel analiz — kolon dokusunda RNA ifadesi
python alphagenome_cli.py variant chr22:36201698:A>C --tissue colon

# Çoklu modality ile analiz
python alphagenome_cli.py variant chr22:36201698:A>C \
  --tissue colon --modality RNA_SEQ CAGE ATAC

# Detaylı skorlama + dosyaya kaydetme
python alphagenome_cli.py variant chr22:36201698:A>C \
  --tissue colon --score --output sonuc.tsv

# Ontoloji terimi ile doğrudan belirtme (tek doku)
python alphagenome_cli.py variant chr1:155235218:G>T \
  --tissue UBERON:0002107 --modality RNA_SEQ
```

#### Çıktı Açıklaması

```
📊 Sonuçlar:
  RNA_SEQ:
    REF aralığı : [0.0000, 43.0000]    ← Referans aleldeki değer aralığı
    ALT aralığı : [0.0000, 43.0000]    ← Alternatif aleldeki değer aralığı
    Max |etki|  : 1.500000              ← Maksimum mutlak etki büyüklüğü
    Ort |etki|  : 0.000399              ← Ortalama mutlak etki

🏆 En yüksek skorlar:                   (--score verilmişse)
  gene_name  scorer_name  score
  APOL4      RNA_SEQ      -0.423
```

---

### 3.2 `vcf` — VCF'den Toplu Analiz

Standart VCF formatındaki dosyalardan tüm varyantları toplu olarak skorlar.

#### Sözdizimi

```
python alphagenome_cli.py vcf VCF_DOSYASI [seçenekler]
```

#### Parametreler

| Parametre | Kısa | Zorunlu | Açıklama |
|-----------|-------|---------|----------|
| `VCF_DOSYASI` | — | ✅ | VCF dosya yolu (.vcf veya .vcf.gz) |
| `--tissue` | `-t` | — | Doku adı |
| `--output` | `-o` | — | Çıktı dosyası (varsayılan: `<vcf>_scores.tsv`) |
| `--workers` | `-w` | — | Paralel iş parçacığı sayısı (varsayılan: 5) |

#### Örnekler

```bash
# Temel toplu analiz
python alphagenome_cli.py vcf my_variants.vcf --tissue brain

# Sıkıştırılmış VCF + özel çıktı
python alphagenome_cli.py vcf sample.vcf.gz \
  --tissue liver --output liver_scores.tsv

# 10 paralel işçi ile hızlı analiz
python alphagenome_cli.py vcf large_cohort.vcf \
  --tissue colon --workers 10 --output cohort_results.tsv
```

#### Desteklenen VCF Formatı

```
#CHROM  POS     ID        REF  ALT  QUAL  FILTER  INFO
chr22   36201698  rs123   A    C    .     PASS    .
chr22   36265860  .       G    A    .     .       .
```

- `#` ile başlayan satırlar otomatik atlanır
- Çoklu ALT alleller (`,` ile ayrılmış) ayrı ayrı işlenir
- `chr` ön eki opsiyoneldir — otomatik eklenir

---

### 3.3 `region` — Bölge Tahmini

Belirtilen genomik bölge için ileri yönlü (forward) tahmin yapar. Varyant karşılaştırması olmadan, sadece referans genom üzerinde modelin ne öngördüğünü gösterir.

#### Sözdizimi

```
python alphagenome_cli.py region BÖLGE [seçenekler]
```

#### Parametreler

| Parametre | Kısa | Zorunlu | Açıklama |
|-----------|-------|---------|----------|
| `BÖLGE` | — | ✅ | Genomik bölge (chr:start-end) |
| `--tissue` | `-t` | — | Doku adı |
| `--modality` | `-m` | — | Çıktı tipleri (varsayılan: `RNA_SEQ`) |
| `--output` | `-o` | — | Özet çıktı dosyası |

> [!IMPORTANT]
> Bölge genişliği AlphaGenome modeli tarafından **1,048,576 bp (1 Mb)** olarak sabitlenmiştir. Farklı genişlikteki bölgeler otomatik ayarlanır.

#### Örnekler

```bash
# Tek modality ile bölge tahmini
python alphagenome_cli.py region chr22:35677410-36725986 --tissue liver

# Çoklu modality
python alphagenome_cli.py region chr22:35677410-36725986 \
  --tissue colon --modality RNA_SEQ ATAC DNASE

# Tüm ana modaliteler ile kapsamlı analiz
python alphagenome_cli.py region chr22:35677410-36725986 \
  --tissue colon \
  --modality RNA_SEQ CAGE ATAC DNASE SPLICE_SITES \
  --output region_full.tsv
```

---

### 3.4 `list` — Dosyadan Toplu Analiz

TSV veya CSV formatındaki varyant listelerini toplu olarak skorlar.

#### Sözdizimi

```
python alphagenome_cli.py list DOSYA [seçenekler]
```

#### Parametreler

| Parametre | Kısa | Zorunlu | Açıklama |
|-----------|-------|---------|----------|
| `DOSYA` | — | ✅ | Varyant listesi dosyası (.tsv veya .csv) |
| `--tissue` | `-t` | — | Doku adı |
| `--output` | `-o` | — | Çıktı dosyası (varsayılan: `<dosya>_scores.tsv`) |
| `--workers` | `-w` | — | Paralel iş parçacığı sayısı (varsayılan: 5) |

#### Dosya Formatı

**TSV (Tab ile ayrılmış):**
```
chr	pos	ref	alt
chr22	36201698	A	C
chr22	36265860	G	A
chr1	155235218	G	T
```

**CSV (Virgül ile ayrılmış):**
```
chr,pos,ref,alt
chr22,36201698,A,C
chr22,36265860,G,A
```

> [!TIP]
> - Başlık satırı opsiyoneldir — otomatik algılanır
> - `chr` ön eki opsiyoneldir
> - `#` ile başlayan satırlar yorum olarak atlanır

#### Örnekler

```bash
# Temel toplu analiz
python alphagenome_cli.py list my_variants.tsv --tissue colon

# CSV dosyasından
python alphagenome_cli.py list my_variants.csv \
  --tissue brain --output brain_scores.csv

# Yüksek paralellik ile büyük listeler
python alphagenome_cli.py list large_list.tsv \
  --tissue liver --workers 10 --output liver_scores.tsv
```

#### Çıktı Tablosu

Skorlama sonuçları detaylı bir tablo olarak kaydedilir:

| Sütun | Açıklama |
|-------|----------|
| `variant_id` | Varyant tanımlayıcısı |
| `gene_name` | Etkilenen gen adı |
| `scorer_name` | Kullanılan skorlama yöntemi |
| `score` | Etki skoru (negatif = azalma, pozitif = artış) |
| `output_type` | Modalite tipi |
| `strand` | DNA iplikçiği (+/-) |

---

### 3.5 `tissues` — Doku/Ontoloji Arama

AlphaGenome'da mevcut olan tüm doku ve hücre tiplerini listeler veya arar. Diğer komutlarda `--tissue` parametresine ne vereceğinizi öğrenmek için kullanın.

#### Sözdizimi

```
python alphagenome_cli.py tissues [seçenekler]
```

#### Parametreler

| Parametre | Kısa | Açıklama |
|-----------|-------|----------|
| `--search` | `-s` | Aranacak kelime (büyük/küçük harf duyarsız) |

#### Örnekler

```bash
# Karaciğer dokularını ara
python alphagenome_cli.py tissues --search liver

# Beyin bölgelerini listele
python alphagenome_cli.py tissues --search brain

# Böbrek dokuları
python alphagenome_cli.py tissues --search kidney

# Hücre tipleri ara
python alphagenome_cli.py tissues --search fibroblast

# Tüm dokuları listele (arama yapmadan)
python alphagenome_cli.py tissues
```

#### Örnek Çıktı

```
'liver' araması — 5 sonuç:

  UBERON:0001115       left lobe of liver       (tissue, adult)
  UBERON:0002107       liver                    (tissue, embryonic)
  UBERON:0002107       liver                    (tissue, child)
  UBERON:0002107       liver                    (tissue, child,adult)
  UBERON:0001114       right lobe of liver      (tissue, adult)
```

> [!TIP]
> Aramada bulunan ontoloji terimlerini (`UBERON:XXXXXXX`) doğrudan `--tissue` parametresi olarak kullanabilirsiniz:
> ```bash
> python alphagenome_cli.py variant chr22:36201698:A>C --tissue UBERON:0002107
> ```

---

### 3.6 `visualize` — Görselleştirme

AlphaGenome tahminlerini GENCODE gen anotasyonları ile birlikte yayın kalitesinde grafik olarak görselleştirir.

#### Sözdizimi

```
python alphagenome_cli.py visualize MOD HEDEF [seçenekler]
```

#### Modlar

| Mod | Açıklama |
|-----|----------|
| `variant` | Varyant etkisi — REF vs ALT karşılaştırmalı |
| `region` | Bölge tahmini — referans genom üzerinde |

#### Parametreler

| Parametre | Kısa | Zorunlu | Açıklama |
|-----------|-------|---------|----------|
| `MOD` | — | ✅ | `variant` veya `region` |
| `HEDEF` | — | ✅ | Varyant veya bölge tanımı |
| `--tissue` | `-t` | — | Doku adı |
| `--modality` | `-m` | — | Çıktı tipleri (varsayılan: `RNA_SEQ`) |
| `--output` | `-o` | — | Çıktı dosyası (.png/.pdf/.svg) |
| `--zoom` | `-z` | — | Zoom penceresi boyutu (bp) |

#### Örnekler

```bash
# Varyant etkisi — RNA ifadesi
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --output variant_effect.png

# Çoklu modality ile görselleştirme
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --modality RNA_SEQ CAGE --output variant_multi.png

# Yakınlaştırılmış görünüm (8kb pencere)
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --zoom 8000 --output variant_zoomed.png

# Bölge tahmini görselleştirme
python alphagenome_cli.py visualize region chr22:35677410-36725986 \
  --tissue colon --modality RNA_SEQ ATAC --output region.png

# PDF olarak kaydetme (vektörel kalite)
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --output variant_hires.pdf

# SVG olarak kaydetme (düzenlenebilir)
python alphagenome_cli.py visualize region chr22:35677410-36725986 \
  --tissue liver --output region.svg
```

#### Grafik Özellikleri

**Variant modu:**
- 🧬 GENCODE v46 gen anotasyonları (protein-kodlayan genler)
- 📊 REF (gri) ve ALT (kırmızı) sinyal karşılaştırması
- 📍 Varyant konumu turuncu dikey çizgi ile işaretlenir
- 🔍 Varsayılan zoom: ±16 kb (32,768 bp)

**Region modu:**
- 🧬 GENCODE v46 gen anotasyonları
- 📊 Her modality için ayrı sinyal track'leri
- 🌈 Viridis renk paleti ile renklendirilmiş track'ler

---

## 4. Desteklenen Modaliteler

AlphaGenome 11 farklı genomik modaliteyi tahmin edebilir:

| Modalite | Açıklama | Çözünürlük |
|----------|----------|------------|
| `RNA_SEQ` | RNA ekspresyonu (gen ifadesi) | 1 bp |
| `CAGE` | Transkripsiyon başlangıç bölgeleri | 1 bp |
| `ATAC` | Açık kromatin bölgeleri (ATAC-seq) | 1 bp |
| `DNASE` | DNase aşırı duyarlı bölgeler | 1 bp |
| `PROCAP` | PRO-cap transkripsiyon başlangıcı | 1 bp |
| `CHIP_HISTONE` | Histon modifikasyonları (H3K4me3 vb.) | 128 bp |
| `CHIP_TF` | Transkripsiyon faktörü bağlanması | 128 bp |
| `SPLICE_SITES` | Splice (eklem) bölgeleri | 1 bp |
| `SPLICE_JUNCTIONS` | Splice bileşimleri (junction'lar) | — |
| `SPLICE_SITE_USAGE` | Splice bölgesi kullanım oranları | 1 bp |
| `CONTACT_MAPS` | 3D kromatin temas haritaları | — |

> [!NOTE]
> - `SPLICE_SITES` doku bağımsızdır — `--tissue` parametresi etki etmez
> - `CHIP_HISTONE` ve `CHIP_TF` 128 bp çözünürlükte çalışır
> - Çoklu modality belirtmek için: `--modality RNA_SEQ CAGE ATAC`

### Modalite Seçim Rehberi

| Araştırma Sorusu | Önerilen Modaliteler |
|-----------------|---------------------|
| Gen ifadesi değişimi | `RNA_SEQ`, `CAGE` |
| Düzenleyici element aktifliği | `ATAC`, `DNASE`, `CHIP_HISTONE` |
| Splice bölgesi etkisi | `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS` |
| Promotör aktivitesi | `CAGE`, `PROCAP`, `CHIP_HISTONE` |
| Transkripsiyon faktörü bağlanması | `CHIP_TF` |
| Kromatin yapısı | `CONTACT_MAPS`, `ATAC` |
| Kapsamlı varyant analizi | `RNA_SEQ CAGE ATAC SPLICE_SITES` |

---

## 5. Doku / Ontoloji Sistemi

### Otomatik Çoklu Doku Eşleşmesi

`--tissue` parametresi ile yapılan arama sonucu **tüm eşleşen dokular** otomatik olarak kullanılır. Her doku ayrı track olarak gösterilir:

```bash
python alphagenome_cli.py variant chr22:36201698:A>C --tissue colon
```

```
✅ 5 doku eşleşmesi bulundu — tümü kullanılacak:
ℹ️    → colonic mucosa (UBERON:0000317)
ℹ️    → transverse colon (UBERON:0001157)
ℹ️    → sigmoid colon (UBERON:0001159)
ℹ️    → mucosa of descending colon (UBERON:0004992)
ℹ️    → left colon (UBERON:0008971)
```

> [!TIP]
> Tek bir doku istiyorsanız doğrudan UBERON kodu ile belirtin:
> ```bash
> --tissue UBERON:0001159   # Sadece sigmoid kolon
> ```

### Doku Belirtme Yöntemleri

```bash
# 1. Doğal dil ile (tüm eşleşmeler kullanılır)
--tissue colon        # → 5 kolon dokusu
--tissue liver        # → 5 karaciğer dokusu
--tissue brain        # → 16 beyin bölgesi

# 2. UBERON ontoloji terimi ile (tek doku, kesin)
--tissue UBERON:0002107      # Karaciğer
--tissue UBERON:0001159      # Kolon - sigmoid
--tissue UBERON:0000955      # Beyin

# 3. Hücre tipi ontoloji terimi ile
--tissue CL:0000057          # Fibroblast
```

### Sık Kullanılan Dokular

| Doku | Ontoloji Terimi | Kısa Ad | Eşleşme Sayısı |
|------|-----------------|---------|----------------|
| Karaciğer | `UBERON:0002107` | `liver` | 5 |
| Kolon (sigmoid) | `UBERON:0001159` | `colon` | 5 |
| Beyin | `UBERON:0000955` | `brain` | 16 |
| Akciğer | `UBERON:0002048` | `lung` | — |
| Kalp | `UBERON:0000948` | `heart` | — |
| Böbrek | `UBERON:0002113` | `kidney` | — |
| Kas | `UBERON:0001134` | `muscle` | — |

---

## 6. Giriş Formatları

### Varyant Formatları

```
chr22:36201698:A>C         # Standart
chr22:36201698:A:C         # : ayırıcı
22:36201698:A>C            # chr ön eki olmadan
chrX:15561033:G>A          # X kromozomu
chr1:155235218:ATCG>A      # Delesyon (çoklu baz REF)
chr7:117559590:A>ATCG      # İnsersiyon (çoklu baz ALT)
```

### Bölge Formatları

```
chr22:35677410-36725986    # Standart
22:35677410-36725986       # chr ön eki olmadan
```

### TSV Dosya Formatı

```tsv
chr	pos	ref	alt
chr22	36201698	A	C
chr22	36265860	G	A
chr1	155235218	G	T
```

Kabul edilen başlık isimleri: `chr`, `chrom`, `chromosome`, `#chrom`, `#chr`

### VCF Dosya Formatı

Standart VCF v4.x formatı desteklenir:
- `.vcf` veya `.vcf.gz` (sıkıştırılmış)
- Çoklu ALT alleller (virgülle ayrılmış) ayrı varyantlar olarak işlenir
- `*` ve `.` ALT allelleri atlanır

---

## 7. Çıktı Formatları

### Tablo Çıktıları (.tsv / .csv)

```bash
--output sonuclar.tsv    # Tab ile ayrılmış
--output sonuclar.csv    # Virgül ile ayrılmış
```

### Görsel Çıktılar (.png / .pdf / .svg)

```bash
--output grafik.png      # Raster (200 DPI)
--output grafik.pdf      # Vektörel, yayın kalitesi
--output grafik.svg      # Düzenlenebilir vektörel
```

> [!NOTE]
> `visualize` komutuna `.tsv` gibi tablo uzantısı verilirse otomatik olarak `.png` eklenir.

---

## 8. Gelişmiş Örnek Çalıştırmalar

Aşağıdaki tüm komutlar test edilmiş ve doğrulanmıştır.

---

### 8.1 Tek Varyant — Temel Analiz

APOL4 genindeki bir varyantın tüm kolon dokularındaki etkisi:

```bash
python alphagenome_cli.py variant chr22:36201698:A>C --tissue colon
```

**Beklenen çıktı:**
```
✅ 5 doku eşleşmesi bulundu — tümü kullanılacak:
ℹ️    → colonic mucosa (UBERON:0000317)
ℹ️    → transverse colon (UBERON:0001157)
ℹ️    → sigmoid colon (UBERON:0001159)
ℹ️    → mucosa of descending colon (UBERON:0004992)
ℹ️    → left colon (UBERON:0008971)

📊 Sonuçlar:
  RNA_SEQ:
    REF aralığı : [0.0000, 43.0000]
    ALT aralığı : [0.0000, 43.0000]
    Max |etki|  : 1.500000
    Ort |etki|  : 0.000399
✅ Tekli varyant analizi tamamlandı!
```

---

### 8.2 Tek Varyant — Detaylı Skorlama + Dosyaya Kaydetme

Gen bazlı skorlama ile analiz etme ve TSV'ye kaydetme:

```bash
python alphagenome_cli.py variant chr22:36201698:A>C \
  --tissue colon --score --output apol4_scores.tsv
```

**Beklenen çıktı:**
```
📊 Sonuçlar:
  RNA_SEQ:
    Max |etki|  : 1.500000

🏆 En yüksek skorlar:
  gene_name    scorer_name    score
  APOL4        RNA_SEQ       -0.423
  APOL1        RNA_SEQ        0.112
✅ Sonuçlar kaydedildi: apol4_scores.tsv
```

---

### 8.3 Tek Varyant — Çoklu Modality

RNA ifadesi, CAGE ve kromatin erişilebilirliği üzerindeki etkiyi aynı anda görmek:

```bash
python alphagenome_cli.py variant chr22:36201698:A>C \
  --tissue colon --modality RNA_SEQ CAGE ATAC
```

---

### 8.4 Tek Varyant — Tek Doku (UBERON Kodu ile)

Sadece sigmoid kolon dokusunda, diğer eşleşmeleri dahil etmeden analiz:

```bash
python alphagenome_cli.py variant chr22:36201698:A>C \
  --tissue UBERON:0001159 --modality RNA_SEQ
```

---

### 8.5 Bölge Tahmini — Çoklu Modality

APOL gen kümesinin bulunduğu 1 Mb bölgede RNA ve ATAC tahminleri:

```bash
python alphagenome_cli.py region chr22:35677410-36725986 \
  --tissue colon --modality RNA_SEQ ATAC
```

**Beklenen çıktı:**
```
📊 Sonuçlar:
  RNA_SEQ:
    Çıktı aralığı: chr22:35677410-36725986:.
    Veri boyutu  : (1048576, 2)
    Değer aralığı: [0.0000, 21.5000]
    Ortalama     : 0.0733
  ATAC:
    Çıktı aralığı: chr22:35677410-36725986:.
    Veri boyutu  : (1048576, 1)
    Değer aralığı: [0.0000, 41.5000]
    Ortalama     : 0.0619
✅ Bölge analizi tamamlandı!
```

---

### 8.6 Bölge Tahmini — Kapsamlı Modaliteler

Karaciğerde tüm temel modaliteleri tek seferde taramak:

```bash
python alphagenome_cli.py region chr22:35677410-36725986 \
  --tissue liver \
  --modality RNA_SEQ CAGE ATAC DNASE SPLICE_SITES \
  --output region_comprehensive.tsv
```

---

### 8.7 Toplu Analiz — TSV Dosyasından

`example_variants.tsv` dosyasından 3 varyantın toplu skorlaması:

```bash
python alphagenome_cli.py list example_variants.tsv \
  --tissue colon --output results.tsv
```

**Beklenen çıktı:**
```
✅ 3 varyant okundu
ℹ️  Varyantlar skorlanıyor...
100%|████████████████| 3/3 [00:02<00:00]
✅ Toplam skor: 138686 satır × 25 sütun

📊 Özet:
  Benzersiz varyant sayısı: 3
  Etkilenen gen sayısı: 103
✅ Sonuçlar kaydedildi: results.tsv
```

---

### 8.8 Toplu Analiz — VCF Dosyasından

VCF dosyasından paralel skorlama (10 işçi ile):

```bash
python alphagenome_cli.py vcf patient_variants.vcf \
  --tissue brain --workers 10 --output brain_scores.tsv
```

---

### 8.9 Doku Arama — Karaciğer

```bash
python alphagenome_cli.py tissues --search liver
```

**Beklenen çıktı:**
```
'liver' araması — 5 sonuç:

  UBERON:0001115       left lobe of liver       (tissue, adult)
  UBERON:0002107       liver                    (tissue, embryonic)
  UBERON:0002107       liver                    (tissue, child)
  UBERON:0002107       liver                    (tissue, child,adult)
  UBERON:0001114       right lobe of liver      (tissue, adult)
```

---

### 8.10 Doku Arama — Beyin

```bash
python alphagenome_cli.py tissues --search brain
```

**Beklenen çıktı:**
```
'brain' araması — 16 sonuç:

  UBERON:0001954       Ammon's horn                        (tissue, adult)
  UBERON:0006469       C1 segment of cervical spinal cord  (tissue, adult)
  UBERON:0001876       amygdala                            (tissue, adult)
  UBERON:0009835       anterior cingulate cortex           (tissue, adult)
  UBERON:0000955       brain                               (tissue, adult)
  UBERON:0001873       caudate nucleus                     (tissue, adult)
  ...
```

---

### 8.11 Görselleştirme — Varyant Etkisi (REF vs ALT)

Gen anotasyonlu, REF/ALT karşılaştırmalı varyant etkisi grafiği:

```bash
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --modality RNA_SEQ CAGE \
  --output variant_effect.png
```

**Beklenen çıktı:**
```
✅ 5 doku eşleşmesi bulundu — tümü kullanılacak:
ℹ️  GENCODE gen anotasyonları yükleniyor...
✅ 11 transkript bulundu
✅ Tahmin tamamlandı!
✅ Görsel kaydedildi: variant_effect.png
```

**Grafik içeriği:**
- Üst panel: GENCODE gen anotasyonları (APOL4 ve çevresi)
- Orta/Alt paneller: Her modality için REF (gri) vs ALT (kırmızı) overlay
- Turuncu dikey çizgi: varyant konumu işaretçisi

---

### 8.12 Görselleştirme — Yakınlaştırılmış Görünüm

Varyant çevresine 8 kb zoom yaparak detaylı inceleme:

```bash
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --zoom 8000 --output zoomed_view.png
```

---

### 8.13 Görselleştirme — Bölge Tahmini

1 Mb bölgede gen ifadesi ve kromatin erişilebilirliği haritası:

```bash
python alphagenome_cli.py visualize region chr22:35677410-36725986 \
  --tissue colon --modality RNA_SEQ ATAC \
  --output region_prediction.png
```

**Beklenen çıktı:**
```
✅ 11 transkript bulundu
✅ Tahmin tamamlandı!
✅ Görsel kaydedildi: region_prediction.png
```

**Grafik içeriği:**
- Gen yapıları: RBFOX2, APOL5, APOL3, APOL4, APOL2, APOL1, MYH9, FOXRED2, CACNG2...
- RNA_SEQ track'leri: + ve - iplikçik ayrımı ile
- ATAC track'leri: açık kromatin bölgeleri

---

### 8.14 Görselleştirme — Yayın Kalitesinde PDF/SVG

```bash
# PDF — dergi yayını kalitesinde vektörel grafik
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --modality RNA_SEQ CAGE ATAC DNASE \
  --output publication_figure.pdf

# SVG — Adobe Illustrator / Inkscape ile düzenlenebilir
python alphagenome_cli.py visualize region chr22:35677410-36725986 \
  --tissue liver --output editable_figure.svg
```

---

### 8.15 Pipeline: VCF → Skorlama → Görselleştirme

Klinik iş akışı — VCF'den en etkili varyantı bulup görselleştirmek:

```bash
# Adım 1: Toplu skorlama
python alphagenome_cli.py vcf patient.vcf \
  --tissue colon --output patient_scores.tsv --workers 10

# Adım 2: En yüksek skorlu varyantları bul
sort -t$'\t' -k4 -rn patient_scores.tsv | head -10

# Adım 3: İlginç varyantı görselleştir
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon --modality RNA_SEQ CAGE ATAC \
  --zoom 16000 --output top_hit.pdf
```

---

### 8.16 Çoklu Doku Karşılaştırması

Aynı varyantın farklı dokulardaki etkisini ayrı ayrı analiz etmek:

```bash
# UBERON kodu ile tek tek dokuları karşılaştır
for tissue in UBERON:0001159 UBERON:0002107 UBERON:0000955; do
  python alphagenome_cli.py variant chr22:36201698:A>C \
    --tissue $tissue --score \
    --output "scores_${tissue##*:}.tsv"
done

# Her doku için ayrı görselleştirme
for tissue in colon liver brain; do
  python alphagenome_cli.py visualize variant chr22:36201698:A>C \
    --tissue $tissue --modality RNA_SEQ \
    --output "viz_${tissue}.png"
done
```

---

### 8.17 Splice Etkisi Analizi

Bir varyantın mRNA splice mekanizması üzerindeki etkisini incelemek:

```bash
# Splice modaliteleri ile analiz
python alphagenome_cli.py variant chr22:36201698:A>C \
  --tissue colon \
  --modality SPLICE_SITES SPLICE_SITE_USAGE

# Splice etkisini görselleştir
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon \
  --modality RNA_SEQ SPLICE_SITES SPLICE_SITE_USAGE \
  --zoom 16000 --output splice_effect.png
```

---

### 8.18 Histon Modifikasyonu ve TF Bağlanma Analizi

Düzenleyici bölgelerdeki epigenetik değişiklik tahminleri:

```bash
# Histon modifikasyonları
python alphagenome_cli.py region chr22:35677410-36725986 \
  --tissue colon --modality CHIP_HISTONE \
  --output histone_marks.tsv

# Histon + TF birlikte görselleştir
python alphagenome_cli.py visualize region chr22:35677410-36725986 \
  --tissue colon --modality CHIP_HISTONE CHIP_TF \
  --output epigenetic_landscape.png
```

---

### 8.19 Kapsamlı Varyant Profili — Tüm Ana Modaliteler

Tek bir varyantın tüm ana genomik katmanlarındaki etkisini görmek:

```bash
python alphagenome_cli.py visualize variant chr22:36201698:A>C \
  --tissue colon \
  --modality RNA_SEQ CAGE ATAC DNASE CHIP_HISTONE SPLICE_SITES \
  --zoom 16000 \
  --output comprehensive_variant_profile.pdf
```

---

### 8.20 Büyük Ölçekli Toplu İşlem — Parçalama Stratejisi

Binlerce varyant içeren büyük dosyaları parçalayarak verimli işlemek:

```bash
# Varyant listesini 100'lük parçalara böl
split -l 100 large_variants.tsv chunk_

# Her parçayı yüksek paralellikle işle
for f in chunk_*; do
  python alphagenome_cli.py list "$f" \
    --tissue colon --workers 10 \
    --output "${f}_scores.tsv"
done

# Tüm sonuçları tek dosyada birleştir
head -1 chunk_aa_scores.tsv > all_scores.tsv
tail -n +2 -q chunk_*_scores.tsv >> all_scores.tsv
echo "Toplam satır: $(wc -l < all_scores.tsv)"
```

---

### 8.21 Bash Kısayolları

Sık kullanılan komutlar için `~/.zshrc` dosyanıza ekleyin:

```bash
alias ag='conda activate alphagenome && python ~/alphagenome/alphagenome_cli.py'
alias ag-var='ag variant'
alias ag-vis='ag visualize'
alias ag-tis='ag tissues'
alias ag-reg='ag region'
alias ag-vcf='ag vcf'
alias ag-lst='ag list'
```

Kullanım:
```bash
ag-var chr22:36201698:A>C --tissue colon
ag-vis variant chr22:36201698:A>C --tissue colon --output result.png
ag-tis --search liver
ag-reg chr22:35677410-36725986 --tissue liver --modality RNA_SEQ ATAC
```

---

## 9. Sorun Giderme

### Sık Karşılaşılan Hatalar

| Hata | Sebep | Çözüm |
|------|-------|-------|
| `Varyant formatı anlaşılamadı` | Yanlış format | `chr22:36201698:A>C` formatını kullanın |
| `Bölge formatı anlaşılamadı` | Yanlış bölge | `chr22:35677410-36725986` formatını kullanın |
| `Bilinmeyen modality` | Yanlış modality adı | Modalite adlarını kontrol edin (Bölüm 4) |
| `Eşleşen doku bulunamadı` | Arama terimi yetersiz | `tissues --search` ile doğru terimi arayın |
| `Format 'tsv' is not supported` | `visualize` komutuna tablo uzantısı | `.png`, `.pdf` veya `.svg` kullanın |
| `VCF dosyası bulunamadı` | Yanlış dosya yolu | Tam dosya yolunu kontrol edin |
| `API bağlantı hatası` | İnternet/API sorunu | İnternet ve API anahtarını kontrol edin |

### Conda Ortamı Sorunları

```bash
# Ortam aktif değilse
conda activate alphagenome

# Paket eksikse
pip install alphagenome tqdm matplotlib pandas numpy
```

### Bellek Sorunları

```bash
# İşçi sayısını azaltın
python alphagenome_cli.py list my_variants.tsv --tissue colon --workers 2

# Büyük dosyaları parçalayın (Bölüm 8.20'ye bakın)
```

---

## 10. API Referansı

### Temel AlphaGenome Python API

| Metot | Açıklama |
|-------|----------|
| `predict_interval()` | Bölge bazlı ileri tahmin |
| `predict_variant()` | Tek varyant etkisi tahmini |
| `predict_variants()` | Toplu varyant etkisi tahmini |
| `score_variant()` | Tek varyant skorlaması |
| `score_variants()` | Toplu varyant skorlaması |
| `output_metadata()` | Mevcut modalite ve doku bilgileri |

### Doğrudan Python Kullanımı

CLI aracı yeterli gelmediğinde, doğrudan Python API kullanarak özel iş akışları oluşturabilirsiniz:

```python
from alphagenome.data import genome
from alphagenome.models import dna_client

# Model oluştur
model = dna_client.create("YOUR_API_KEY")

# Bölge tahmini — çoklu doku
interval = genome.Interval("chr22", 35677410, 36725986)
output = model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.RNA_SEQ},
    ontology_terms=[
        "UBERON:0001159",  # Sigmoid kolon
        "UBERON:0001157",  # Transvers kolon
    ],
)

# Varyant etkisi — çoklu modality
variant = genome.Variant.from_str("chr22:36201698:A>C")
variant_output = model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs={
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.CAGE,
        dna_client.OutputType.ATAC,
    },
    ontology_terms=["UBERON:0001159"],
)

# Varyant skorlama
from alphagenome.models.variant_scorers import get_recommended_scorers, tidy_scores
scorers = get_recommended_scorers(dna_client.Organism.HOMO_SAPIENS)
scores = model.score_variant(interval, variant, variant_scorers=scorers)
df = tidy_scores(scores)
print(df.sort_values("score", ascending=False).head(10))

# Görselleştirme — Python API ile
from alphagenome.visualization import plot_components
fig = plot_components.plot(
    components=[
        plot_components.OverlaidTracks(
            tdata={"REF": variant_output.reference.rna_seq,
                   "ALT": variant_output.alternate.rna_seq},
            colors={"REF": "dimgrey", "ALT": "red"},
        ),
    ],
    interval=interval,
    title="Variant Effect: chr22:36201698 A>C",
    annotations=[plot_components.VariantAnnotation([variant])],
)
fig.savefig("custom_plot.png", dpi=200, bbox_inches="tight")
```

---

## 11. AlphaGenome Atlas (9 Milyar SNV Veritabanı)

### 11.1 AlphaGenome Atlas Nedir?

**AlphaGenome Atlas**, tüm insan genomundaki (~3.2 milyar baz çifti) olası tüm tek nükleotid değişimleri (SNV) için **önceden hesaplanmış 9 milyar varyantlık** devasa (1 Petabayt) düzenleyici genomik etki veritabanıdır.

| Özellik | Standart Model Çıkarımı | AlphaGenome Atlas |
|---------|-------------------------|-------------------|
| **Veri Kaynağı** | API üzerinden anlık derin öğrenme inferansı | Önceden hesaplanmış 1 PB bulut veritabanı |
| **Hız** | ~1-3 saniye / varyant | Milisaniyeler mertebesinde (anında) |
| **Kapsam** | İnsan & Fare | İnsan genomundaki tüm 9 milyar olası SNV |
| **Skor Türü** | Ham fark/oran skorları | Ham skorlar + **AVI Yüzdelik (Quantile)** skorları |
| **Bölge Taraması** | Her pencere için model çağrısı | Genomik aralıktaki tüm varyantları tek sorguda çekme |

### 11.2 Atlas Komutları

#### Tek Varyant Sorgulama (`atlas variant`)
```bash
# Temel sorgu
python alphagenome_cli.py atlas variant chr22:36201698:A>C

# Doku ve hedef gen filtreli
python alphagenome_cli.py atlas variant chr22:36201698:A>C --tissue colon --gene APOL4

# Klasik komuta --atlas bayrağı ile
python alphagenome_cli.py variant chr22:36201698:A>C --atlas
```

#### Genomik Aralık Sorgulama (`atlas interval`)
Belirtilen penceredeki tüm olası tek nükleotid varyantlarının Atlas skorlarını topluca çeker:
```bash
python alphagenome_cli.py atlas interval chr22:35677410-35678410 --tissue liver --output liver_region_atlas.tsv
```

#### Dosyadan Toplu Atlas Sorgusu (`atlas list` & `atlas vcf`)
VCF veya TSV/CSV dosyalarınızdaki varyantları model bekleme süresi olmadan doğrudan Atlas üzerinden paralel olarak skorlar:
```bash
# TSV/CSV listesi
python alphagenome_cli.py atlas list variants.tsv --tissue colon --output atlas_scores.tsv --workers 10

# VCF dosyası
python alphagenome_cli.py atlas vcf variants.vcf --output vcf_atlas_scores.tsv --workers 10

# Mevcut komutlarda --atlas bayrağı kullanımı:
python alphagenome_cli.py vcf variants.vcf --atlas --output scores.tsv
python alphagenome_cli.py list variants.tsv --atlas --output scores.tsv
```

#### Scorer Kataloğunu İnceleme (`atlas scorers`)
Atlas'ta tanımlı 30'u aşkın düzenleyici etki scorer'ını ve track sayılarını listeler:
```bash
python alphagenome_cli.py atlas scorers
python alphagenome_cli.py atlas scorers --search splice
```

---

## 12. AlphaGenome Studio Web Arayüzü (Drag & Drop)

### 12.1 Web Sunucusunu Başlatma

Modern ve kolay kullanımlı web arayüzünü iki farklı yöntemle başlatabilirsiniz:

```bash
# Yöntem 1: CLI üzerinden
python alphagenome_cli.py ui --port 8000

# Yöntem 2: Doğrudan web_app betiği ile
python web_app.py --port 8000
```

Varsayılan olarak `http://127.0.0.1:8000` adresinde açılır ve tarayıcınızı otomatik olarak yönlendirir.

### 12.2 Öne Çıkan Özellikler

1. **Sürükle-Bırak (Drag & Drop) Yükleme**:
   - `.vcf`, `.vcf.gz`, `.tsv`, `.csv` ve `.txt` formatındaki dosyalarınızı doğrudan sürükleyip bırakabilirsiniz.
   - Dosya yüklendiğinde anında ilk 5-10 satır önizleme tablosu ve tespit edilen geçerli varyant sayısı görüntülenir.
   - Hızlı deneme için **Örnek TSV Yükle** ve **Örnek VCF Yükle** butonları mevcuttur.

2. **Atlas ve Model İkili Motor Seçimi**:
   - **⚡ AlphaGenome Atlas**: 9 milyar SNV veritabanından anında sorgulama.
   - **🧬 Model Doğrudan İnferans**: Canlı model çalıştırarak skorlama.

3. **Biyolojik ve Klinik Karar Destek Özeti**:
   - Varyantın tahmini düzenleyici etki seviyesini (Belirgin Etki, Orta Düzey Etki, Zayıf Etki, İhmal Edilebilir) canlı olarak etiketler.
   - En çok etkilenen genleri, mekanizmaları (Splice kavşakları, Gene Mask LFC, Promotor merkez maskesi) özetler.

4. **İnteraktif Sonuç Tablosu**:
   - Kolon başlıklarına tıklayarak sıralama (Ham Skor, Yüzdelik AVI, Gen, Varyant).
   - Metin araması ve etki seviyesi filtrelemesi (🔴 Belirgin, 🟠 Orta, 🟡 Zayıf).
   - Tek tıkla **TSV İndir** veya **CSV İndir** dışa aktarımı.

5. **Atlas Scorer Kataloğu & Canlı Ayarlar**:
   - Tüm Atlas scorer'larını, yönlerini ve açıklamalarını aratabileceğiniz görsel katalog.
   - Sağ üstteki ayarlar simgesiyle API anahtarını çalışma anında güncelleme imkanı.

---

> **Sürüm:** AlphaGenome CLI & Web Studio v1.2 | AlphaGenome Atlas v0.9.0 | Python 3.11
>
> **Tarih:** 8 Eylül 2026

