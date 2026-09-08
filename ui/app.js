/**
 * AlphaGenome Studio — Modern Web UI Client
 * Multi-language (TR / EN), Atlas & Model querying, interactive table sorting,
 * filtering, biological interpretation, smart error auto-correction, and data export.
 */

// ─── Internationalization (i18n) Dictionary ────────────────────────────────
const i18n = {
  tr: {
    brand_subtitle: "Google DeepMind Düzenleyici Genomik Analiz Stüdyosu",
    status_online: "Sistem Hazır",
    status_offline: "Bağlantı Yok",
    status_connecting: "Bağlanıyor...",
    settings_title: "Ayarlar ve API Anahtarı",
    
    // Tabs
    tab_upload: "Sürükle-Bırak & Toplu Dosya (VCF/TSV)",
    tab_variant: "Tekli Varyant Analizi",
    tab_interval: "Genomik Aralık (Interval)",
    tab_scorers: "Atlas Scorer Kataloğu",

    // Tab 1: Upload
    upload_title: "VCF / TSV / CSV Sürükle & Bırak Toplu Analiz",
    upload_desc: "Varyant dosyalarınızı sürükleyin. 9 milyar önceden hesaplanmış AlphaGenome Atlas veritabanından saniyeler içinde zenginleştirilmiş skorları alın.",
    drop_title: "Dosyanızı buraya sürükleyip bırakın",
    drop_subtitle: "veya bilgisayarınızdan dosya seçmek için tıklayın",
    drop_quick_test: "Hızlı Test:",
    btn_sample_tsv: "Örnek TSV Yükle (3 varyant)",
    btn_sample_vcf: "Örnek VCF Yükle (FAS / APOL4)",
    file_preview_title: "Dosya Başlangıç Önizlemesi:",
    label_engine: "Analiz Motoru",
    label_tissue_opt: "Doku / Ontoloji (İsteğe bağlı)",
    label_gene_opt: "Gen Filtresi (İsteğe bağlı)",
    btn_start_batch: "Analizi Başlat",

    // Tab 2: Single Variant
    var_title: "Tekli Varyant Etki Analizi",
    var_desc: "HG38 koordinatlarında bir tek nükleotid varyantı girin. Atlas üzerinden anında AVI ve düzenleyici etki skorlarını görüntüleyin.",
    var_input_label: "Varyant (chr:pos:ref>alt)",
    var_engine_label: "Arama Modu",
    engine_atlas: "⚡ AlphaGenome Atlas (Anında 9B SNV)",
    engine_model: "🧬 Model Doğrudan İnferans",
    tissue_label: "Doku / Hücre Tipi",
    tissue_placeholder: "ör. colon, liver, brain",
    gene_label: "Hedef Gen (Odak)",
    gene_placeholder: "ör. APOL4, FAS, LDLR",
    btn_query_scores: "Skorları Sorgula",
    input_clear_hint: "Tıklayınca temizlenir",

    // Tab 3: Interval
    int_title: "Genomik Bölge / Aralık Atlas Sorgusu",
    int_desc: "Belirli bir genomik penceredeki tüm olası tek nükleotid varyantlarının Atlas skorlarını çekin ve en yüksek etkiye sahip varyantları keşfedin.",
    int_input_label: "Genomik Aralık (chr:start-end)",
    int_hint: "İpucu: Performans için ilk etapta 1-10 kb pencereler önerilir.",
    btn_query_interval: "Bölgedeki Tüm Varyantları Çek",

    // Tab 4: Scorers
    scorers_title: "AlphaGenome Atlas Scorer Kataloğu",
    scorers_desc: "Atlas veritabanında önceden hesaplanmış tüm düzenleyici etki scorer'larını inceleyin.",
    scorers_search: "Scorer veya doku ara...",

    // KPIs
    kpi_total: "Toplam Skor Satırı",
    kpi_variants: "Benzersiz Varyant",
    kpi_genes: "Etkilenen Genler",
    kpi_max_score: "Maksimum |Ham Skor|",

    // Interpretation
    interpret_title: "Biyolojik & Klinik Karar Destek Özeti",
    interpret_target: "Analiz Hedefi:",
    interpret_high_badge: "BELİRGİN DÜZENLEYİCİ ETKİ",
    interpret_med_badge: "ORTA SEVİYE DÜZENLEYİCİ ETKİ",
    interpret_low_badge: "ZAYIF ETKİ",
    interpret_negligible_badge: "İHMAL EDİLEBİLİR DÜZENLEYİCİ ETKİ",
    interpret_notes_title: "Klinik ve Biyoenformatik Karar Destek Notları:",
    note_1: "Skor büyüklüğü görecelidir; klinik sınır değeri değildir. Aynı gendeki patojenik/benign varyantlarla karşılaştırılmalıdır.",
    note_2: "Splice bölgesi varyantları için SpliceAI ve RNA-seq transkript modelleri ile çapraz doğrulama önerilir.",
    note_3: "Koordinatların GRCh38/hg38 referansı olduğunu doğrulayın.",

    // Table & Filters
    table_search_placeholder: "Gen, Scorer veya Doku ara...",
    filter_all: "Tüm Etki Seviyeleri",
    filter_high: "🔴 Belirgin Etki (Yüksek)",
    filter_med: "🟠 Orta Etki",
    filter_low: "🟡 Zayıf Etki",
    btn_export_tsv: "TSV İndir",
    btn_export_csv: "CSV İndir",
    col_variant: "Varyant ⬍",
    col_gene: "Gen ⬍",
    col_scorer: "Scorer ⬍",
    col_track: "Track / Doku ⬍",
    col_raw: "Ham Skor ⬍",
    col_quantile: "Yüzdelik (AVI) ⬍",
    col_impact: "Etki Seviyesi",
    pagination_showing: "sonuç gösteriliyor",
    btn_prev: "Önceki",
    btn_next: "Sonraki",

    // Official Atlas Portal
    badge_variant_text: "Varyant",
    atlas_official_api: "AlphaGenome Atlas (DeepMind)",
    atlas_avi_title: "AVI score",
    atlas_phred_unit: "Phred Score",
    atlas_fi_heading: "Feature importance",
    atlas_col_biosample: "Doku / Hücre Tipi",
    atlas_col_delta: "Δ Skor",
    atlas_pos_label: "Pozitif",
    atlas_neg_label: "Negatif",
    atlas_feature_breakdown: "Özellik Dağılımı",
    atlas_track_added_filter: "Tablo bu dokuya göre filtrelendi",

    // Modals
    settings_heading: "Sistem ve API Ayarları",
    label_api_key: "AlphaGenome API Anahtarı",
    api_key_hint: "Anahtar sunucu oturumunda geçerli kalır veya ALPHAGENOME_API_KEY ortam değişkeninden okunur.",
    label_cache_dir: "Yerel Önbellek Dizini",
    btn_save_settings: "Kaydet ve Bağlantıyı Güncelle",
    btn_close: "Kapat",
    btn_apply: "Uygula",

    // Errors & Autofix
    error_title: "Sorgu Uyarısı",
    autofix_heading: "Önerilen Düzeltme (hg38 Uyumlu):",
    btn_autofix_apply: "Doğru Varyantı Uygula ve Yeniden Sorgula",
    error_no_file: "Lütfen önce bir varyant dosyası yükleyin.",
    error_no_results: "Kriterlere uygun skor sonucu bulunamadı.",
    error_no_export: "Dışa aktarılacak veri yok.",
    error_api_key_empty: "Lütfen geçerli bir API anahtarı girin.",
    toast_key_saved: "API anahtarı başarıyla güncellendi.",
  },
  en: {
    brand_subtitle: "Google DeepMind Regulatory Genomics Analysis Studio",
    status_online: "System Ready",
    status_offline: "Offline",
    status_connecting: "Connecting...",
    settings_title: "Settings & API Key",

    // Tabs
    tab_upload: "Drag & Drop & Batch Files (VCF/TSV)",
    tab_variant: "Single Variant Analysis",
    tab_interval: "Genomic Interval",
    tab_scorers: "Atlas Scorer Catalog",

    // Tab 1: Upload
    upload_title: "VCF / TSV / CSV Drag & Drop Batch Analysis",
    upload_desc: "Drag and drop your variant files. Retrieve enriched regulatory effect scores in seconds from the 9-billion precomputed AlphaGenome Atlas database.",
    drop_title: "Drag & drop your file here",
    drop_subtitle: "or click to select a file from your computer",
    drop_quick_test: "Quick Test:",
    btn_sample_tsv: "Load Sample TSV (3 variants)",
    btn_sample_vcf: "Load Sample VCF (FAS / APOL4)",
    file_preview_title: "File Preview (First 5 rows):",
    label_engine: "Analysis Engine",
    label_tissue_opt: "Tissue / Ontology (Optional)",
    label_gene_opt: "Gene Filter (Optional)",
    btn_start_batch: "Run Batch Analysis",

    // Tab 2: Single Variant
    var_title: "Single Variant Regulatory Impact Analysis",
    var_desc: "Enter a single nucleotide variant in HG38 coordinates. Instantly view precomputed AVI and regulatory effect scores across tissues.",
    var_input_label: "Variant (chr:pos:ref>alt)",
    var_engine_label: "Search Engine",
    engine_atlas: "⚡ AlphaGenome Atlas (Instant 9B SNVs)",
    engine_model: "🧬 Direct Model Inference",
    tissue_label: "Tissue / Cell Type",
    tissue_placeholder: "e.g. colon, liver, brain",
    gene_label: "Target Gene (Focus)",
    gene_placeholder: "e.g. APOL4, FAS, LDLR",
    btn_query_scores: "Query Scores",
    input_clear_hint: "Click to clear",

    // Tab 3: Interval
    int_title: "Genomic Interval Atlas Query",
    int_desc: "Query all possible single nucleotide variants within a genomic window and discover high-impact regulatory disruptions.",
    int_input_label: "Genomic Interval (chr:start-end)",
    int_hint: "Tip: 1-10 kb windows are recommended for optimal performance.",
    btn_query_interval: "Fetch All Variants in Region",

    // Tab 4: Scorers
    scorers_title: "AlphaGenome Atlas Scorer Catalog",
    scorers_desc: "Browse all precomputed regulatory functional scorers available in the Atlas database.",
    scorers_search: "Search scorers or tissues...",

    // KPIs
    kpi_total: "Total Score Records",
    kpi_variants: "Unique Variants",
    kpi_genes: "Affected Genes",
    kpi_max_score: "Max |Raw Score|",

    // Interpretation
    interpret_title: "Biological & Clinical Decision Support Summary",
    interpret_target: "Target:",
    interpret_high_badge: "DISTINCT REGULATORY EFFECT",
    interpret_med_badge: "MODERATE REGULATORY EFFECT",
    interpret_low_badge: "WEAK EFFECT",
    interpret_negligible_badge: "NEGLIGIBLE REGULATORY EFFECT",
    interpret_notes_title: "Clinical and Bioinformatics Decision Notes:",
    note_1: "Score magnitude is relative; it does not represent an absolute clinical diagnostic threshold. Compare with known pathogenic/benign variants.",
    note_2: "Cross-validation with SpliceAI and RNA-seq transcript coverage is strongly advised for splice-altering candidates.",
    note_3: "Confirm that genomic coordinates adhere to the GRCh38/hg38 reference assembly.",

    // Table & Filters
    table_search_placeholder: "Search gene, scorer, tissue...",
    filter_all: "All Impact Tiers",
    filter_high: "🔴 Distinct Effect (High)",
    filter_med: "🟠 Moderate Effect",
    filter_low: "🟡 Weak Effect",
    btn_export_tsv: "Download TSV",
    btn_export_csv: "Download CSV",
    col_variant: "Variant ⬍",
    col_gene: "Gene ⬍",
    col_scorer: "Scorer ⬍",
    col_track: "Track / Tissue ⬍",
    col_raw: "Raw Score ⬍",
    col_quantile: "Quantile (AVI) ⬍",
    col_impact: "Impact Tier",
    pagination_showing: "results displayed",
    btn_prev: "Previous",
    btn_next: "Next",

    // Official Atlas Portal
    badge_variant_text: "Variant",
    atlas_official_api: "AlphaGenome Atlas (DeepMind)",
    atlas_avi_title: "AVI score",
    atlas_phred_unit: "Phred Score",
    atlas_fi_heading: "Feature importance",
    atlas_col_biosample: "Biosample / Tissue",
    atlas_col_delta: "Δ Score",
    atlas_pos_label: "Positive",
    atlas_neg_label: "Negative",
    atlas_feature_breakdown: "Feature breakdown",
    atlas_track_added_filter: "Table filtered to this biosample",

    // Modals
    settings_heading: "System & API Settings",
    label_api_key: "AlphaGenome API Key",
    api_key_hint: "Stored for the current server session or read from ALPHAGENOME_API_KEY environment variable.",
    label_cache_dir: "Local Disk Cache Directory",
    btn_save_settings: "Save and Reconnect",
    btn_close: "Close",
    btn_apply: "Apply",

    // Errors & Autofix
    error_title: "Query Warning",
    autofix_heading: "Suggested Correction (hg38 Reference):",
    btn_autofix_apply: "Apply Correct Variant & Re-query",
    error_no_file: "Please upload a variant file first.",
    error_no_results: "No score records found matching criteria.",
    error_no_export: "No data available to export.",
    error_api_key_empty: "Please enter a valid API key.",
    toast_key_saved: "API key successfully updated.",
  }
};

// ─── State ─────────────────────────────────────────────────────────────────
let currentLang = localStorage.getItem("alphagenome_lang") || "tr";
let currentFile = null;
let currentFileData = null;
let currentResults = [];
let filteredResults = [];
let sortCol = "quantile_score";
let sortAsc = false;
let currentPage = 1;
const pageSize = 25;
let timerInterval = null;
let startTime = 0;
let lastQueryParams = null;

// ─── DOM Initialization ───────────────────────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
  initLanguage();
  initTabs();
  initDragAndDrop();
  initForms();
  initSearchInputBehavior();
  initTableActions();
  initSettingsModal();
  initErrorModal();
  checkSystemStatus();
  loadScorersList();
});

/* ─── Language / Localization (i18n) ─────────────────────────────────────── */
function initLanguage() {
  const btnTr = document.getElementById("langTrBtn");
  const btnEn = document.getElementById("langEnBtn");

  if (btnTr && btnEn) {
    btnTr.addEventListener("click", () => applyLanguage("tr"));
    btnEn.addEventListener("click", () => applyLanguage("en"));
  }

  applyLanguage(currentLang);
}

function t(key) {
  return (i18n[currentLang] && i18n[currentLang][key]) || (i18n.tr[key] || key);
}

function applyLanguage(lang) {
  currentLang = lang;
  localStorage.setItem("alphagenome_lang", lang);
  document.documentElement.lang = lang;

  // Update switcher buttons
  const btnTr = document.getElementById("langTrBtn");
  const btnEn = document.getElementById("langEnBtn");
  if (btnTr) btnTr.classList.toggle("active", lang === "tr");
  if (btnEn) btnEn.classList.toggle("active", lang === "en");

  // Update data-i18n texts
  document.querySelectorAll("[data-i18n]").forEach(el => {
    const key = el.dataset.i18n;
    if (i18n[lang] && i18n[lang][key]) {
      el.textContent = i18n[lang][key];
    }
  });

  // Update data-i18n-placeholder
  document.querySelectorAll("[data-i18n-placeholder]").forEach(el => {
    const key = el.dataset.i18nPlaceholder;
    if (i18n[lang] && i18n[lang][key]) {
      el.placeholder = i18n[lang][key];
    }
  });

  // Update data-i18n-title
  document.querySelectorAll("[data-i18n-title]").forEach(el => {
    const key = el.dataset.i18nTitle;
    if (i18n[lang] && i18n[lang][key]) {
      el.title = i18n[lang][key];
    }
  });

  // Update table headers
  updateTableHeadersLanguage();

  // If results are visible, re-render interpretation and table
  if (currentResults && currentResults.length) {
    renderResults(currentResults, lastQueryParams ? (lastQueryParams.variant || lastQueryParams.interval || lastQueryParams.fileName) : "", lastQueryParams ? lastQueryParams.gene : "");
  }
}

function updateTableHeadersLanguage() {
  const mapping = [
    { sort: "variant", key: "col_variant" },
    { sort: "gene_name", key: "col_gene" },
    { sort: "variant_scorer", key: "col_scorer" },
    { sort: "track_name", key: "col_track" },
    { sort: "raw_score", key: "col_raw" },
    { sort: "quantile_score", key: "col_quantile" },
  ];

  mapping.forEach(m => {
    const th = document.querySelector(`#resultsTable th[data-sort="${m.sort}"]`);
    if (th) th.textContent = t(m.key);
  });
}

/* ─── Search Bar Click-to-Clear Behavior ─────────────────────────────────── */
function initSearchInputBehavior() {
  const variantInput = document.getElementById("variantInput");
  if (variantInput) {
    // "varyant arama barı tıklayınca boşalsın"
    variantInput.addEventListener("click", function () {
      this.value = "";
    });
  }

  const intervalInput = document.getElementById("intervalInput");
  if (intervalInput) {
    intervalInput.addEventListener("click", function () {
      this.value = "";
    });
  }
}

/* ─── Navigation Tabs ─────────────────────────────────────────────────── */
function initTabs() {
  const tabs = document.querySelectorAll(".nav-tab");
  tabs.forEach(tab => {
    tab.addEventListener("click", () => {
      tabs.forEach(t => t.classList.remove("active"));
      document.querySelectorAll(".tab-pane").forEach(p => p.classList.remove("active"));

      tab.classList.add("active");
      const targetId = tab.dataset.tab;
      const pane = document.getElementById(targetId);
      if (pane) pane.classList.add("active");
    });
  });
}

/* ─── Drag & Drop File Upload ─────────────────────────────────────────── */
function initDragAndDrop() {
  const dropZone = document.getElementById("dropZone");
  const fileInput = document.getElementById("fileInput");
  const removeFileBtn = document.getElementById("removeFileBtn");
  const loadSampleTsvBtn = document.getElementById("loadSampleTsvBtn");
  const loadSampleVcfBtn = document.getElementById("loadSampleVcfBtn");
  const startBatchBtn = document.getElementById("startBatchBtn");

  dropZone.addEventListener("click", (e) => {
    if (e.target.tagName.toLowerCase() === "button") return;
    fileInput.click();
  });

  dropZone.addEventListener("dragover", (e) => {
    e.preventDefault();
    dropZone.classList.add("dragover");
  });

  dropZone.addEventListener("dragleave", () => {
    dropZone.classList.remove("dragover");
  });

  dropZone.addEventListener("drop", (e) => {
    e.preventDefault();
    dropZone.classList.remove("dragover");
    if (e.dataTransfer.files.length) {
      handleFile(e.dataTransfer.files[0]);
    }
  });

  fileInput.addEventListener("change", () => {
    if (fileInput.files.length) {
      handleFile(fileInput.files[0]);
    }
  });

  removeFileBtn.addEventListener("click", () => {
    currentFile = null;
    currentFileData = null;
    fileInput.value = "";
    document.getElementById("filePreviewCard").classList.add("hidden");
    dropZone.classList.remove("hidden");
  });

  // Sample data loaders (using verified hg38 reference coordinates)
  loadSampleTsvBtn.addEventListener("click", (e) => {
    e.stopPropagation();
    const sampleTsv = `chr\tpos\tref\talt\nchr22\t36201698\tA\tC\nchr10\t89010000\tG\tA\nchr19\t11100000\tT\tC\n`;
    const blob = new Blob([sampleTsv], { type: "text/tab-separated-values" });
    const file = new File([blob], "ornek_varyantlar.tsv", { type: "text/tab-separated-values" });
    handleFile(file);
  });

  loadSampleVcfBtn.addEventListener("click", (e) => {
    e.stopPropagation();
    const sampleVcf = `##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr22\t36201698\trs123\tA\tC\t.\tPASS\t.\nchr10\t89010000\trs456\tG\tA\t.\tPASS\t.\nchr19\t11100000\trs789\tT\tC\t.\tPASS\t.\n`;
    const blob = new Blob([sampleVcf], { type: "text/plain" });
    const file = new File([blob], "ornek_klinik.vcf", { type: "text/plain" });
    handleFile(file);
  });

  startBatchBtn.addEventListener("click", executeBatchAnalysis);
}

function handleFile(file) {
  currentFile = file;
  const reader = new FileReader();

  reader.onload = (e) => {
    currentFileData = e.target.result;
    parseAndPreviewFile(file.name, currentFileData, file.size);
  };

  reader.readAsText(file);
}

function parseAndPreviewFile(fileName, content, fileSize) {
  const lines = content.split(/\r?\n/).filter(l => l.trim().length > 0);
  const ext = fileName.split(".").pop().toUpperCase();

  document.getElementById("fileExtBadge").textContent = ext;
  document.getElementById("fileName").textContent = fileName;

  const sizeKb = (fileSize / 1024).toFixed(1);
  const dataRows = lines.filter(l => !l.startsWith("#"));
  const varCount = Math.max(0, dataRows.length - (lines[0].includes("chr") || lines[0].includes("pos") ? 1 : 0));

  document.getElementById("fileDetails").textContent = `${sizeKb} KB · ${varCount} ${currentLang === "tr" ? "varyant tespit edildi" : "variants detected"}`;

  // Generate 5-row preview
  const headRow = document.getElementById("previewTableHead");
  const bodyRows = document.getElementById("previewTableBody");
  headRow.innerHTML = "";
  bodyRows.innerHTML = "";

  const previewLines = dataRows.slice(0, 6);
  if (!previewLines.length) return;

  const delimiter = fileName.toLowerCase().endsWith(".csv") ? "," : "\t";
  const headers = previewLines[0].split(delimiter);

  headers.forEach(h => {
    const th = document.createElement("th");
    th.textContent = h.trim();
    headRow.appendChild(th);
  });

  previewLines.slice(1).forEach(line => {
    const tr = document.createElement("tr");
    const cols = line.split(delimiter);
    cols.forEach(c => {
      const td = document.createElement("td");
      td.textContent = c.trim();
      tr.appendChild(td);
    });
    bodyRows.appendChild(tr);
  });

  document.getElementById("dropZone").classList.add("hidden");
  document.getElementById("filePreviewCard").classList.remove("hidden");
}

/* ─── Query Forms ────────────────────────────────────────────────────────── */
function initForms() {
  // Single Variant Form
  const singleForm = document.getElementById("singleVariantForm");
  singleForm.addEventListener("submit", (e) => {
    e.preventDefault();
    const variant = document.getElementById("variantInput").value.trim();
    const engine = document.getElementById("variantEngine").value;
    const tissue = document.getElementById("variantTissue").value.trim();
    const gene = document.getElementById("variantGene").value.trim();

    executeSingleVariantQuery({ variant, engine, tissue, gene });
  });

  // Quick Examples Pills
  document.querySelectorAll(".btn-chip").forEach(btn => {
    btn.addEventListener("click", () => {
      document.getElementById("variantInput").value = btn.dataset.var;
      document.getElementById("variantGene").value = btn.dataset.gene || "";
    });
  });

  // Interval Form
  const intervalForm = document.getElementById("intervalForm");
  intervalForm.addEventListener("submit", (e) => {
    e.preventDefault();
    const interval = document.getElementById("intervalInput").value.trim();
    const tissue = document.getElementById("intervalTissue").value.trim();
    const gene = document.getElementById("intervalGene").value.trim();

    executeIntervalQuery({ interval, tissue, gene });
  });
}

async function executeSingleVariantQuery(params) {
  lastQueryParams = params;
  showLoading(
    currentLang === "tr" ? "AlphaGenome Atlas Sorgulanıyor..." : "Querying AlphaGenome Atlas...",
    `${currentLang === "tr" ? "Varyant" : "Variant"}: ${params.variant} (${params.engine === "atlas" ? "Atlas 9B SNV" : "Model Direct"})`
  );

  try {
    const res = await fetch("/api/query_variant", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(params),
    });

    const data = await res.json();
    hideLoading();

    if (!res.ok) {
      showErrorModal(data);
      return;
    }

    renderResults(data.scores || [], params.variant, params.gene, data.avi_summary);
  } catch (err) {
    hideLoading();
    showErrorModal({ error: err.message });
  }
}

async function executeIntervalQuery(params) {
  lastQueryParams = params;
  showLoading(
    currentLang === "tr" ? "Genomik Aralık Taranıyor..." : "Scanning Genomic Interval...",
    `${currentLang === "tr" ? "Bölge" : "Region"}: ${params.interval}`
  );

  try {
    const res = await fetch("/api/query_interval", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(params),
    });

    const data = await res.json();
    hideLoading();

    if (!res.ok) {
      showErrorModal(data);
      return;
    }

    renderResults(data.scores || [], params.interval, params.gene, data.avi_summary);
  } catch (err) {
    hideLoading();
    showErrorModal({ error: err.message });
  }
}

async function executeBatchAnalysis() {
  if (!currentFile || !currentFileData) {
    showToast(t("error_no_file"), "error");
    return;
  }

  const engine = document.getElementById("batchEngine").value;
  const tissue = document.getElementById("batchTissue").value.trim();
  const gene = document.getElementById("batchGene").value.trim();

  lastQueryParams = { fileName: currentFile.name, engine, tissue, gene };
  showLoading(
    currentLang === "tr" ? "Toplu Varyant Analizi Yürütülüyor..." : "Running Batch Variant Analysis...",
    `${currentLang === "tr" ? "Dosya" : "File"}: ${currentFile.name}`
  );

  try {
    const res = await fetch("/api/query_batch", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        fileName: currentFile.name,
        fileContent: currentFileData,
        engine,
        tissue,
        gene,
      }),
    });

    const data = await res.json();
    hideLoading();

    if (!res.ok) {
      showErrorModal(data);
      return;
    }

    renderResults(data.scores || [], currentFile.name, gene, data.avi_summary);
  } catch (err) {
    hideLoading();
    showErrorModal({ error: err.message });
  }
}

/* ─── Results Dashboard Rendering ────────────────────────────────────────── */
function renderResults(scores, queryLabel, focusGene, aviSummary) {
  currentResults = scores || [];
  filteredResults = [...currentResults];
  currentPage = 1;

  const section = document.getElementById("resultsSection");
  if (!currentResults.length) {
    showToast(t("error_no_results"), "info");
    section.classList.add("hidden");
    return;
  }

  section.classList.remove("hidden");
  section.scrollIntoView({ behavior: "smooth" });

  // Official Atlas AVI Portal Card
  const atlasPortalCard = document.getElementById("atlasPortalCard");
  if (aviSummary && aviSummary.has_avi) {
    renderAtlasPortalCard(aviSummary, queryLabel, focusGene);
    atlasPortalCard.classList.remove("hidden");
  } else {
    atlasPortalCard.classList.add("hidden");
  }

  // KPIs
  const totalRows = currentResults.length;
  const uniqueVariants = new Set(currentResults.map(r => r.variant)).size;
  const uniqueGenes = new Set(currentResults.map(r => r.gene_name).filter(Boolean)).size;
  const maxRaw = currentResults.reduce((max, r) => Math.max(max, Math.abs(r.raw_score || 0)), 0);

  document.getElementById("kpiTotal").textContent = totalRows.toLocaleString();
  document.getElementById("kpiVariants").textContent = uniqueVariants.toLocaleString();
  document.getElementById("kpiGenes").textContent = uniqueGenes.toLocaleString();
  document.getElementById("kpiMaxScore").textContent = maxRaw.toFixed(4);

  // Biological Interpretation
  renderInterpretation(currentResults, queryLabel, focusGene);

  // Table
  applyTableFilterAndSort();
}

function renderAtlasPortalCard(summary, queryLabel, focusGene) {
  document.getElementById("atlasPortalVariant").textContent = queryLabel;
  document.getElementById("atlasPortalGeneTag").textContent = focusGene || "Locus";

  // Phred Score
  const phredEl = document.getElementById("atlasPhredVal");
  if (summary.phred_score !== null && summary.phred_score !== undefined) {
    phredEl.textContent = Number(summary.phred_score).toFixed(1);
  } else {
    phredEl.textContent = "-";
  }

  // Top Percentile Badge
  const percentileText = document.getElementById("atlasPercentileText");
  const topPct = summary.top_percentile !== null && summary.top_percentile !== undefined
    ? Number(summary.top_percentile).toFixed(1)
    : "-";
  percentileText.textContent = currentLang === "tr"
    ? `Tüm genomdaki varyantlar arasında ilk %${topPct} etki diliminde`
    : `Top ${topPct}% predicted impact of all genome-wide SNVs`;

  // Feature Importance List
  const fiList = document.getElementById("atlasFiList");
  fiList.innerHTML = "";

  const features = summary.features || [];
  if (!features.length) {
    fiList.innerHTML = `<div class="text-xs text-dim">Özellik bilgisi bulunamadı</div>`;
    return;
  }

  const positiveFeatures = features.filter(f => f.direction === "positive");
  const negativeFeatures = features.filter(f => f.direction === "negative");

  function renderGroup(groupFeatures, label, icon) {
    if (!groupFeatures.length) return;
    const groupHeading = document.createElement("div");
    groupHeading.className = "atlas-fi-group-label";
    groupHeading.innerHTML = `<span>${icon}</span> <span>${label}</span>`;
    fiList.appendChild(groupHeading);

    groupFeatures.forEach(feat => {
      const item = document.createElement("div");
      item.className = "atlas-fi-item";
      item.dataset.featureId = feat.id;
      item.innerHTML = `
        <div class="atlas-fi-left">
          <span class="atlas-fi-dot" style="background: ${feat.color};"></span>
          <span class="atlas-fi-name">${feat.name}</span>
        </div>
        <div class="atlas-fi-right">
          <span class="atlas-fi-pct">${feat.percentage.toFixed(2)}%</span>
          <span class="atlas-fi-arrow">›</span>
        </div>
      `;
      item.addEventListener("click", () => {
        selectAtlasFeature(feat, queryLabel);
      });
      fiList.appendChild(item);
    });
  }

  renderGroup(positiveFeatures, t("atlas_pos_label"), "▲");
  renderGroup(negativeFeatures, t("atlas_neg_label"), "▼");

  // Select default feature: prioritize one with tracks, or the highest feature (DNASE-seq if present)
  let defaultFeature = features.find(f => f.name.includes("DNASE")) ||
                       features.find(f => f.tracks && f.tracks.length > 0) ||
                       features[0];
  selectAtlasFeature(defaultFeature, queryLabel);
}

function selectAtlasFeature(feat, queryLabel) {
  // Active state on items
  document.querySelectorAll(".atlas-fi-item").forEach(el => {
    if (el.dataset.featureId === feat.id) {
      el.classList.add("active");
    } else {
      el.classList.remove("active");
    }
  });

  // Header
  const dot = document.getElementById("atlasFeatureDot");
  dot.style.background = feat.color;
  dot.style.boxShadow = `0 0 12px ${feat.color}`;

  const title = document.getElementById("atlasBreakdownTitle");
  title.textContent = `${feat.name} (${feat.percentage.toFixed(2)}%) ${t("atlas_feature_breakdown") || "Feature breakdown"}`;

  const varName = document.getElementById("atlasBreakdownVarName");
  varName.textContent = queryLabel;

  // Tracks
  const container = document.getElementById("atlasTracksItems");
  container.innerHTML = "";

  const tracks = feat.tracks || [];
  if (!tracks.length) {
    const isConservation = feat.id.includes("CACTUS") || feat.id.includes("PHASTCONS");
    const noteText = isConservation
      ? (currentLang === "tr"
          ? `Bu özellik (${feat.name}), dokuya özel epigenomik izlerden ziyade tüm genom genelindeki evrimsel dizi korumasını (conservation) temsil eder. Ham model katkı skoru: ${feat.raw_val}`
          : `This feature (${feat.name}) represents genome-wide multi-species evolutionary conservation rather than tissue-specific assays. Raw model contribution: ${feat.raw_val}`)
      : (currentLang === "tr"
          ? `Bu özellik için detaylı doku izi bulunamadı. Ham skor: ${feat.raw_val}`
          : `No tissue-specific tracks available for this feature. Raw score: ${feat.raw_val}`);

    container.innerHTML = `
      <div class="p-4 text-center text-dim text-sm font-mono" style="background: rgba(255,255,255,0.02); border-radius: 6px;">
        ${noteText}
      </div>
    `;
    return;
  }

  // Calculate dynamic axis range
  let maxAbs = tracks.reduce((m, t) => Math.max(m, Math.abs(t.score || 0)), 0);
  if (maxAbs < 0.1) maxAbs = 0.1;
  const step = maxAbs;

  document.getElementById("scaleMinTick").textContent = `-${step.toFixed(2)}`;
  document.getElementById("scaleMidNegTick").textContent = `-${(step / 2).toFixed(2)}`;
  document.getElementById("scaleMidPosTick").textContent = `+${(step / 2).toFixed(2)}`;
  document.getElementById("scaleMaxTick").textContent = `+${step.toFixed(2)}`;

  tracks.forEach(tr => {
    const row = document.createElement("div");
    row.className = "atlas-track-row";

    const score = tr.score || 0;
    const isPos = score >= 0;
    const widthPct = Math.min(50, (Math.abs(score) / step) * 50);

    row.innerHTML = `
      <div class="atlas-track-bio" title="${tr.biosample} (${tr.name || ''})">
        <span class="atlas-track-name">${tr.biosample}</span>
        ${tr.name && tr.name !== tr.biosample ? `<span class="atlas-track-curie">${tr.name}</span>` : ""}
      </div>
      <div class="atlas-diverging-bar-container">
        <div class="atlas-axis-center"></div>
        <div class="atlas-diverging-bar ${isPos ? "positive" : "negative"}" style="width: ${widthPct.toFixed(1)}%;"></div>
      </div>
      <div class="atlas-track-val-cell">
        <span class="atlas-track-score ${isPos ? "pos" : "neg"}">${score > 0 ? `+${score.toFixed(2)}` : score.toFixed(2)}</span>
        <button type="button" class="atlas-btn-add-track" title="Tabloda bu dokuyu filtrele" data-bname="${tr.biosample}">+</button>
      </div>
    `;

    const addBtn = row.querySelector(".atlas-btn-add-track");
    addBtn.addEventListener("click", () => {
      const filterInput = document.getElementById("tableFilterInput");
      filterInput.value = tr.biosample;
      applyTableFilterAndSort();
      showToast(`${t("atlas_track_added_filter")}: ${tr.biosample}`, "info");
      document.getElementById("resultsTable").scrollIntoView({ behavior: "smooth" });
    });

    container.appendChild(row);
  });
}

function renderInterpretation(scores, queryLabel, focusGene) {
  const interpretVariantText = document.getElementById("interpretVariantText");
  const interpretBody = document.getElementById("interpretBody");
  interpretVariantText.textContent = `${t("interpret_target")} ${queryLabel}${focusGene ? ` · ${focusGene}` : ""}`;

  let maxQuantile = 0;
  let maxScoreRow = null;
  for (const r of scores) {
    const q = r.quantile_score || 0;
    if (q > maxQuantile) {
      maxQuantile = q;
      maxScoreRow = r;
    }
  }

  let badgeClass = "badge-impact negligible";
  let badgeText = t("interpret_negligible_badge");
  let desc = currentLang === "tr"
    ? "Bu varyantın incelenen dokularda gen ekspresyonu veya kromatin erişilebilirliğinde anlamlı bir sapmaya yol açmadığı tahmin edilmektedir."
    : "This variant is predicted to have no significant perturbation on gene expression or chromatin accessibility across evaluated tissues.";

  if (maxQuantile >= 0.98 || (maxScoreRow && Math.abs(maxScoreRow.raw_score) >= 1.5)) {
    badgeClass = "badge-impact high";
    badgeText = t("interpret_high_badge");
    desc = currentLang === "tr"
      ? `AlphaGenome Atlas bu varyantı düzenleyici varyantlar arasında en üst %2'lik dilimde (AVI: ${(maxQuantile * 100).toFixed(1)}%) derecelendirdi. Hedef genin ekspresyonunda veya bağlanma profilinde kayda değer değişiklik beklenmektedir.`
      : `AlphaGenome Atlas ranks this variant in the top 2% of regulatory effects (AVI: ${(maxQuantile * 100).toFixed(1)}%). Substantial alterations in target gene expression or binding affinity are predicted.`;
  } else if (maxQuantile >= 0.90 || (maxScoreRow && Math.abs(maxScoreRow.raw_score) >= 0.8)) {
    badgeClass = "badge-impact med";
    badgeText = t("interpret_med_badge");
    desc = currentLang === "tr"
      ? `Bu varyant hedef gendeki transkripsiyon veya kromatin düzenlemesinde orta düzeyde etki yaratmaktadır (Yüzdelik: ${(maxQuantile * 100).toFixed(1)}%).`
      : `This variant exhibits a moderate regulatory impact on target gene transcription or chromatin accessibility (Quantile: ${(maxQuantile * 100).toFixed(1)}%).`;
  } else if (maxQuantile >= 0.70) {
    badgeClass = "badge-impact low";
    badgeText = t("interpret_low_badge");
    desc = currentLang === "tr"
      ? `Belirtilen varyant zayıf bir etki göstermektedir. Biyolojik önemi diğer klinik ve moleküler kanıtlarla desteklenmelidir.`
      : `This variant displays a borderline or minor effect. Biological relevance should be corroborated with clinical context.`;
  }

  interpretBody.innerHTML = `
    <div class="mb-2">
      <span class="${badgeClass}">${badgeText}</span>
      ${maxScoreRow ? `<span class="text-dim text-xs font-mono ml-2">${currentLang === "tr" ? "En Yüksek Scorer" : "Top Scorer"}: ${maxScoreRow.variant_scorer} (${maxScoreRow.track_name || "Global"})</span>` : ""}
    </div>
    <p class="mb-3 text-secondary">${desc}</p>
    <div class="interpret-notes">
      <strong>${t("interpret_notes_title")}</strong>
      <ul class="mt-1 text-xs" style="padding-left: 20px;">
        <li>${t("note_1")}</li>
        <li>${t("note_2")}</li>
        <li>${t("note_3")}</li>
      </ul>
    </div>
  `;
}

/* ─── Table Actions, Filters & Pagination ────────────────────────────────── */
function initTableActions() {
  const filterInput = document.getElementById("tableFilterInput");
  const impactFilter = document.getElementById("impactFilterSelect");
  const prevBtn = document.getElementById("prevPageBtn");
  const nextBtn = document.getElementById("nextPageBtn");

  filterInput.addEventListener("input", applyTableFilterAndSort);
  impactFilter.addEventListener("change", applyTableFilterAndSort);

  document.querySelectorAll("#resultsTable th[data-sort]").forEach(th => {
    th.addEventListener("click", () => {
      const col = th.dataset.sort;
      if (sortCol === col) {
        sortAsc = !sortAsc;
      } else {
        sortCol = col;
        sortAsc = false;
      }
      applyTableFilterAndSort();
    });
  });

  prevBtn.addEventListener("click", () => {
    if (currentPage > 1) {
      currentPage--;
      renderTablePage();
    }
  });

  nextBtn.addEventListener("click", () => {
    const totalPages = Math.ceil(filteredResults.length / pageSize) || 1;
    if (currentPage < totalPages) {
      currentPage++;
      renderTablePage();
    }
  });

  document.getElementById("exportTsvBtn").addEventListener("click", () => exportData("tsv"));
  document.getElementById("exportCsvBtn").addEventListener("click", () => exportData("csv"));
}

function applyTableFilterAndSort() {
  const query = document.getElementById("tableFilterInput").value.toLowerCase().trim();
  const impact = document.getElementById("impactFilterSelect").value;

  filteredResults = currentResults.filter(r => {
    const textMatch = !query ||
      (r.variant && r.variant.toLowerCase().includes(query)) ||
      (r.gene_name && r.gene_name.toLowerCase().includes(query)) ||
      (r.variant_scorer && r.variant_scorer.toLowerCase().includes(query)) ||
      (r.biosample_name && r.biosample_name.toLowerCase().includes(query)) ||
      (r.track_name && r.track_name.toLowerCase().includes(query)) ||
      (r.ontology_curie && r.ontology_curie.toLowerCase().includes(query));

    if (!textMatch) return false;

    const q = r.quantile_score || 0;
    if (impact === "high") return q >= 0.95 || Math.abs(r.raw_score || 0) >= 1.0;
    if (impact === "med") return q >= 0.80 && q < 0.95;
    if (impact === "low") return q < 0.80;

    return true;
  });

  filteredResults.sort((a, b) => {
    let valA = a[sortCol];
    let valB = b[sortCol];

    if (typeof valA === "number" || typeof valB === "number") {
      valA = Math.abs(valA || 0);
      valB = Math.abs(valB || 0);
    } else {
      valA = String(valA || "").toLowerCase();
      valB = String(valB || "").toLowerCase();
    }

    if (valA < valB) return sortAsc ? -1 : 1;
    if (valA > valB) return sortAsc ? 1 : -1;
    return 0;
  });

  currentPage = 1;
  renderTablePage();
}

function renderTablePage() {
  const tbody = document.getElementById("resultsTableBody");
  tbody.innerHTML = "";

  const totalPages = Math.ceil(filteredResults.length / pageSize) || 1;
  const startIdx = (currentPage - 1) * pageSize;
  const pageRows = filteredResults.slice(startIdx, startIdx + pageSize);

  pageRows.forEach(row => {
    const tr = document.createElement("tr");

    const q = row.quantile_score || 0;
    let tierBadge = `<span class="badge-tier low">${t("interpret_low_badge")}</span>`;
    if (q >= 0.95 || Math.abs(row.raw_score || 0) >= 1.0) {
      tierBadge = `<span class="badge-tier high">${t("interpret_high_badge")}</span>`;
    } else if (q >= 0.80) {
      tierBadge = `<span class="badge-tier med">${t("interpret_med_badge")}</span>`;
    }

    const mainLabel = row.biosample_name || row.track_name || "-";
    const subLabel = (row.biosample_name && row.track_name && row.track_name !== row.biosample_name)
      ? `<span class="text-dim text-xs block font-mono">${row.track_name}</span>`
      : "";

    tr.innerHTML = `
      <td class="font-mono font-bold">${row.variant || "-"}</td>
      <td>
        <span class="gene-tag">${row.gene_name || "-"}</span>
        ${row.gene_id ? `<span class="text-dim text-xs font-mono block">${row.gene_id}</span>` : ""}
      </td>
      <td>
        <span class="scorer-tag">${row.variant_scorer || "-"}</span>
      </td>
      <td>
        <span class="text-sm font-semibold">${mainLabel}</span>
        ${subLabel}
        ${row.ontology_curie ? `<span class="curie-tag font-mono block">${row.ontology_curie}</span>` : ""}
      </td>
      <td class="font-mono">${row.raw_score !== null && row.raw_score !== undefined ? Number(row.raw_score).toFixed(4) : "-"}</td>
      <td class="font-mono font-bold ${q >= 0.95 ? "text-rose" : q >= 0.80 ? "text-cyan" : ""}">
        ${row.quantile_score !== null && row.quantile_score !== undefined ? Number(row.quantile_score).toFixed(4) : "-"}
      </td>
      <td>${tierBadge}</td>
    `;
    tbody.appendChild(tr);
  });

  // Pagination UI
  const total = filteredResults.length;
  document.getElementById("paginationInfo").textContent = `${total.toLocaleString()} ${t("pagination_showing")}`;
  document.getElementById("pageIndicator").textContent = `${currentPage} / ${totalPages}`;
  document.getElementById("prevPageBtn").disabled = currentPage <= 1;
  document.getElementById("nextPageBtn").disabled = currentPage >= totalPages;
}

function exportData(format) {
  if (!currentResults.length) {
    showToast(t("error_no_export"), "info");
    return;
  }

  const cols = ["variant", "gene_name", "gene_id", "variant_scorer", "track_name", "ontology_curie", "raw_score", "quantile_score"];
  const delimiter = format === "csv" ? "," : "\t";

  const headerRow = cols.join(delimiter);
  const rows = currentResults.map(r => {
    return cols.map(c => {
      let val = r[c] === null || r[c] === undefined ? "" : String(r[c]);
      if (format === "csv" && val.includes(",")) val = `"${val}"`;
      return val;
    }).join(delimiter);
  });

  const content = [headerRow, ...rows].join("\n");
  const blob = new Blob([content], { type: format === "csv" ? "text/csv;charset=utf-8;" : "text/tab-separated-values;charset=utf-8;" });
  const url = URL.createObjectURL(blob);

  const a = document.createElement("a");
  a.href = url;
  a.download = `alphagenome_scores_${Date.now()}.${format}`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

/* ─── Loading Overlay ─────────────────────────────────────────────────────── */
function showLoading(title, desc) {
  document.getElementById("loadingTitle").textContent = title;
  document.getElementById("loadingDesc").textContent = desc;
  document.getElementById("loadingModal").classList.remove("hidden");

  startTime = Date.now();
  const timer = document.getElementById("loadingTimer");
  timer.textContent = `${currentLang === "tr" ? "Geçen Süre" : "Elapsed"}: 0.0s`;

  if (timerInterval) clearInterval(timerInterval);
  timerInterval = setInterval(() => {
    const elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
    timer.textContent = `${currentLang === "tr" ? "Geçen Süre" : "Elapsed"}: ${elapsed}s`;
  }, 100);
}

function hideLoading() {
  if (timerInterval) clearInterval(timerInterval);
  document.getElementById("loadingModal").classList.add("hidden");
}

/* ─── Modern Error Dialog & Auto-Fix ─────────────────────────────────────── */
function initErrorModal() {
  const closeBtn = document.getElementById("closeErrorModalBtn");
  const dismissBtn = document.getElementById("dismissErrorBtn");
  const applyFixBtn = document.getElementById("applyAutoFixBtn");

  if (closeBtn) closeBtn.addEventListener("click", closeErrorModal);
  if (dismissBtn) dismissBtn.addEventListener("click", closeErrorModal);

  if (applyFixBtn) {
    applyFixBtn.addEventListener("click", () => {
      const fixedVar = applyFixBtn.dataset.suggestedVariant;
      if (fixedVar) {
        document.getElementById("variantInput").value = fixedVar;
        closeErrorModal();
        showToast(
          currentLang === "tr" ? `Varyant ${fixedVar} olarak düzeltildi. Analiz başlatılıyor...` : `Variant updated to ${fixedVar}. Starting query...`,
          "success"
        );
        const engine = document.getElementById("variantEngine").value;
        const tissue = document.getElementById("variantTissue").value.trim();
        const gene = document.getElementById("variantGene").value.trim();
        executeSingleVariantQuery({ variant: fixedVar, engine, tissue, gene });
      }
    });
  }
}

function showErrorModal(data) {
  const modal = document.getElementById("errorModal");
  if (!modal) return;

  const titleEl = document.getElementById("errorModalTitle");
  const subtitleEl = document.getElementById("errorModalSubtitle");
  const messageEl = document.getElementById("errorModalMessage");
  const autoFixContainer = document.getElementById("autoFixContainer");
  const autoFixVarText = document.getElementById("autoFixVarText");
  const applyAutoFixBtn = document.getElementById("applyAutoFixBtn");

  titleEl.textContent = t("error_title");
  
  const rawMsg = (currentLang === "en" && data.error_en) ? data.error_en : (data.error || "Bilinmeyen bir hata oluştu.");
  subtitleEl.textContent = data.error_type === "ref_mismatch" ? (currentLang === "tr" ? "Referans Baz Uyumsuzluğu" : "Reference Base Mismatch") : "";
  messageEl.textContent = rawMsg;

  // Reference base mismatch auto-fix suggestion
  if (data.error_type === "ref_mismatch" && data.suggested_variant) {
    autoFixContainer.classList.remove("hidden");
    autoFixVarText.textContent = data.suggested_variant;
    applyAutoFixBtn.dataset.suggestedVariant = data.suggested_variant;
  } else {
    autoFixContainer.classList.add("hidden");
  }

  modal.classList.remove("hidden");
}

function closeErrorModal() {
  const modal = document.getElementById("errorModal");
  if (modal) modal.classList.add("hidden");
}

/* ─── Toast Notifications ─────────────────────────────────────────────────── */
function showToast(message, type = "info") {
  const container = document.getElementById("toastContainer");
  if (!container) return;

  const toast = document.createElement("div");
  toast.className = `toast toast-${type}`;
  toast.innerHTML = `
    <div class="toast-content">${message}</div>
  `;

  container.appendChild(toast);
  setTimeout(() => toast.classList.add("show"), 10);

  setTimeout(() => {
    toast.classList.remove("show");
    setTimeout(() => toast.remove(), 300);
  }, 4000);
}

/* ─── Settings Modal ─────────────────────────────────────────────────────── */
function initSettingsModal() {
  const modal = document.getElementById("settingsModal");
  const openBtn = document.getElementById("quickSettingsBtn");
  const closeBtn = document.getElementById("closeSettingsBtn");
  const saveBtn = document.getElementById("saveApiKeyBtn");
  const apiKeyInput = document.getElementById("apiKeyInput");

  openBtn.addEventListener("click", () => {
    modal.classList.remove("hidden");
    apiKeyInput.focus();
  });

  closeBtn.addEventListener("click", () => modal.classList.add("hidden"));

  saveBtn.addEventListener("click", async () => {
    const key = apiKeyInput.value.trim();
    if (!key) {
      showToast(t("error_api_key_empty"), "error");
      return;
    }

    try {
      const res = await fetch("/api/set_key", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ api_key: key }),
      });
      const data = await res.json();
      if (res.ok) {
        showToast(t("toast_key_saved"), "success");
        modal.classList.add("hidden");
        checkSystemStatus();
      } else {
        showToast("Hata: " + (data.error || "Kaydedilemedi"), "error");
      }
    } catch (e) {
      showToast("Hata: " + e.message, "error");
    }
  });
}

/* ─── System Status & Scorer Catalog ─────────────────────────────────────── */
async function checkSystemStatus() {
  const statusEl = document.getElementById("connectionStatus");
  const textEl = document.getElementById("statusText");

  try {
    const res = await fetch("/api/status");
    const data = await res.json();

    if (res.ok) {
      statusEl.className = "status-indicator online";
      textEl.textContent = t("status_online");
      if (data.cache_dir) {
        const cacheInput = document.getElementById("cacheDirInput");
        if (cacheInput) cacheInput.value = data.cache_dir;
      }
    } else {
      statusEl.className = "status-indicator offline";
      textEl.textContent = t("status_offline");
    }
  } catch {
    statusEl.className = "status-indicator offline";
    textEl.textContent = t("status_offline");
  }
}

async function loadScorersList() {
  const container = document.getElementById("scorersList");
  try {
    const res = await fetch("/api/scorers");
    const data = await res.json();

    if (!res.ok || !data.scorers || !data.scorers.length) {
      container.innerHTML = `<div class="text-dim p-4">Scorer listesi alınamadı.</div>`;
      return;
    }

    window.allScorers = data.scorers;
    renderScorers(data.scorers);

    // Scorer search filter
    const searchInput = document.getElementById("scorerSearchInput");
    if (searchInput) {
      searchInput.addEventListener("input", (e) => {
        const q = e.target.value.toLowerCase().trim();
        const filtered = window.allScorers.filter(s =>
          s.name.toLowerCase().includes(q) || (s.description && s.description.toLowerCase().includes(q))
        );
        renderScorers(filtered);
      });
    }
  } catch (err) {
    container.innerHTML = `<div class="text-dim p-4">Scorer bağlantı hatası: ${err.message}</div>`;
  }
}

function renderScorers(scorers) {
  const container = document.getElementById("scorersList");
  container.innerHTML = "";

  scorers.forEach(s => {
    const card = document.createElement("div");
    card.className = "scorer-card";
    card.innerHTML = `
      <div class="scorer-card-header">
        <div class="scorer-name">${s.name}</div>
        <span class="pill ${s.is_signed ? "pill-cyan" : "pill-purple"}">${s.is_signed ? (currentLang === "tr" ? "Yönlü (Signed)" : "Signed") : (currentLang === "tr" ? "Yönsüz" : "Unsigned")}</span>
      </div>
      <div class="scorer-desc">${s.description || (currentLang === "tr" ? "Fonksiyonel varyant etki skoru modeli." : "Functional variant effect scorer model.")}</div>
      <div class="scorer-meta">
        <span>${currentLang === "tr" ? "İz Sayısı" : "Track Count"}: <strong>${s.track_count || "Çoklu"}</strong></span>
      </div>
    `;
    container.appendChild(card);
  });
}
