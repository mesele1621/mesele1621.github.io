/* =========================================================================
   Advanced Primer Designer Pro — script.js
   Modular, no dependencies. Organized as small namespaces so later modules
   (qPCR, Gibson, mutagenesis, etc.) can be added without touching this file.
   ========================================================================= */
'use strict';

/* ---------------------------------------------------------------------
   0. APP STATE
   --------------------------------------------------------------------- */
const AppState = {
  sequence: '',      // cleaned, uppercase, DNA-alphabet sequence currently loaded
  header: '',
  isRNA: false,
  results: []         // last generated primer pairs
};

/* ---------------------------------------------------------------------
   1. SEQUENCE PARSING  (Module 1)
   --------------------------------------------------------------------- */
const SeqParser = {
  /** Detects format and returns { header, sequence, isRNA } */
  parse(raw, forcedFormat) {
    const text = raw.trim();
    if (!text) throw new Error('No sequence provided.');

    let format = forcedFormat && forcedFormat !== 'auto' ? forcedFormat : this.detectFormat(text);
    let header = '', body = '';

    if (format === 'genbank') {
      const originMatch = text.match(/ORIGIN([\s\S]*?)(?:\/\/|$)/i);
      const defMatch = text.match(/DEFINITION\s+(.*)/i);
      header = defMatch ? defMatch[1].trim() : 'GenBank record';
      body = originMatch ? originMatch[1].replace(/[0-9\s]/g, '') : '';
      if (!body) throw new Error('Could not find an ORIGIN block in this GenBank record.');
    } else if (format === 'fasta') {
      const lines = text.split(/\r?\n/);
      if (lines[0].startsWith('>')) {
        header = lines[0].slice(1).trim();
        body = lines.slice(1).join('');
      } else {
        header = 'Untitled sequence';
        body = lines.join('');
      }
    } else {
      // raw DNA or RNA
      header = 'Untitled sequence';
      body = text.replace(/\r?\n/g, '');
    }

    body = body.replace(/[^A-Za-z]/g, '').toUpperCase();
    if (!body.length) throw new Error('No valid bases found after parsing.');

    const isRNA = format === 'rna' || (body.includes('U') && !body.includes('T'));
    const sequence = isRNA ? body.replace(/U/g, 'T') : body; // internally always work in DNA space

    return { header, sequence, isRNA };
  },

  detectFormat(text) {
    if (/^LOCUS\s/im.test(text) || /ORIGIN/i.test(text)) return 'genbank';
    if (text.startsWith('>')) return 'fasta';
    if (/U/i.test(text) && !/T/i.test(text.replace(/^>.*$/m, ''))) return 'rna';
    return 'raw';
  }
};

/* ---------------------------------------------------------------------
   2. SEQUENCE STATISTICS  (Module 1 display)
   --------------------------------------------------------------------- */
const SeqStats = {
  compute(seq) {
    const len = seq.length;
    const counts = { A: 0, T: 0, G: 0, C: 0, N: 0 };
    for (const ch of seq) {
      if (counts[ch] !== undefined) counts[ch]++;
      else counts.N++;
    }
    const gc = counts.G + counts.C;
    const at = counts.A + counts.T;
    return {
      length: len,
      counts,
      gcPercent: len ? (gc / len * 100) : 0,
      atPercent: len ? (at / len * 100) : 0,
      nCount: counts.N
    };
  }
};

/* ---------------------------------------------------------------------
   3. CORE ALGORITHMS: reverse complement, Tm, GC
   --------------------------------------------------------------------- */
const Bio = {
  COMPLEMENT: { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' },

  reverseComplement(seq) {
    let out = '';
    for (let i = seq.length - 1; i >= 0; i--) {
      out += this.COMPLEMENT[seq[i]] || 'N';
    }
    return out;
  },

  gcPercent(seq) {
    if (!seq.length) return 0;
    let gc = 0;
    for (const ch of seq) if (ch === 'G' || ch === 'C') gc++;
    return (gc / seq.length) * 100;
  },

  /** Wallace rule — quick estimate, used only as a sanity fallback for <14 nt. */
  tmWallace(seq) {
    let a = 0, t = 0, g = 0, c = 0;
    for (const ch of seq) {
      if (ch === 'A') a++; else if (ch === 'T') t++;
      else if (ch === 'G') g++; else if (ch === 'C') c++;
    }
    return 2 * (a + t) + 4 * (g + c);
  },

  // Unified SantaLucia (1998) nearest-neighbor ΔH (kcal/mol) and ΔS (cal/mol·K)
  NN_PARAMS: {
    AA: { h: -7.9, s: -22.2 }, TT: { h: -7.9, s: -22.2 },
    AT: { h: -7.2, s: -20.4 }, TA: { h: -7.2, s: -21.3 },
    CA: { h: -8.5, s: -22.7 }, TG: { h: -8.5, s: -22.7 },
    GT: { h: -8.4, s: -22.4 }, AC: { h: -8.4, s: -22.4 },
    CT: { h: -7.8, s: -21.0 }, AG: { h: -7.8, s: -21.0 },
    GA: { h: -8.2, s: -22.2 }, TC: { h: -8.2, s: -22.2 },
    CG: { h: -10.6, s: -27.2 }, GC: { h: -9.8, s: -24.4 },
    GG: { h: -8.0, s: -19.9 }, CC: { h: -8.0, s: -19.9 }
  },
  // Initiation terms (unified SantaLucia 1998)
  INIT_TERMINAL_GC: { h: 0.1, s: -2.8 },
  INIT_TERMINAL_AT: { h: 2.3, s: 4.1 },

  /**
   * Nearest-neighbor Tm with Na+ salt correction (Owczarzy-style log correction).
   * primerConc in M (default 250nM), saltConc in mM.
   */
  tmNearestNeighbor(seq, primerConcM = 2.5e-7, saltConcMM = 50) {
    if (seq.length < 2) return this.tmWallace(seq);
    if (seq.length < 14) {
      // Nearest-neighbor is noisy on very short primers; blend with Wallace.
      return this.tmWallace(seq);
    }

    let dH = 0, dS = 0;
    for (let i = 0; i < seq.length - 1; i++) {
      const pair = seq[i] + seq[i + 1];
      const p = this.NN_PARAMS[pair];
      if (!p) continue; // skip ambiguous bases
      dH += p.h;
      dS += p.s;
    }
    // terminal correction based on first/last base identity
    const endBases = [seq[0], seq[seq.length - 1]];
    for (const b of endBases) {
      const term = (b === 'G' || b === 'C') ? this.INIT_TERMINAL_GC : this.INIT_TERMINAL_AT;
      dH += term.h;
      dS += term.s;
    }

    const R = 1.987; // cal/(mol*K)
    // Tm (Kelvin) for a non-self-complementary duplex: divide primerConc by 4
    const Tm_K = (dH * 1000) / (dS + R * Math.log(primerConcM / 4)) ;
    let Tm_C = Tm_K - 273.15;

    // Salt correction (Owczarzy 2004 simplified log form)
    const saltM = saltConcMM / 1000;
    Tm_C = Tm_C + 16.6 * Math.log10(saltM / (1 + 0.7 * saltM)) - 16.6 * Math.log10(0.05);

    return Tm_C;
  },

  /** Longest run of a single repeated base, e.g. AAAAA -> 5 */
  longestHomopolymer(seq) {
    let longest = 1, run = 1;
    for (let i = 1; i < seq.length; i++) {
      if (seq[i] === seq[i - 1]) { run++; longest = Math.max(longest, run); }
      else run = 1;
    }
    return longest;
  },

  /** Longest run of a repeated dinucleotide, e.g. ATATAT -> 3 repeats */
  longestDinucRepeat(seq) {
    let longest = 1;
    for (let start = 0; start < 2; start++) {
      let run = 1;
      for (let i = start; i + 4 <= seq.length; i += 2) {
        if (seq.slice(i, i + 2) === seq.slice(i + 2, i + 4)) run++;
        else { longest = Math.max(longest, run); run = 1; }
      }
      longest = Math.max(longest, run);
    }
    return longest;
  },

  gcClampCount(seq, window = 5) {
    const tail = seq.slice(-window);
    let n = 0;
    for (const ch of tail) if (ch === 'G' || ch === 'C') n++;
    return n;
  },

  /**
   * Very lightweight self-complementarity scan used for hairpin / self-dimer /
   * cross-dimer screening. Returns the length of the longest complementary
   * run found between seqA and reverse-complement of seqB (self-dimer when
   * seqB===seqA). This is a heuristic stand-in for full thermodynamic folding
   * (which will be added as its own module later) but is enough to reject
   * obviously bad primers.
   */
  longestComplementaryRun(seqA, seqB) {
    const revB = this.reverseComplement(seqB);
    let best = 0;
    for (let i = 0; i < seqA.length; i++) {
      for (let j = 0; j < revB.length; j++) {
        let run = 0;
        while (i + run < seqA.length && j + run < revB.length && seqA[i + run] === revB[j + run]) run++;
        if (run > best) best = run;
      }
    }
    return best;
  },

  /** Hairpin heuristic: does the primer fold back on itself (min loop 3nt)? */
  hairpinRisk(seq) {
    let best = 0;
    const n = seq.length;
    for (let loop = 3; loop < n - 4; loop++) {
      for (let stem = 4; stem <= (n - loop) / 2; stem++) {
        const left = seq.slice(0, stem);
        const rightStart = stem + loop;
        if (rightStart + stem > n) continue;
        const right = seq.slice(rightStart, rightStart + stem);
        const rc = this.reverseComplement(right);
        let match = 0;
        for (let k = 0; k < stem; k++) if (left[stem - 1 - k] === rc[stem - 1 - k]) match++;
        if (match > best) best = match;
      }
    }
    return best; // longer matched stem = higher risk
  }
};

/* ---------------------------------------------------------------------
   4. PRIMER DESIGN ENGINE  (Module 2 + Module 3 conventional PCR)
   --------------------------------------------------------------------- */
const PrimerEngine = {
  /**
   * Generates candidate primers along `seq` in the given direction ('fwd' scans
   * left-to-right taking the window as-is; 'rev' scans and returns the
   * reverse complement of the window so it reads 5'->3' as a reverse primer).
   */
  scanCandidates(seq, params, direction, regionStart, regionEnd) {
    const candidates = [];
    for (let len = params.lenMin; len <= params.lenMax; len++) {
      for (let start = regionStart; start + len <= regionEnd; start++) {
        const windowDNA = seq.slice(start, start + len);
        if (/[^ATGC]/.test(windowDNA)) continue; // skip Ns / ambiguous for now

        const primerSeq = direction === 'fwd' ? windowDNA : Bio.reverseComplement(windowDNA);
        const gc = Bio.gcPercent(primerSeq);
        if (gc < params.gcMin || gc > params.gcMax) continue;

        const tm = Bio.tmNearestNeighbor(primerSeq, 2.5e-7, params.saltConc);
        if (tm < params.tmMin || tm > params.tmMax) continue;

        const homopolymer = Bio.longestHomopolymer(primerSeq);
        if (homopolymer > params.maxHomopolymer) continue;

        const repeats = Bio.longestDinucRepeat(primerSeq);
        if (repeats > params.maxRepeats) continue;

        const clampCount = Bio.gcClampCount(primerSeq);
        if (clampCount < params.gcClampMin) continue;

        if (params.rejectHairpin && Bio.hairpinRisk(primerSeq) >= Math.min(6, len - 2)) continue;
        if (params.rejectSelfDimer && Bio.longestComplementaryRun(primerSeq, primerSeq) >= Math.min(6, len - 2)) continue;

        candidates.push({
          sequence: primerSeq,
          start,               // 0-based, in the ORIGINAL sequence coordinate
          end: start + len,    // exclusive
          length: len,
          gc, tm,
          direction,
          score: this.scoreCandidate({ gc, tm, length: len }, params)
        });
      }
    }
    // best candidates first
    candidates.sort((a, b) => b.score - a.score);
    return candidates;
  },

  scoreCandidate(c, params) {
    // 0-100: penalize distance from optimal Tm/GC/length, each weighted.
    const tmPenalty = Math.abs(c.tm - params.tmOpt) * 3;
    const gcPenalty = Math.abs(c.gc - params.gcOpt) * 1.2;
    const lenPenalty = Math.abs(c.length - params.lenOpt) * 1.5;
    let score = 100 - tmPenalty - gcPenalty - lenPenalty;
    return Math.max(0, Math.min(100, Math.round(score)));
  },

  /** Pairs forward + reverse candidates into valid product-size pairs, ranked. */
  buildPairs(fwdCandidates, revCandidates, params, maxPairs = 25) {
    const pairs = [];
    // limit search space for performance: take best N of each side
    const fwdTop = fwdCandidates.slice(0, 60);
    const revTop = revCandidates.slice(0, 60);

    for (const f of fwdTop) {
      for (const r of revTop) {
        if (r.start <= f.start) continue; // reverse primer must sit downstream
        const productSize = r.end - f.start;
        if (productSize < params.prodMin || productSize > params.prodMax) continue;

        if (params.rejectCrossDimer &&
            Bio.longestComplementaryRun(f.sequence, r.sequence) >= Math.min(6, Math.min(f.length, r.length) - 2)) {
          continue;
        }

        const tmDiff = Math.abs(f.tm - r.tm);
        const pairScore = Math.round((f.score + r.score) / 2 - tmDiff * 2);
        pairs.push({ forward: f, reverse: r, productSize, tmDiff, pairScore: Math.max(0, Math.min(100, pairScore)) });
      }
    }
    pairs.sort((a, b) => b.pairScore - a.pairScore);
    return pairs.slice(0, maxPairs);
  },

  design(seq, params) {
    let regionStart = 0, regionEnd = seq.length;
    if (params.targetRegion) {
      const m = params.targetRegion.match(/(\d+)\s*-\s*(\d+)/);
      if (m) {
        regionStart = Math.max(0, parseInt(m[1], 10) - 1);
        regionEnd = Math.min(seq.length, parseInt(m[2], 10));
      }
    }
    // forward primers are picked from the first half of the region, reverse from the second half,
    // but we allow overlap and let product-size filtering do the real work.
    const fwdCandidates = this.scanCandidates(seq, params, 'fwd', regionStart, regionEnd);
    const revCandidates = this.scanCandidates(seq, params, 'rev', regionStart, regionEnd);
    return this.buildPairs(fwdCandidates, revCandidates, params);
  }
};

/* ---------------------------------------------------------------------
   5. UI: navigation + theme
   --------------------------------------------------------------------- */
const UI = {
  init() {
    document.querySelectorAll('.nav-item[data-page]').forEach(btn => {
      btn.addEventListener('click', () => this.goTo(btn.dataset.page));
    });
    document.getElementById('themeBtn').addEventListener('click', () => this.toggleTheme());
    document.getElementById('settingsThemeCheck').addEventListener('change', e => {
      this.setTheme(e.target.checked ? 'light' : 'dark');
    });
  },

  goTo(page) {
    document.querySelectorAll('.page').forEach(p => p.style.display = 'none');
    document.getElementById('page-' + page).style.display = 'block';
    document.querySelectorAll('.nav-item[data-page]').forEach(b => b.classList.toggle('active', b.dataset.page === page));
  },

  toggleTheme() {
    const current = document.documentElement.getAttribute('data-theme');
    this.setTheme(current === 'dark' ? 'light' : 'dark');
  },

  setTheme(theme) {
    document.documentElement.setAttribute('data-theme', theme);
    document.getElementById('themeBtn').textContent = theme === 'dark' ? '🌙' : '☀️';
    document.getElementById('themeLabel').textContent = theme === 'dark' ? 'Dark mode' : 'Light mode';
    document.getElementById('settingsThemeCheck').checked = theme === 'light';
  },

  toast(msg) {
    const t = document.createElement('div');
    t.className = 'toast';
    t.textContent = msg;
    document.body.appendChild(t);
    setTimeout(() => t.remove(), 2600);
  }
};

/* ---------------------------------------------------------------------
   6. MODULE 1 CONTROLLER: input, drag/drop, stats rendering
   --------------------------------------------------------------------- */
const InputController = {
  selectedFormat: 'auto',

  init() {
    document.querySelectorAll('#formatChips .chip').forEach(chip => {
      chip.addEventListener('click', () => {
        document.querySelectorAll('#formatChips .chip').forEach(c => c.classList.remove('active'));
        chip.classList.add('active');
        this.selectedFormat = chip.dataset.format;
      });
    });

    const dz = document.getElementById('dropzone');
    const fileInput = document.getElementById('fileInput');
    dz.addEventListener('click', () => fileInput.click());
    dz.addEventListener('dragover', e => { e.preventDefault(); dz.classList.add('drag-over'); });
    dz.addEventListener('dragleave', () => dz.classList.remove('drag-over'));
    dz.addEventListener('drop', e => {
      e.preventDefault();
      dz.classList.remove('drag-over');
      if (e.dataTransfer.files.length) this.handleFile(e.dataTransfer.files[0]);
    });
    fileInput.addEventListener('change', e => {
      if (e.target.files.length) this.handleFile(e.target.files[0]);
    });

    document.getElementById('parseBtn').addEventListener('click', () => this.parseFromTextarea());
    document.getElementById('clearBtn').addEventListener('click', () => {
      document.getElementById('seqInput').value = '';
      document.getElementById('statsCard').style.display = 'none';
    });
    document.getElementById('sampleBtn').addEventListener('click', () => {
      document.getElementById('seqInput').value =
        '>example_Soymovirus_partial_ORF\n' +
        'ATGGCTAGCAAAGGTGAAGAACTGTTCACTGGCGTGGTGCCTATTCTGGTGGAACTGGATGGTGATGTG' +
        'AACGGTCATAAATTTAGCGTGAGCGGTGAAGGTGAAGGTGATGCTACCTATGGTAAACTGACCCTGAAA' +
        'TTCATTTGCACCACCGGTAAACTGCCTGTGCCTTGGCCTACCCTGGTGACCACCCTGACCTATGGTGTG' +
        'CAGTGCTTTAGCCGTTATCCTGATCATATGAAACGTCATGATTTCTTCAAAAGCGCTATGCCT';
      this.parseFromTextarea();
    });

    document.getElementById('revCompBtn').addEventListener('click', () => {
      if (!AppState.sequence) return;
      const rc = Bio.reverseComplement(AppState.sequence);
      navigator.clipboard?.writeText(rc).catch(() => {});
      UI.toast('Reverse complement copied to clipboard (' + rc.length + ' bp).');
    });
    document.getElementById('goDesignBtn').addEventListener('click', () => UI.goTo('design'));
  },

  handleFile(file) {
    const reader = new FileReader();
    reader.onload = e => {
      document.getElementById('seqInput').value = e.target.result;
      let fmt = this.selectedFormat;
      if (fmt === 'auto') {
        const name = file.name.toLowerCase();
        if (name.endsWith('.gb') || name.endsWith('.gbk')) fmt = 'genbank';
        else if (name.endsWith('.fa') || name.endsWith('.fasta')) fmt = 'fasta';
      }
      this.parseFromTextarea(fmt);
    };
    reader.readAsText(file);
  },

  parseFromTextarea(forcedFormat) {
    const raw = document.getElementById('seqInput').value;
    try {
      const { header, sequence, isRNA } = SeqParser.parse(raw, forcedFormat || this.selectedFormat);
      AppState.sequence = sequence;
      AppState.header = header;
      AppState.isRNA = isRNA;
      this.renderStats(sequence);
      UI.toast('Parsed ' + sequence.length + ' bp' + (isRNA ? ' (RNA → converted to cDNA space)' : '') + '.');
    } catch (err) {
      UI.toast('Error: ' + err.message);
    }
  },

  renderStats(seq) {
    const stats = SeqStats.compute(seq);
    document.getElementById('statsCard').style.display = 'block';
    document.getElementById('statLength').textContent = stats.length.toLocaleString() + ' bp';
    document.getElementById('statGC').textContent = stats.gcPercent.toFixed(1) + '%';
    document.getElementById('statAT').textContent = stats.atPercent.toFixed(1) + '%';
    document.getElementById('statN').textContent = stats.nCount;

    const gcBar = document.getElementById('gcBar');
    gcBar.innerHTML = `<span class="gc" style="width:${stats.gcPercent}%"></span><span class="at" style="width:${stats.atPercent}%"></span>`;

    document.getElementById('baseA').textContent = stats.counts.A.toLocaleString();
    document.getElementById('baseT').textContent = stats.counts.T.toLocaleString();
    document.getElementById('baseG').textContent = stats.counts.G.toLocaleString();
    document.getElementById('baseC').textContent = stats.counts.C.toLocaleString();

    this.renderRuler(seq);
  },

  /** Renders a binned base-composition ruler so even long sequences stay readable. */
  renderRuler(seq) {
    const ruler = document.getElementById('baseRuler');
    ruler.innerHTML = '';
    const bins = 120;
    const binSize = Math.max(1, Math.ceil(seq.length / bins));
    for (let i = 0; i < seq.length; i += binSize) {
      const chunk = seq.slice(i, i + binSize);
      const counts = { A: 0, T: 0, G: 0, C: 0, N: 0 };
      for (const ch of chunk) counts[ch] !== undefined ? counts[ch]++ : counts.N++;
      const dominant = Object.entries(counts).sort((a, b) => b[1] - a[1])[0][0];
      const seg = document.createElement('i');
      seg.className = dominant;
      seg.title = `pos ${i + 1}-${Math.min(i + binSize, seq.length)}: A${counts.A} T${counts.T} G${counts.G} C${counts.C}`;
      ruler.appendChild(seg);
    }
  }
};

/* ---------------------------------------------------------------------
   7. MODULE 2/3 CONTROLLER: reading params + triggering design
   --------------------------------------------------------------------- */
const DesignController = {
  init() {
    document.getElementById('designBtn').addEventListener('click', () => this.run());
  },

  readParams() {
    const num = id => parseFloat(document.getElementById(id).value);
    return {
      lenMin: num('lenMin'), lenOpt: num('lenOpt'), lenMax: num('lenMax'),
      gcMin: num('gcMin'), gcOpt: num('gcOpt'), gcMax: num('gcMax'),
      tmMin: num('tmMin'), tmOpt: num('tmOpt'), tmMax: num('tmMax'),
      prodMin: num('prodMin'), prodMax: num('prodMax'),
      gcClampMin: num('gcClamp'),
      maxHomopolymer: num('maxHomopolymer'),
      maxRepeats: num('maxRepeats'),
      maxNs: num('maxNs'),
      saltConc: num('saltConc'),
      rejectHairpin: document.getElementById('rejectHairpin').checked,
      rejectSelfDimer: document.getElementById('rejectSelfDimer').checked,
      rejectCrossDimer: document.getElementById('rejectCrossDimer').checked,
      targetRegion: document.getElementById('targetRegion').value.trim()
    };
  },

  run() {
    if (!AppState.sequence) {
      UI.toast('Load a sequence in Module 1 first.');
      UI.goTo('input');
      return;
    }
    const status = document.getElementById('designStatus');
    status.textContent = 'Designing…';
    // yield to the browser so the status text paints before the (synchronous) search runs
    setTimeout(() => {
      const params = this.readParams();
      const t0 = performance.now();
      const pairs = PrimerEngine.design(AppState.sequence, params);
      const elapsed = (performance.now() - t0).toFixed(0);
      AppState.results = pairs;
      ResultsController.render(pairs);
      status.textContent = pairs.length
        ? `Found ${pairs.length} candidate pair(s) in ${elapsed} ms.`
        : `No pairs matched these constraints (${elapsed} ms) — try loosening Tm/GC/product-size ranges.`;
      UI.goTo('results');
    }, 30);
  }
};

/* ---------------------------------------------------------------------
   8. MODULE 17/18 CONTROLLER: results table + export
   --------------------------------------------------------------------- */
const ResultsController = {
  init() {
    document.getElementById('exportCsvBtn').addEventListener('click', () => this.exportCsv());
    document.getElementById('exportFastaBtn').addEventListener('click', () => this.exportFasta());
    document.getElementById('exportTxtBtn').addEventListener('click', () => this.exportTxt());
    document.getElementById('copyBtn').addEventListener('click', () => this.copyToClipboard());
  },

  scoreClass(score) {
    if (score >= 80) return 'green';
    if (score >= 60) return 'yellow';
    if (score >= 40) return 'orange';
    return 'red';
  },

  markClamp(seq) {
    // visually highlight the last 5 bases (GC clamp region)
    const head = seq.slice(0, -5);
    const tail = seq.slice(-5);
    return `${head}<span class="clamp">${tail}</span>`;
  },

  render(pairs) {
    const body = document.getElementById('resultsBody');
    if (!pairs.length) {
      body.innerHTML = `<tr><td colspan="8"><div class="empty-state"><div class="icon">🧬</div>No primer pairs matched — loosen your constraints and try again.</div></td></tr>`;
      return;
    }
    let rows = '';
    pairs.forEach((pair, idx) => {
      const pairNum = idx + 1;
      [['F', pair.forward], ['R', pair.reverse]].forEach(([label, p]) => {
        const cls = this.scoreClass(p.score);
        rows += `<tr>
          <td>#${pairNum}</td>
          <td>${label} · pos ${p.start + 1}-${p.end}</td>
          <td class="seq-cell">${this.markClamp(p.sequence)}</td>
          <td>${p.length}</td>
          <td>${p.gc.toFixed(1)}</td>
          <td>${p.tm.toFixed(1)}</td>
          <td>${label === 'F' ? pair.productSize : ''}</td>
          <td><span class="score-pill ${cls}">${p.score}</span></td>
        </tr>`;
      });
    });
    document.getElementById('resultsBody').innerHTML = rows;
  },

  toRows() {
    const rows = [['pair', 'primer', 'sequence_5to3', 'length', 'gc_percent', 'tm_c', 'position_start', 'position_end', 'product_size', 'score']];
    AppState.results.forEach((pair, idx) => {
      [['F', pair.forward], ['R', pair.reverse]].forEach(([label, p]) => {
        rows.push([idx + 1, label, p.sequence, p.length, p.gc.toFixed(1), p.tm.toFixed(1), p.start + 1, p.end, label === 'F' ? pair.productSize : '', p.score]);
      });
    });
    return rows;
  },

  downloadBlob(content, filename, type) {
    const blob = new Blob([content], { type });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    URL.revokeObjectURL(url);
  },

  exportCsv() {
    if (!AppState.results.length) return UI.toast('No results to export yet.');
    const csv = this.toRows().map(r => r.map(v => `"${v}"`).join(',')).join('\n');
    this.downloadBlob(csv, 'primer_results.csv', 'text/csv');
  },

  exportFasta() {
    if (!AppState.results.length) return UI.toast('No results to export yet.');
    let fasta = '';
    AppState.results.forEach((pair, idx) => {
      fasta += `>pair${idx + 1}_forward pos=${pair.forward.start + 1}-${pair.forward.end} Tm=${pair.forward.tm.toFixed(1)} GC=${pair.forward.gc.toFixed(1)}\n${pair.forward.sequence}\n`;
      fasta += `>pair${idx + 1}_reverse pos=${pair.reverse.start + 1}-${pair.reverse.end} Tm=${pair.reverse.tm.toFixed(1)} GC=${pair.reverse.gc.toFixed(1)}\n${pair.reverse.sequence}\n`;
    });
    this.downloadBlob(fasta, 'primer_results.fasta', 'text/plain');
  },

  exportTxt() {
    if (!AppState.results.length) return UI.toast('No results to export yet.');
    const txt = this.toRows().map(r => r.join('\t')).join('\n');
    this.downloadBlob(txt, 'primer_results.txt', 'text/plain');
  },

  copyToClipboard() {
    if (!AppState.results.length) return UI.toast('No results to export yet.');
    const txt = this.toRows().map(r => r.join('\t')).join('\n');
    navigator.clipboard?.writeText(txt).then(() => UI.toast('Results copied to clipboard.'))
      .catch(() => UI.toast('Clipboard copy failed — your browser may block it on this page.'));
  }
};

/* ---------------------------------------------------------------------
   9. BOOTSTRAP
   --------------------------------------------------------------------- */
document.addEventListener('DOMContentLoaded', () => {
  UI.init();
  InputController.init();
  DesignController.init();
  ResultsController.init();
  UI.setTheme('dark');
});
