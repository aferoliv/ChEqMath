const FRONT_VERSION = "0.3.1"; // mesma versão de APP_VERSION em backend/app/main.py
const $ = (id) => document.getElementById(id);
const OPENED_AS_FILE = location.protocol === "file:";
const ALL = "__all__";
const PROTON = "0002";
const WATER = "0001";
const WATER_CONC = 55.6; // mol L⁻¹, água pura (usada só para exibição)
const nuOf = {}; // código da espécie → [[componente, ν], ...]
let names = {};

function status(msg, isError = false) {
  const el = $("status");
  el.textContent = msg;
  el.className = isError ? "error" : "";
}

async function api(path, body) {
  const opts = body
    ? { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) }
    : {};
  const r = await fetch(path, opts);
  if (!r.ok) {
    const e = await r.json().catch(() => ({}));
    throw new Error(typeof e.detail === "string" ? e.detail : r.statusText);
  }
  return r.json();
}

// Números nas tabelas: 4 algarismos significativos. |x| ≥ 0,001 em decimal; abaixo, em notação científica.
const SUP = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };

function fmt(x) {
  if (x === null || x === undefined || Number.isNaN(x)) return "—";
  if (x === 0) return "0";
  if (Math.abs(x) >= 0.001) return x.toPrecision(4).replace(".", ",");
  const [m, e] = x.toExponential(3).split("e");
  const exp = [...String(parseInt(e, 10))].map((c) => SUP[c]).join("");
  return `${m.replace(".", ",")} × 10${exp}`;
}

function activityOpts() {
  if (!$("dh-on").checked) return { activity: "ideal" };
  const o = { activity: "dh_ext", ion_size: parseFloat($("dh-a").value) || 4.5 };
  const I = $("dh-I").value;
  if (I !== "") o.ionic_strength = parseFloat(I);
  return o;
}

function kinds() {
  return [...document.querySelectorAll("input[name=kind]:checked")].map((x) => x.value);
}

// Hidroxila: a concentração informada entra como T_H⁺ negativo (OH⁻ = −H⁺ no balanço de prótons).
const OH = "__OH__";

function totals() {
  const out = {};
  document.querySelectorAll("#comp-list .row").forEach((row) => {
    const num = row.querySelector("input[type=number]");
    if (num.value === "") return;
    const v = parseFloat(num.value);
    if (row.dataset.code === OH) out[PROTON] = (out[PROTON] || 0) - v;
    else out[row.dataset.code] = (out[row.dataset.code] || 0) + v;
  });
  return out;
}

// ---------- categorias (lista suspensa com seleção múltipla) ----------

function selectedCategories() {
  const boxes = [...document.querySelectorAll("#cat-box input[type=checkbox]")];
  const all = boxes.find((b) => b.value === ALL);
  if (!boxes.some((b) => b.checked) || (all && all.checked)) return [];
  return boxes.filter((b) => b.checked && b.value !== ALL).map((b) => b.value);
}

function updateCategorySummary() {
  const sel = selectedCategories();
  $("cat-summary").textContent = sel.length ? sel.join(", ") : "Todas as categorias";
}

function buildCategoryBox(categories) {
  const box = $("cat-box");
  box.textContent = "";
  [ALL, ...categories].forEach((value) => {
    const label = document.createElement("label");
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.value = value;
    cb.checked = value === ALL;
    label.append(cb, document.createTextNode(value === ALL ? " Todas as categorias" : " " + value));
    box.append(label);
    cb.addEventListener("change", () => onCategoryChange(cb));
  });
}

function onCategoryChange(changed) {
  const boxes = [...document.querySelectorAll("#cat-box input[type=checkbox]")];
  const all = boxes.find((b) => b.value === ALL);
  if (changed.value === ALL && changed.checked) {
    boxes.forEach((b) => { if (b !== all) b.checked = false; });
  } else if (changed.checked && all) {
    all.checked = false;
  }
  if (!boxes.some((b) => b.checked) && all) all.checked = true;
  updateCategorySummary();
  loadComponents().catch((e) => status(e.message, true));
}

// ---------- componentes (lista suspensa filtrada) ----------

async function loadComponents() {
  const params = selectedCategories().map((c) => `categories=${encodeURIComponent(c)}`).join("&");
  const items = (await api(`/api/components${params ? "?" + params : ""}`)).items;
  const sel = $("comp");
  sel.textContent = "";
  const ph = document.createElement("option");
  ph.value = "";
  ph.textContent = "— escolha um componente —";
  sel.append(ph);
  const add = (value, text) => {
    const o = document.createElement("option");
    o.value = value;
    o.textContent = text;
    sel.append(o);
  };
  items.forEach((c) => { names[c.code] = c.name; });
  // H⁺ e OH⁻ sempre no início; demais componentes em ordem alfabética (vinda da API)
  const proton = items.find((c) => c.code === PROTON);
  if (proton) add(PROTON, proton.name);
  add(OH, "Hidroxila (OH⁻)");
  items.filter((c) => !c.is_solvent && c.code !== PROTON).forEach((c) => add(c.code, c.name));
}

function addSelected() {
  const sel = $("comp");
  const code = sel.value;
  if (!code) return status("Escolha um componente na lista.", true);
  if (document.querySelector(`#comp-list .row[data-code="${code}"]`)) {
    return status("Este componente já está na lista de totais.", true);
  }
  const row = document.createElement("div");
  row.className = "row";
  row.dataset.code = code;
  const label = document.createElement("span");
  label.textContent = code === OH ? "Hidroxila (OH⁻)" : names[code] || code;
  const num = document.createElement("input");
  num.type = "number";
  num.step = "any";
  num.placeholder = "mol/L";
  num.addEventListener("input", refreshRefSelect);
  const rm = document.createElement("button");
  rm.type = "button";
  rm.className = "rm";
  rm.textContent = "✕";
  rm.title = "remover";
  rm.addEventListener("click", () => { row.remove(); refreshRefSelect(); });
  row.append(label, num, rm);
  $("comp-list").append(row);
  num.focus();
  sel.value = "";
  refreshRefSelect();
}

function refreshRefSelect() {
  const sel = $("ref");
  const keys = Object.keys(totals()).filter((k) => k !== "0002");
  const current = sel.value;
  sel.textContent = "";
  keys.forEach((k) => {
    const o = document.createElement("option");
    o.value = k;
    o.textContent = names[k] || k;
    sel.append(o);
  });
  if (keys.includes(current)) sel.value = current;
}

// ---------- cálculos ----------

let curveData = null; // último resultado de curva

function plotCurve() {
  if (!curveData) return;
  const d = curveData;
  const log = $("curve-log").checked;
  const series = [
    { name: `${names[d.ref] || d.ref} livre`, y: d.free },
    ...Object.entries(d.fractions).map(([code, y]) => ({ name: names[code] || code, y })),
  ];
  const traces = series
    .filter((s) => s.y.some((v) => v > 1e-8))
    .map((s) => ({ x: d.ph, y: s.y, mode: "lines", name: s.name }));
  const yaxis = { title: "fração do total", type: log ? "log" : "linear" };
  if (!log) yaxis.range = [0, 1.02];
  const xaxis = { title: "pH" };
  const x0 = parseFloat($("curve-xmin").value), x1 = parseFloat($("curve-xmax").value);
  if (Number.isFinite(x0) && Number.isFinite(x1)) {
    if (x1 > x0) xaxis.range = [x0, x1];
    else return status("Intervalo de pH inválido: pH máx deve ser maior que pH mín.", true);
  }
  const y0 = parseFloat($("curve-ymin").value), y1 = parseFloat($("curve-ymax").value);
  if (Number.isFinite(y0) && Number.isFinite(y1)) {
    if (y1 > y0 && (!log || y0 > 0)) yaxis.range = log ? [Math.log10(y0), Math.log10(y1)] : [y0, y1];
    else return status("Intervalo de y inválido: y máx deve ser maior que y mín (e y mín > 0 no eixo log).", true);
  }
  Plotly.newPlot("plot-dist", traces, {
    title: `Fração de ${names[d.ref] || d.ref} por espécie`,
    xaxis,
    yaxis,
    legend: { orientation: "h", y: -0.2 },   // legenda embaixo: a área do gráfico não é espremida
    margin: { t: 40, b: 90 },
  }, { responsive: true });
}

async function runCurve() {
  const t = totals();
  const ligands = Object.keys(t).filter((k) => k !== "0002");
  if (!ligands.length) return status("Adicione ao menos um componente (exceto H⁺) com total.", true);
  status("calculando curva…");
  const data = await api("/api/curve", {
    totals: t,
    kinds: kinds(),
    ph_min: parseFloat($("ph-min").value),
    ph_max: parseFloat($("ph-max").value),
    ph_points: parseInt($("ph-points").value, 10),
    ref: $("ref").value || ligands[0],
    ...activityOpts(),
  });
  Object.assign(names, data.names);
  curveData = data;
  plotCurve();
  Plotly.newPlot("plot-ht", [{ x: data.ph, y: data.H_total, mode: "lines", name: "T_H⁺ implícito" }], {
    title: "Prótons adicionados para manter o pH (negativo = base)",
    xaxis: { title: "pH" },
    yaxis: { title: "T_H⁺ (mol L⁻¹)" },
    margin: { t: 40 },
  }, { responsive: true });
  let msg = `curva pronta: ${data.ph.length} pontos.`;
  if (data.activity_model === "dh_ext") {
    const Is = data.ionic_strength;
    msg += ` Debye–Hückel: I de ${fmt(Math.min(...Is))} a ${fmt(Math.max(...Is))} mol L⁻¹.`;
  }
  status(msg);
}

function table(headers, rows) {
  const t = document.createElement("table");
  const thead = t.createTHead().insertRow();
  headers.forEach((h) => { const th = document.createElement("th"); th.textContent = h; thead.append(th); });
  const tb = t.createTBody();
  rows.forEach((cells) => {
    const tr = tb.insertRow();
    cells.forEach((c) => { const td = tr.insertCell(); td.textContent = c; });
  });
  return t;
}

// Agrupamento na tabela e no gráfico: família = primeiro componente (exceto H⁺ e água);
// dentro da família, componente livre primeiro e depois espécies por número crescente de ligantes.
function familyOf(code) {
  const nu = nuOf[code];
  if (!nu) return code;
  const comps = nu.map(([j]) => j).filter((j) => j !== WATER && j !== PROTON).sort();
  return comps.length ? comps[0] : PROTON;
}

function orderKey(e) {
  if (e.free) return -1;
  const fam = familyOf(e.code);
  const key = (nuOf[e.code] || [])
    .filter(([j]) => j !== fam)
    .reduce((s, [, n]) => s + Math.abs(n), 0);
  return key + (e.solid ? 0.5 : 0);
}

let barData = null; // último resultado de ponto: {labels, values}

function plotBars() {
  if (!barData) return;
  const log = $("bar-log").checked;
  const lo = parseFloat($("bar-min").value);
  const hi = parseFloat($("bar-max").value);
  const yaxis = { title: "concentração (mol L⁻¹)", type: log ? "log" : "linear", exponentformat: "e" };
  if (Number.isFinite(lo) && Number.isFinite(hi)) {
    if (hi <= lo || (log && lo <= 0)) {
      status("Intervalo inválido: y máx deve ser maior que y mín (e y mín > 0 no eixo log).", true);
    } else {
      // no eixo log, o Plotly usa log10 dos limites
      yaxis.range = log ? [Math.log10(lo), Math.log10(hi)] : [lo, hi];
    }
  }
  Plotly.newPlot("plot-bar", [{ x: barData.labels, y: barData.values, type: "bar", name: "concentração" }], {
    title: "Concentrações no equilíbrio",
    xaxis: { tickangle: -45, automargin: true },
    yaxis,
    margin: { t: 40, b: 120 },
  }, { responsive: true });
}

async function runPoint() {
  const body = { totals: totals(), kinds: kinds(), ...activityOpts() };
  status("calculando…");
  const r = await api("/api/speciate", body);
  Object.assign(names, r.names);
  const out = $("point-result");
  out.textContent = "";
  const h = document.createElement("p");
  const dh = r.activity_model === "dh_ext";
  h.textContent = `pH = ${fmt(r.ph)}   ·   I = ${fmt(r.ionic_strength)} mol L⁻¹` +
    (dh ? ` (Debye–Hückel, ${r.I_iterations} it. em I)` : " (γ = 1)") +
    `   ·   resíduo máx. = ${fmt(r.residual)}   ·   iterações = ${r.iterations}`;
  out.append(h);
  (r.warnings || []).forEach((w) => {
    const p = document.createElement("p");
    p.className = "warn-msg";
    p.textContent = w;
    out.append(p);
  });
  const entries = [];
  Object.entries(r.free).forEach(([code, c]) => { if (c > 1e-15) entries.push({ code, c, free: true }); });
  Object.entries(r.species).forEach(([code, c]) => { if (c > 1e-15) entries.push({ code, c }); });
  Object.entries(r.solids).forEach(([code, m]) => { if (m > 1e-15) entries.push({ code, c: m, solid: true }); });
  entries.sort((a, b) => {
    const fa = familyOf(a.code), fb = familyOf(b.code);
    if (fa !== fb) return fa < fb ? -1 : 1;
    return orderKey(a) - orderKey(b);
  });
  const labelOf = (e) => (e.solid ? `${names[e.code] || e.code} (sólido)` : names[e.code] || e.code);
  // atividade a = γ·c (γ = 1 no modelo ideal); sólidos e totais analíticos não têm atividade
  const rows = entries.map((e) => {
    if (e.solid) return [labelOf(e), fmt(e.c), "—", "—"];
    const lg = (r.log_gamma || {})[e.code] ?? 0;
    return [labelOf(e), fmt(e.c), lg.toFixed(4), fmt(e.c * 10 ** lg)];
  });
  // concentrações analíticas (totais informados), logo acima da água
  Object.entries(body.totals).forEach(([code, v]) => rows.push([`${names[code] || code} (analítica)`, fmt(v), "—", "—"]));
  rows.push(["Água", fmt(WATER_CONC), "—", "1 (convenção)"]);
  const tbl = table(["espécie", "concentração (mol L⁻¹)", "log γ", "atividade / (mol L⁻¹)"], rows);
  if (!r.converged) {
    // aviso no início da tabela; os valores do último iterado continuam exibidos
    const tr = tbl.tBodies[0].insertRow(0);
    tr.className = "warn";
    const td = tr.insertCell();
    td.colSpan = 4;
    td.textContent = r.message;
  }
  out.append(tbl);
  barData = { labels: entries.map(labelOf), values: entries.map((e) => e.c) };
  plotBars();
  status("ponto calculado.");
}

async function runMC() {
  status("Monte Carlo…");
  const r = await api("/api/mc", {
    totals: totals(),
    kinds: kinds(),
    sd: parseFloat($("mc-sd").value),
    n: 300,
    seed: 0,
    ...activityOpts(),
  });
  Object.assign(names, r.names);
  const out = $("point-result");
  out.textContent = "";
  if (r.pH) {
    const p = document.createElement("p");
    p.textContent =
      `pH = ${fmt(r.pH.mean)} ± ${fmt(r.pH.std)} (1σ, σ log β = ${r.sd_logbeta}, ` +
      `${r.n_ok} de ${r.n_ok + r.n_fail} simulações convergiram)`;
    out.append(p);
  }
  const rows = Object.entries(r.species)
    .sort((a, b) => b[1].mean - a[1].mean)
    .map(([code, s]) => [names[code] || code, fmt(s.mean), fmt(s.std)]);
  out.append(table(["espécie", "média (mol L⁻¹)", "σ (mol L⁻¹)"], rows));
  status("Monte Carlo pronto.");
}

async function showVersion() {
  let back = "?";
  try { back = (await api("/api/version")).version; } catch (e) { back = "indisponível"; }
  const tag = $("version");
  tag.textContent = back === FRONT_VERSION ? `v${FRONT_VERSION}` : `frontend v${FRONT_VERSION} ≠ backend v${back}`;
  if (back !== FRONT_VERSION) tag.style.color = "#fca5a5";
}

async function init() {
  showVersion();
  if (OPENED_AS_FILE) {
    status("Esta pagina precisa ser aberta pelo Docker: http://localhost:8090 (nao pelo arquivo index.html).", true);
    return;
  }
  try {
    // nomes de todos os componentes (para tabelas de resultado, mesmo fora do filtro atual)
    (await api("/api/components")).items.forEach((c) => { names[c.code] = c.name; });
    (await api("/api/species")).items.forEach((s) => { nuOf[s.code] = s.nu; });
    buildCategoryBox((await api("/api/categories")).items);
    updateCategorySummary();
    await loadComponents();
    $("btn-add").addEventListener("click", addSelected);
    document.querySelectorAll("input[name=act-model]").forEach((el) =>
      el.addEventListener("change", () => {
        const on = $("dh-on").checked;
        $("dh-a").disabled = !on;
        $("dh-I").disabled = !on;
      }));
    ["bar-log", "bar-min", "bar-max"].forEach((id) => $(id).addEventListener("change", plotBars));
    ["curve-log", "curve-xmin", "curve-xmax", "curve-ymin", "curve-ymax"].forEach((id) =>
      $(id).addEventListener("change", plotCurve));
    // trocar o componente recalcula a curva (as frações dependem do total de referência)
    $("ref").addEventListener("change", () => {
      if (curveData) runCurve().catch((e) => status(e.message, true));
    });
    $("btn-curve").addEventListener("click", () => runCurve().catch((e) => status(e.message, true)));
    $("btn-point").addEventListener("click", () => runPoint().catch((e) => status(e.message, true)));
    $("btn-mc").addEventListener("click", () => runMC().catch((e) => status(e.message, true)));
    status("pronto: escolha categorias, adicione componentes e informe os totais.");
  } catch (e) {
    status(`não foi possível carregar a API: ${e.message}`, true);
  }
}

init();
