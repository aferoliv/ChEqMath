const $ = (id) => document.getElementById(id);
const OPENED_AS_FILE = location.protocol === "file:";
let components = [];
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

function kinds() {
  return [...document.querySelectorAll("input[name=kind]:checked")].map((x) => x.value);
}

function totals() {
  const out = {};
  document.querySelectorAll("#comp-list .row").forEach((row) => {
    const cb = row.querySelector("input[type=checkbox]");
    const num = row.querySelector("input[type=number]");
    if (cb.checked && num.value !== "") out[cb.dataset.code] = parseFloat(num.value);
  });
  return out;
}

function buildList() {
  const box = $("comp-list");
  box.textContent = "";
  components.filter((c) => !c.is_solvent).forEach((c) => {
    const row = document.createElement("div");
    row.className = "row";
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.dataset.code = c.code;
    const label = document.createElement("span");
    label.textContent = `${c.name} (${c.code})`;
    const num = document.createElement("input");
    num.type = "number";
    num.step = "any";
    num.placeholder = "mol/L";
    num.disabled = true;
    cb.addEventListener("change", () => { num.disabled = !cb.checked; });
    row.append(cb, label, num);
    box.append(row);
  });
  $("filter").addEventListener("input", (e) => {
    const q = e.target.value.toLowerCase();
    document.querySelectorAll("#comp-list .row").forEach((r) => {
      r.style.display = r.textContent.toLowerCase().includes(q) ? "" : "none";
    });
  });
}

function refreshRefSelect() {
  const sel = $("ref");
  const keys = Object.keys(totals()).filter((k) => k !== "0002");
  const current = sel.value;
  sel.textContent = "";
  keys.forEach((k) => {
    const o = document.createElement("option");
    o.value = k;
    o.textContent = `${names[k] || k} (${k})`;
    sel.append(o);
  });
  if (keys.includes(current)) sel.value = current;
}

async function runCurve() {
  const t = totals();
  const ligands = Object.keys(t).filter((k) => k !== "0002");
  if (!ligands.length) return status("Marque ao menos um componente (exceto H⁺) com total.", true);
  status("calculando curva…");
  const data = await api("/api/curve", {
    totals: t,
    kinds: kinds(),
    ph_min: parseFloat($("ph-min").value),
    ph_max: parseFloat($("ph-max").value),
    ph_points: parseInt($("ph-points").value, 10),
    ref: $("ref").value || ligands[0],
  });
  Object.assign(names, data.names);
  const traces = Object.entries(data.fractions)
    .filter(([, arr]) => arr.some((v) => v > 1e-8))
    .map(([code, arr]) => ({ x: data.ph, y: arr, mode: "lines", name: names[code] || code }));
  Plotly.newPlot("plot-dist", traces, {
    title: `Fração de ${names[data.ref] || data.ref} por espécie`,
    xaxis: { title: "pH" },
    yaxis: { title: "fração do total", range: [0, 1.02] },
    margin: { t: 40 },
  }, { responsive: true });
  Plotly.newPlot("plot-ht", [{ x: data.ph, y: data.H_total, mode: "lines", name: "T_H⁺ implícito" }], {
    title: "Prótons adicionados para manter o pH (negativo = base)",
    xaxis: { title: "pH" },
    yaxis: { title: "T_H⁺ (mol L⁻¹)" },
    margin: { t: 40 },
  }, { responsive: true });
  status(`curva pronta: ${data.ph.length} pontos.`);
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

async function runPoint() {
  const body = { totals: totals(), kinds: kinds() };
  status("calculando…");
  const r = await api("/api/speciate", body);
  Object.assign(names, r.names);
  const out = $("point-result");
  out.textContent = "";
  const ph = r.ph === null ? "—" : r.ph.toFixed(3);
  const h = document.createElement("p");
  h.textContent = `pH = ${ph}   ·   resíduo máx. = ${r.residual.toExponential(1)}   ·   iterações = ${r.iterations}`;
  out.append(h);
  const rows = Object.entries(r.species)
    .filter(([, c]) => c > 1e-15)
    .sort((a, b) => b[1] - a[1])
    .map(([code, c]) => [names[code] || code, code, c.toExponential(3)]);
  rows.push(...Object.entries(r.solids).map(([code, m]) => [names[code] || code, code + " (sólido)", m.toExponential(3)]));
  out.append(table(["espécie", "código", "concentração (mol L⁻¹)"], rows));
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
  });
  Object.assign(names, r.names);
  const out = $("point-result");
  out.textContent = "";
  if (r.pH) {
    const p = document.createElement("p");
    p.textContent =
      `pH = ${r.pH.mean.toFixed(3)} ± ${r.pH.std.toFixed(3)} (1σ, σ log β = ${r.sd_logbeta}, ` +
      `${r.n_ok} de ${r.n_ok + r.n_fail} simulações convergiram)`;
    out.append(p);
  }
  const rows = Object.entries(r.species)
    .sort((a, b) => b[1].mean - a[1].mean)
    .map(([code, s]) => [names[code] || code, code, s.mean.toExponential(3), s.std.toExponential(2)]);
  out.append(table(["espécie", "código", "média (mol L⁻¹)", "σ (mol L⁻¹)"], rows));
  status("Monte Carlo pronto.");
}

async function init() {
  if (OPENED_AS_FILE) {
    status("Esta pagina precisa ser aberta pelo Docker: http://localhost:8090 (nao pelo arquivo index.html).", true);
    return;
  }
  try {
    components = (await api("/api/components")).items;
    components.forEach((c) => { names[c.code] = c.name; });
    buildList();
    $("comp-list").addEventListener("change", refreshRefSelect);
    $("btn-curve").addEventListener("click", () => runCurve().catch((e) => status(e.message, true)));
    $("btn-point").addEventListener("click", () => runPoint().catch((e) => status(e.message, true)));
    $("btn-mc").addEventListener("click", () => runMC().catch((e) => status(e.message, true)));
    status("pronto: marque componentes e informe os totais.");
  } catch (e) {
    status(`não foi possível carregar a API: ${e.message}`, true);
  }
}

init();
