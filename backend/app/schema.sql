-- Esquema do banco de especiação (gerado por app/importer.py a partir da planilha ChEq).
-- Modelo de Teoria.docx: componentes j, espécies i, estequiometria ν_ij.

CREATE TABLE components (
    code       TEXT PRIMARY KEY,           -- '0002'
    name       TEXT NOT NULL,
    smiles     TEXT,
    is_solvent INTEGER NOT NULL DEFAULT 0  -- 1 para a água (a = 1)
);

CREATE TABLE species (
    code         TEXT PRIMARY KEY,         -- 'SP0001' (ID próprio; nomes podem repetir)
    name         TEXT NOT NULL,
    kind         TEXT NOT NULL CHECK (kind IN ('ABB', 'complex', 'S0')),
    log_beta     REAL NOT NULL,            -- valor usado no cálculo (log K ou log So para S0)
    log_beta_raw REAL,                     -- valor original da planilha
    cap_of       TEXT REFERENCES species(code),  -- S0: espécie aquosa limitada
    water_ref    INTEGER NOT NULL DEFAULT 0 CHECK (water_ref IN (0, 1)),  -- 1: constante referida a água 55,6 mol/L
    source       TEXT NOT NULL,            -- 'planilha:Species' | 'gerada:pKa'
    note         TEXT                      -- correções e observações
);

CREATE TABLE stoich (
    species_code TEXT NOT NULL REFERENCES species(code),
    comp_code    TEXT NOT NULL REFERENCES components(code),
    nu           INTEGER NOT NULL CHECK (nu <> 0),
    PRIMARY KEY (species_code, comp_code)
);

CREATE TABLE pka (
    comp_code TEXT NOT NULL REFERENCES components(code),
    idx       INTEGER NOT NULL,            -- 1 = menor pKa
    value     REAL NOT NULL,
    carga_max INTEGER,
    PRIMARY KEY (comp_code, idx)
);

CREATE TABLE import_log (
    id      INTEGER PRIMARY KEY AUTOINCREMENT,
    level   TEXT NOT NULL,                 -- info | correção | alerta | erro
    message TEXT NOT NULL
);
