# ChEqMath - Chemical Equilibria Mathematics

Especiação em solução aquosa.

Motor de equilíbrio químico baseado no modelo **componentes–espécies** do `Teoria.docx`,
com banco gerado a partir da planilha `componentes_e_especies_ChEq_Math.xlsx`.

## Estrutura

```
compose.yaml              orquestra api + web
data/                     planilha de entrada e speciation.sqlite (gerado)
backend/
  app/
    engine.py             motor: solver de Newton, sólidos S0, curva de pH, Monte Carlo
    importer.py           planilha → SQLite (regras e correções registradas em import_log)
    schema.sql            esquema do banco
    db.py                 leitura do banco → objetos Species
    main.py               API FastAPI (/api/*)
  tests/                  motor, importador e API
frontend/                 HTML/JS + Plotly; nginx faz proxy de /api/
```

## Rodar com Docker

```bash
docker compose up --build
# abra http://localhost:8090
```

## Desenvolvimento local

```bash
cd backend
pip install -r requirements-dev.txt
python -m app.importer --xlsx ../data/componentes_e_especies_ChEq_Math.xlsx --db ../data/speciation.sqlite
pytest -q
uvicorn app.main:app --reload --port 8000
```

## Convenções

- `a_i = 10^log_β · Π a_j^ν_ij`, com constantes termodinâmicas (I = 0) a 25 °C.
- Atividade: `a_i = γ_i c_i`. Padrão `ideal` (γ = 1). Opção `dh_ext` (Debye–Hückel estendida):
  `log γ_i = −A z_i² √I / (1 + B å_i √I)`, A = 0,5091 e B = 0,3283 Å⁻¹ (25 °C), å padrão 4,5 Å
  (ajustável; `ion_sizes` por código). `I = ½ Σ c_k z_k²` é resolvida por ponto fixo externo ao
  Newton, ou fixada pelo usuário (`ionic_strength`, meio iônico). Aviso quando I > 0,1 mol L⁻¹.
- `z_i = Σ ν_ij z_j`. **Coluna 'Maximum Electrical Charge' da aba components** (achada pelo
  cabeçalho; as colunas 'code' e 'category …' também) = carga da espécie totalmente protonada:
  Cd(II) `2+`, glicina `+`, acetato `0`. Como o componente é a forma totalmente desprotonada,
  `z_j = q_max − n_H`, com `n_H` = maior entre o nº de pKa da aba pKa e o maior k das espécies
  HₖX da aba Species (metais: `n_H = 0`). Célula vazia: metal pelo nome, com alerta; demais
  ficam sem carga e o cálculo com Debye–Hückel os recusa.
- Íons espectadores (Na⁺, Cl⁻ etc.) só entram em I se forem componentes do sistema; caso
  contrário, use I fixa.
- pH = −log a_H+; pH fixo (curva) fixa a **atividade** de H⁺.
- Balanço de massa: `T_j = [C_j] + Σ ν_ij c_i + Σ ν_sj m_s`.
- O componente de um ligante é a **forma totalmente desprotonada** (ex.: `0003` = acetato).
- `0001` é a água (a = 1); `0002` é H⁺, cujo total pode ser negativo (adição de base).
- Sólidos S0 limitam a espécie aquosa a `So = 10^log_β(S0)`; o excesso precipita.
- Incerteza: Monte Carlo sobre log β com σ escolhido pelo usuário (padrão 0,05). Esse σ não vem da planilha.

## Correções aplicadas à planilha (ver `import_log`)

| Espécie | Planilha | Usado | Motivo |
|---|---|---|---|
| hidróxido (OH⁻) | +15,75 | −14,00 | Sinal inconsistente com a estequiometria; −pKw a 25 °C |
| Cd(OH)ₙ, n = 1–4 | −10,08; −10,27; −12,95; −14,05 (stepwise) | −10,08; −20,35; −33,30; −47,35 (cumulativo) | O Teoria define K como constante de formação a partir dos componentes |

Confirme as correções de cádmio antes de usar os resultados.
