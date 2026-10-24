FROM python:3.11-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

ENV PYTHONUNBUFFERED=1 \
    DB_PATH=/app/data/speciation.sqlite \
    XLSX_PATH=/app/data/componentes_e_especies_ChEq_Math.xlsx

EXPOSE 8000
# Reconstrói o banco a partir da planilha a cada subida e inicia a API
CMD ["sh", "-c", "python -m app.importer && uvicorn app.main:app --host 0.0.0.0 --port 8000"]
