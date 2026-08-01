# Protein Analysis Tools

Este repositório contém scripts em Python para análise físico-química de sequências de proteínas, com cálculo de propriedades básicas (massa molecular, ponto isoelétrico, índice de instabilidade, coeficiente de extinção) e predição de antigenicidade (VaxiJen).  Os resultados incluem tabelas e gráficos gerados automaticamente, baseados em construções proteicas utilizadas no meu mestrado.

---

##  Objetivos
- Automatizar a análise físico-química de proteínas a partir de arquivos FASTA.
- Calcular propriedades relevantes para caracterização de proteínas recombinantes.
- Integrar dados experimentais (ex.: scores de antigenicidade) com análises computacionais.
- Gerar relatórios e gráficos de forma reprodutível.

---

---

## Tecnologias utilizadas
- **Linguagem:** Python 3.9+
- **Bibliotecas:** Biopython, pandas, matplotlib, seaborn
- **Ambiente:** Linux/Unix, Google Colab ou qualquer IDE Python

---

## Como usar
1. Clone o repositório:
   ```bash
   git clone https://github.com/bsastolfi-maker/protein-analysis-tools.git
   cd protein-analysis-tools

  pip install -r requirements.txt
  python analyzer.py

