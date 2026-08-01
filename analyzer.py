# -*- coding: utf-8 -*-
"""
Protein Physicochemical Analyzer
--------------------------------
Script em Python para análise físico-química de construções proteicas
utilizando Biopython (ProtParam) e visualização com matplotlib/seaborn.

Autor: Bianca S Astolfi
"""

# ==============================
# Dependências
# ==============================
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ==============================
# Funções 
# ==============================
def analyze_protein(sequence: str, protein_id: str) -> dict:
    """Analisa propriedades físico-químicas de uma proteína."""
    analysed = ProteinAnalysis(sequence)

    num_aa = len(sequence)
    mw_kda = round(analysed.molecular_weight() / 1000, 2)
    pi = round(analysed.isoelectric_point(), 2)
    ext_coeff = analysed.molar_extinction_coefficient()[0]
    instability_index = round(analysed.instability_index(), 2)
    stability = "Estável" if instability_index < 40 else "Instável"

    return {
        "Proteína": protein_id,
        "Nº Aminoácidos": num_aa,
        "Massa Molecular (kDa)": mw_kda,
        "pI Teórico": pi,
        "Coef. Extinção (M^-1 cm^-1)": ext_coeff,
        "Índice de Instabilidade": instability_index,
        "Classificação": stability
    }

def run_analysis(fasta_path: str, output_csv: str = "physicochemical_properties.csv") -> pd.DataFrame:
    """Executa análise físico-química a partir de um arquivo FASTA."""
    results = []
    fasta_file = Path(fasta_path)

    if not fasta_file.exists():
        raise FileNotFoundError(f"Arquivo {fasta_path} não encontrado!")

    for record in SeqIO.parse(fasta_file, "fasta"):
        results.append(analyze_protein(str(record.seq), record.id))

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"✅ Análise concluída! Resultados salvos em {output_csv}")
    return df

def plot_vaxijen(df: pd.DataFrame, scores: dict, output_img: str = "vaxijen_antigenicity.png"):
    """Gera gráfico comparativo de antigenicidade (VaxiJen)."""
    df["Score VaxiJen"] = df["Proteína"].map(scores)

    sns.set_theme(style="whitegrid", palette="muted")
    plt.figure(figsize=(10, 5))
    colors = ['#2b5c8f' if x < max(scores.values()) else '#d9534f' for x in df["Score VaxiJen"]]

    ax = sns.barplot(data=df, x="Proteína", y="Score VaxiJen", palette=colors)
    plt.axhline(0.40, color='red', linestyle='--', linewidth=1.5, label='Limiar de Antigenicidade (0.40)')

    for p in ax.patches:
        ax.annotate(f'{p.get_height():.2f}',
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center',
                    xytext=(0, 7), textcoords='offset points',
                    fontsize=10, weight='bold')

    plt.title("Predição de Antigenicidade (VaxiJen v2.0)", fontsize=13, weight='bold')
    plt.xlabel("Proteína", weight='bold')
    plt.ylabel("Score VaxiJen", weight='bold')
    plt.legend(loc='upper left')
    plt.savefig(output_img, dpi=300)
    plt.show()
    print(f"📊 Gráfico salvo em {output_img}")

# ==============================
# Execução principal
# ==============================
if __name__ == "__main__":
    # Caminho do FASTA
    fasta_path = "data/processed/constructs.fasta"

    # Rodar análise
    df = run_analysis(fasta_path)

    # Scores de antigenicidade (exemplo da dissertação)
    vaxijen_scores = {
        "RBD": 0.46,
        "RBD-E": 0.54,
        "RBD-M": 0.58,
        "E-RBD-E": 0.55,
        "E-RBD": 0.51,
        "E-RBD-M": 0.61,
        "M-RBD-M": 0.60,
        "M-RBD": 0.49,
        "M-RBD-E": 0.54
    }

    # Gerar gráfico
    plot_vaxijen(df, vaxijen_scores)
