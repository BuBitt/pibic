import os
import subprocess
from Bio import SeqIO
from tqdm import tqdm

# Diretório onde os arquivos FASTA estão armazenados (do código anterior)
fasta_dir = 'HPV_FASTA'

# Caminho para o diretório de saída dos alinhamentos
aligned_dir = 'HPV_ALIGNED'

# Criar o diretório para os alinhamentos, se não existir
os.makedirs(aligned_dir, exist_ok=True)

# Função para realizar o alinhamento usando o ClustalW2
def align_sequences(fasta_file, output_file):
    try:
        # Comando para chamar o clustalw2 via subprocess, redirecionando stdout e stderr para DEVNULL para ocultar o output
        subprocess.run(['clustalw2', fasta_file, '-OUTFILE=' + output_file, '-OUTPUT=CLUSTAL'],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Erro ao alinhar {fasta_file}: {e}")

# Contar o número total de arquivos FASTA a serem processados
total_files = 0
for type_dir in os.listdir(fasta_dir):
    type_path = os.path.join(fasta_dir, type_dir)
    if os.path.isdir(type_path):
        for lineage_dir in os.listdir(type_path):
            lineage_path = os.path.join(type_path, lineage_dir)
            if os.path.isdir(lineage_path):
                # Contar os arquivos FASTA para os genes E5, E6 e E7
                for file_name in os.listdir(lineage_path):
                    if file_name.endswith(".fasta") and file_name.split("_")[-1].replace(".fasta", "") in ["E5", "E6", "E7"]:
                        total_files += 1

# Criar a barra de progresso para todo o processo
with tqdm(total=total_files, desc="Processando alinhamentos", unit="arquivo") as pbar:
    # Percorre as pastas de linhagem dentro de HPV_FASTA
    for type_dir in os.listdir(fasta_dir):
        type_path = os.path.join(fasta_dir, type_dir)
        if os.path.isdir(type_path):
            for lineage_dir in os.listdir(type_path):
                lineage_path = os.path.join(type_path, lineage_dir)
                if os.path.isdir(lineage_path):
                    # Para cada linhagem, buscar os arquivos FASTA para os genes E5, E6 e E7
                    gene_files = {}
                    for file_name in os.listdir(lineage_path):
                        if file_name.endswith(".fasta"):
                            gene = file_name.split("_")[-1].replace(".fasta", "")
                            if gene in ["E5", "E6", "E7"]:
                                gene_files[gene] = os.path.join(lineage_path, file_name)

                    # Para cada gene, alinhar as sequências e salvar o arquivo
                    for gene, fasta_file in gene_files.items():
                        # Define o caminho do arquivo de saída de alinhamento com extensão .aln
                        aligned_lineage_path = os.path.join(aligned_dir, type_dir, lineage_dir)
                        os.makedirs(aligned_lineage_path, exist_ok=True)
                        aligned_file = os.path.join(aligned_lineage_path, f"{type_dir}_{lineage_dir}_{gene}_aligned.aln")

                        # Realiza o alinhamento se houver sequências
                        if os.path.exists(fasta_file):
                            align_sequences(fasta_file, aligned_file)
                            pbar.update(1)  # Atualiza a barra de progresso após processar cada arquivo

print("Processo de alinhamento concluído. Alinhamentos salvos em 'HPV_ALIGNED'.")
