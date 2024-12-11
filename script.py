import os
import csv
from Bio import Entrez
from Bio import SeqIO
from tqdm import tqdm

# Configure o seu e-mail para o Entrez
Entrez.email = "seu_email@exemplo.com"

# Nome do arquivo CSV
csv_file = 'hpv_lineages.csv'

# Nome do diretório raiz onde as pastas serão criadas
root_dir = 'HPV'

# Criar o diretório raiz se ele não existir
os.makedirs(root_dir, exist_ok=True)

# Lê o arquivo CSV e prepara a lista de linhas
with open(csv_file, mode='r', encoding='utf-8') as file:
    reader = list(csv.DictReader(file))

# Barra de progresso
with tqdm(total=len(reader), desc="Processando dados", unit="linha") as pbar:
    for row in reader:
        # Obtem os valores das colunas 'type', 'lineage' e 'gb_accession_number'
        type_dir = row['type']
        lineage_dir = row['lineage']
        accession_number = row['gb_accession_number']

        # Define o caminho da pasta
        folder_path = os.path.join(root_dir, type_dir, lineage_dir)

        # Cria a hierarquia de pastas
        os.makedirs(folder_path, exist_ok=True)

        # Define o caminho completo para salvar o arquivo GenBank
        genbank_file = os.path.join(folder_path, f"{accession_number}.gb")

        # Baixa o arquivo GenBank se ainda não existir
        if not os.path.exists(genbank_file):
            try:
                with Entrez.efetch(db="nuccore", id=accession_number, rettype="gb", retmode="text") as handle:
                    record = handle.read()

                # Salva o arquivo GenBank
                with open(genbank_file, "w") as output_file:
                    output_file.write(record)
            except Exception as e:
                print(f"Erro ao baixar {accession_number}: {e}")
        
        # Atualiza a barra de progresso
        pbar.update(1)

print(f"Processo concluído. Hierarquia de pastas e arquivos GenBank criados em: {root_dir}")