import os
import csv
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import time

# Configure o seu e-mail para o Entrez
Entrez.email = "seu_email@exemplo.com"

# Nome do arquivo CSV
csv_file = 'hpv_lineages.csv'

# Nome do diretório raiz onde as pastas serão criadas
root_dir = 'HPV'
fasta_dir = 'HPV_FASTA'

# Criar os diretórios raiz se eles não existirem
os.makedirs(root_dir, exist_ok=True)
os.makedirs(fasta_dir, exist_ok=True)

# Lê o arquivo CSV e prepara a lista de linhas
with open(csv_file, mode='r', encoding='utf-8') as file:
    reader = list(csv.DictReader(file))

# Barra de progresso para download dos arquivos GenBank
with tqdm(total=len(reader), desc="Baixando arquivos GenBank", unit="linha", mininterval=1) as pbar:
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
                    with open(genbank_file, "w") as output_file:
                        output_file.write(handle.read())
            except Exception as e:
                print(f"Erro ao baixar {accession_number}: {e}")
        pbar.update(1)
        time.sleep(0.1)  # Pausa breve para dar tempo à atualização da barra

# Barra de progresso para o processamento dos arquivos GenBank
with tqdm(total=len(reader), desc="Processando arquivos GenBank", unit="linha", mininterval=1) as pbar:
    for row in reader:
        type_dir = row['type']
        lineage_dir = row['lineage']
        accession_number = row['gb_accession_number']

        # Define o caminho da pasta
        folder_path = os.path.join(root_dir, type_dir, lineage_dir)

        # Dicionário para armazenar sequências por gene
        gene_sequences = {"E5": [], "E6": [], "E7": []}

        # Processar todos os arquivos GenBank na pasta
        for file_name in os.listdir(folder_path):
            if file_name.endswith(".gb"):
                genbank_file = os.path.join(folder_path, file_name)
                try:
                    with open(genbank_file, "r") as input_handle:
                        gb_record = SeqIO.read(input_handle, "genbank")

                    # Extrair sequências CDS específicas
                    for feature in gb_record.features:
                        if feature.type == "CDS":
                            gene = feature.qualifiers.get("gene", [None])[0]
                            product = feature.qualifiers.get("product", [None])[0]

                            # Verificar se o gene é E5, E6 ou E7
                            if gene and gene in gene_sequences:
                                sequence = feature.location.extract(gb_record).seq
                                record_id = f"{gb_record.id}_{gene}"
                                description = f"{product}" if product else "Sem descrição"
                                gene_sequences[gene].append(SeqRecord(sequence, id=record_id, description=description))
                            elif product and any(x in product.lower() for x in ["e5", "e6", "e7"]):
                                # Se o produto do gene for relacionado a E5, E6 ou E7, atribuí-lo corretamente
                                sequence = feature.location.extract(gb_record).seq
                                record_id = f"{gb_record.id}_{product}"
                                description = f"{product}" if product else "Sem descrição"
                                if "e5" in product.lower():
                                    gene_sequences["E5"].append(SeqRecord(sequence, id=record_id, description=description))
                                elif "e6" in product.lower():
                                    gene_sequences["E6"].append(SeqRecord(sequence, id=record_id, description=description))
                                elif "e7" in product.lower():
                                    gene_sequences["E7"].append(SeqRecord(sequence, id=record_id, description=description))

                except Exception as e:
                    print(f"Erro ao processar {genbank_file}: {e}")

        # Criar hierarquia de pastas para FASTA
        fasta_lineage_path = os.path.join(fasta_dir, type_dir, lineage_dir)
        os.makedirs(fasta_lineage_path, exist_ok=True)

        # Criar arquivos FASTA combinados para cada gene, garantindo que apenas E5, E6 ou E7 sejam incluídos
        for gene, sequences in gene_sequences.items():
            if sequences:
                fasta_file = os.path.join(fasta_lineage_path, f"HPV_{type_dir}_{lineage_dir}_{gene}.fasta")
                with open(fasta_file, "w") as output_handle:
                    SeqIO.write(sequences, output_handle, "fasta")
        pbar.update(1)
        time.sleep(0.1)  # Pausa breve para dar tempo à atualização da barra

print(f"Processo concluído. Arquivos FASTA criados com sequências CDS específicas para cada gene por linhagem.")
