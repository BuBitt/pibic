import os
from Bio import AlignIO
from tqdm import tqdm

# Caminho para o diretório onde os arquivos alinhados estão armazenados
aligned_dir = 'HPV_ALIGNED'
output_dir = 'HPV_DIVERGENCES'

# Função para comparar as sequências e identificar as posições divergentes
def find_divergent_positions(alignment):
    divergent_positions = []
    num_sequences = len(alignment)
    length = alignment.get_alignment_length()

    # Percorrer cada posição do alinhamento
    for i in range(length):
        # Extrair os nucleotídeos em cada sequência naquela posição
        column = [alignment[j, i] for j in range(num_sequences)]
        # Verificar se há divergência (se não são todos iguais)
        if len(set(column)) > 1:  # Se houver pelo menos 2 nucleotídeos diferentes
            divergent_positions.append(i + 1)  # Posições começam de 1
    return divergent_positions, column

# Função para gerar um arquivo com as posições divergentes e suas sequências
def generate_divergence_report(alignment_file, output_file):
    try:
        alignment = AlignIO.read(alignment_file, "clustal")

        # Encontre as posições divergentes
        divergent_positions, _ = find_divergent_positions(alignment)

        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Abrir o arquivo para escrita
        with open(output_file, 'w') as f:
            f.write("Lista de posições divergentes: " + ", ".join(map(str, divergent_positions)) + "\n")
            for pos in divergent_positions:
                f.write(f"Posição {pos}:\n")
                for seq_record in alignment:
                    f.write(f"  {seq_record.id}: {seq_record.seq[pos-1]}\n")

    except ValueError:
        pass  # Não imprime nada para arquivos vazios ou malformados

# Função principal para processar os arquivos de alinhamento
def process_alignments():
    total_files = 0
    files_to_process = []

    # Contar o número total de arquivos de alinhamento
    for type_dir in os.listdir(aligned_dir):
        type_path = os.path.join(aligned_dir, type_dir)
        if os.path.isdir(type_path):
            for lineage_dir in os.listdir(type_path):
                lineage_path = os.path.join(type_path, lineage_dir)
                if os.path.isdir(lineage_path):
                    # Coletar todos os arquivos .aln
                    for file_name in os.listdir(lineage_path):
                        if file_name.endswith(".aln"):
                            alignment_file = os.path.join(lineage_path, file_name)
                            files_to_process.append(alignment_file)
                            total_files += 1

    # Processar todos os arquivos de alinhamento com barra de progresso
    with tqdm(total=total_files, desc="Processando alinhamentos") as pbar:
        for alignment_file in files_to_process:
            # Determinar o caminho de saída correspondente, mantendo a mesma hierarquia
            relative_path = os.path.relpath(alignment_file, aligned_dir)  # Caminho relativo para manter a estrutura
            output_file = os.path.join(output_dir, relative_path.replace(".aln", "_divergences.txt"))

            generate_divergence_report(alignment_file, output_file)
            pbar.update(1)  # Atualiza a barra de progresso após cada arquivo

# Executar o processo
process_alignments()
