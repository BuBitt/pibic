import subprocess
from tqdm import tqdm
import sys

def run_script(script_path):
    try:
        # Executa o script e permite que a saída seja mostrada em tempo real
        with subprocess.Popen(['python', script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
            # Cria uma barra de progresso específica para a execução do script
            with tqdm(total=100, desc=f"Executando {script_path}", ncols=100, file=sys.stdout) as pbar:
                for line in proc.stdout:
                    # Exibe a saída do script em tempo real
                    print(line, end='', flush=True)
                    # Atualiza a barra de progresso, se for aplicável
                    pbar.update(1)
                # Captura qualquer erro
                stderr = proc.stderr.read()
                if stderr:
                    print(f"Erros do script {script_path}:\n{stderr}", file=sys.stderr)
            proc.wait()  # Espera o processo terminar
    except Exception as e:
        print(f"Erro ao executar o script {script_path}: {str(e)}", file=sys.stderr)

def execute_scripts_sequentially(scripts):
    for script in scripts:
        run_script(script)

# Lista de scripts a serem executados
scripts = ['nuccore_getter.py', 'alignment.py', 'divergences.py']

# Executando os scripts um após o outro, exibindo as barras de progresso de cada um
execute_scripts_sequentially(scripts)
