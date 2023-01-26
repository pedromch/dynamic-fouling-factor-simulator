from calculations import calculate
from utils import get_data_to_fit
from plot_data import plot_data
from inputs import *

data = get_data_to_fit(filename = nome_arquivo_de_dados, dirpath = diretorio_arquivo_de_dados)

results = calculate(
    Ds = diametro_casco,
    Dte = diametro_externo_tubo,
    Npt = numero_passes_tubos,
    Rtp = razao_passo_tubos,
    layout = arranjo_do_feixe,
    Ltube = comprimento_tubos,
    esp = espessura_tubos,
    ktube = condutividade_termica_tubos,
    Nc = numero_chicanas,
    fluid_inside_tubes = fluido_lado_dos_tubos,
    data_to_fit = data,
    tempo_final = tempo_final
)

plot_data(*results, dirpath = diretorio_de_saida, filename = nome_imagem_de_saida)