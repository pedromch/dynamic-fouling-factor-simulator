diametro_casco = 938.8e-3 # m
diametro_externo_tubo = 25.4e-3 # m
numero_passes_tubos = 2
razao_passo_tubos = 1.25
arranjo_do_feixe = "triangular" # quadrado ou triangular
comprimento_tubos = 3.0488 # m
espessura_tubos = 2.769e-3 # m
condutividade_termica_tubos = 50 # W/mK
numero_chicanas = 16
fluido_lado_dos_tubos = "frio" # frio ou quente

tempo_final = 1000 # tempo em dias para estimativa de Tho, Tco, U e Q ao longo do tempo, se None, plota apenas as curvas de Rft

diretorio_arquivo_de_dados = None # endereço do diretório que contém arquivo com os dados de deposição, se None, checa em assets
nome_arquivo_de_dados = "HEXdataRev1.xlsx" # nome do arquivo de dados de deposição

diretorio_de_saida = None # endereço do diretório de saída, se None, salva na pasta de assets
nome_imagem_de_saida = "resultados.jpeg" # nome do arquivo de imagem de saída