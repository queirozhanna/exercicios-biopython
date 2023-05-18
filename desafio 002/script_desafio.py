#Baixe a sequência de DNA de 10 genomas completos do NCBI, à sua escolha. Com Biopython, calcule o conteúdo GC de cada um dos genomas e crie um gráfico representando esses valores. E se fossem 500 genomas?

import Bio
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

out = open('gc_content.txt', 'w')

for file in glob.glob('*.fna'):
	for seq_record in SeqIO.parse(file, 'fasta'):
		gc = gc_fraction(seq_record.seq)
		escreva = file + '\t' + str(gc) + '\n'
		out.writelines(escreva)
		
out.close()

#eixos
genomas = ['ecoli1', 'ecoli2', 'ecoli3', 'ecoli4', 'ecoli5', 'mycotube1', 'mycotube2', 'mycotube3', 'mycotube4', 'mycotube5']
gc = [0.5074638178374098, 0.506221885595354, 0.5110640060403953, 0.510370993871227, 0.5102192882964885, 0.6531425521863135, 0.6043227049004922, 0.6557536661961154, 0.639344262295082, 0.6556060619918235]

#estilo do tema

plt.style.use('seaborn-pastel')

#barras na vertical
plt.bar(genomas, gc, color = 'hotpink')

plt.ylabel("% GC")
plt.xticks(rotation=45, fontsize=8)
plt.xlabel("Genomas")
plt.title("Conteúdo GC de genomas de Escherichia coli e Mycobacterium tuberculosis")
plt.show()


