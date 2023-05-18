import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

#1 Monte um arquivo do tipo fasta com a sequência de nucleotídeos de 10 genes, à sua escolha. Com Biopython, gere um segundo arquivo do tipo fasta com a transcrição e um terceiro com a tradução de todos esses genes. E se fossem 100 genes? 

out_trans = open('genes_cruzi_trans.fna', 'w')
out_trad = open('genes_cruzi_trad.faa', 'w')

for seq_record in SeqIO.parse('proteina.fasta', 'fasta'):
	rna_seq = SeqRecord(Seq(seq_record.seq.transcribe()), id=seq_record.id, description=seq_record.description)
	ptna_seq = SeqRecord(Seq(seq_record.seq.translate(table=4)), id=seq_record.id, description=seq_record.description)

	SeqIO.write(rna_seq, out_trans, 'fasta')
	SeqIO.write(ptna_seq, out_trad, 'fasta')

out_trans.close()
out_trad.close()

#parte 2 - gráfico

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#lendo o arquivo do tipo csv
cog_file = pd.read_csv('cog_freq.csv')

#transformando em um dataframe
df = pd.DataFrame(cog_file)

#criando o gráfico
sns.barplot(data=df, x='COG', y='freq', palette='dark')

#título
plt.title('Frequência de COGs')
plt.xlabel('COG')
plt.ylabel('Frequência')

#salvar figura em PNG
plt.savefig('cog.png', dpi=100)