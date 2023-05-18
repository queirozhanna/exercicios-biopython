import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import glob
import os 

#Baixe o arquivo de um genoma do NCBI a sua escolha. Selecione um gene e, com Biopython, recupere a sequencia de nucleotideos, de aminoacidos e a anotacao desse gene. Por fim, crie um arquivo para armazenar a sequencia de 10 nucleotideos, de aminoacidos e para a anotacao desse gene. E se fossem 10 genes? Em 10 genomas diferentes?

out1 = open('nt_gene.fasta', 'w')
out2 = open('aas_gene.fasta', 'w')
out3 = open('anot_gene.gbff', 'w')

for file in glob.glob('*.gbff'):
	for record in SeqIO.parse(file, 'gb'): #iterando sobre o arquivo gbff
		for feature in record.features: #percorre as features do arquivo gbff
			if feature.type == 'CDS': #verifica se a feature eh um gene e se tem o qualifier 'pseudo'
				if 'product' in feature.qualifiers:
					if 'hexokinase' in feature.qualifiers['product']:
						#acessando a sequencia do gene pela coordenada
						location = feature.location
						start = feature.location.start
						end = feature.location.end
						coord = SeqFeature(FeatureLocation(start, end))
						
						#arquivo de nucleotideos
						
						seq = coord.extract(record.seq)
						idee = feature.qualifiers['protein_id'][0]
						desc = feature.qualifiers['product'][0] + f' [{file}]'
						gene = SeqRecord(seq, id=idee, description=desc)
						
						SeqIO.write(gene, out1, 'fasta')
						
						#arquivo de aminoacidos
						
						seq = feature.qualifiers['translation'][0]
						idee = feature.qualifiers['protein_id'][0]
						desc = feature.qualifiers['product'][0] + f' [{file}]'
						proteina = SeqRecord(Seq(seq), id=idee, description=desc)
						
						SeqIO.write(proteina, out2, 'fasta')
						
						#arquivo de anotacao
						
						seq = coord.extract(record.seq)
						location = location = feature.location
						protein_id = feature.qualifiers['protein_id'][0]
						translation = feature.qualifiers['translation'][0]
						feature = SeqFeature(location, type='CDS', qualifiers={'gene':'HK', 'protein_id':protein_id, 'translation':translation})
						record_anot = SeqRecord(Seq(seq), id=file, description='Anotacao HK', annotations={'molecule_type':'DNA'})
						record_anot.features.append(feature)
						
						SeqIO.write(record_anot, out3, 'genbank')
out1.close()
out2.close()
out3.close()
