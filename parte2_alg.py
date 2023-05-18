import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import glob
import os

with open('lista_genes.txt') as arquivo:
    lista = list(arquivo)
    
out_genes = open('nts_lista.fasta', 'w')
out_proteinas = open('ptns_lista.fasta', 'w')
out_anot = open('anot_lista.gbk', 'w')

for item in lista:
    item = str(item.strip('\n'))
    print(lista)

    #iterando sobre o arquivo gbff
    for file in glob.glob('*.gbff'):
        for record in SeqIO.parse(file, 'gb'):
            #percorre as features do arquivo gbff
            for feature in record.features:
                #verifica se a feature é um gene é se possui a qualifier 'pseudo'
                if feature.type == 'CDS':
                    if 'product' in feature.qualifiers: #verificar se tem esse qualifier
                        if item in feature.qualifiers['product']:
                            #acessando a sequencia do gene pela coordenada
                            location = feature.location
                            start = feature.location.start
                            end = feature.location.end
                            coord = SeqFeature(FeatureLocation(start, end))

                            #arquivo de nucleotideos
                            seq = coord.extract(feature.seq)
                            idee = feature.qualifiers['product'][0]
                            desc = feature.qualifiers['product'][0] + f"[{file}]"
                            gene = SeqRecord(Seq(seq), id=idee, description=desc)
                            SeqIO.write(gene, out_genes, 'fasta')

                            #arquivo de aminoacidos
                            seq = feature.qualifiers['translation'][0]
                            idee = feature.qualifiers['protein_id'][0]
                            desc = feature.qualifiers['product'][0] + f"[{file}]"
                            proteina = SeqRecord(Seq(seq), id=idee, description=desc)
                            SeqIO.write(proteina, out_proteinas, 'fasta')

                            #arquivo de anotacao
                            seq = coord.extract(record.seq)
                            location = feature.location
                            protein_id = feature.qualifiers['protein_id'][0]
                            translation = feature.qualifiers['translation'][0]
                            feature = SeqFeature(location, type='CDS', qualifiers={'gene':item, 'protein_id':protein_id, 'translation':translation})

                            record_anot = SeqRecord(Seq(seq), id=file, description=f"Anotação{item}", annotations={'molecule_type':'DNA'})
                            record_anot.features.append(feature)

                            SeqIO.write(record_anot, out_anot, 'genbank')

out_genes.close()
out_proteinas.close()
out_anot.close()