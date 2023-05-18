#desafio parte 2 - gráfico

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#lendo o arquivo do tipo csv
cog_file = pd.read_csv('cog_freq2.csv')

#transformando em um dataframe
df = pd.DataFrame(cog_file)

#criando o gráfico
sns.barplot(data=df, x='COG', y='freq', palette='pastel')

#título
plt.title('Frequência de COGs')
plt.xlabel('COG')
plt.ylabel('Frequência')

#salvar figura em PNG
plt.savefig('cog.png', dpi=100)