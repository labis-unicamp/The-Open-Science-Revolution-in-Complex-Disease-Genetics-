import matplotlib.pyplot as plt
import numpy as np

# Dados convertidos para minutos (exemplo para SRR1279767)
labels = ['real', 'user', 'sys']
bwa_mem = [139.59, 129.48, 19.31]  # SRR1279767
bwa_mem2 = [87.76, 83.98, 12.15]    # SRR1279767
bowtie = [194.98, 392.02, 12.33]    # SRR1279767

angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
angles += angles[:1]  # Fechar o círculo

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'polar': True})

# Adicionar dados
def plot_radar(data, color, label):
    data += data[:1]
    ax.plot(angles, data, color=color, linewidth=2, label=label)
    ax.fill(angles, data, color=color, alpha=0.1)

plot_radar(bwa_mem, 'blue', 'bwa-mem')
plot_radar(bwa_mem2, 'red', 'bwa-mem2')
plot_radar(bowtie, 'green', 'bowtie')

# Configurações
ax.set_xticks(angles[:-1])
ax.set_xticklabels(labels)
ax.set_title('Comparação de Tempos (SRR1279767)', pad=20)
ax.legend(loc='upper right')

plt.show()
