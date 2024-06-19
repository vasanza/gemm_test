import csv
import matplotlib.pyplot as plt

# Leer los datos del archivo CSV
m_values = []
my_gflops = []
ref_gflops = []

with open('results.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    next(csvreader)  # Saltar la cabecera
    for row in csvreader:
        m_values.append(int(row[0]))
        my_gflops.append(float(row[3]))
        ref_gflops.append(float(row[4]))

# Graficar los resultados
plt.plot(m_values, my_gflops, label='Version5')
plt.plot(m_values, ref_gflops, label='Version1')
plt.xlabel('Matrix Dimension (m)')
plt.ylabel('GFLOPS')
plt.title('Performance Comparison')
plt.legend()
plt.grid(True)
plt.savefig('performance_comparison.png')
plt.show()
