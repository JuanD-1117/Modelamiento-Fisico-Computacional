import pandas as pd
import matplotlib.pyplot as plt

# Cargar archivos
try:
    v_data = pd.read_csv('potencial.csv')
    e_levels = pd.read_csv('niveles.csv')
except:
    print("Faltan archivos CSV. Ejecuta primero el código de Fortran.")
    exit()

plt.figure(figsize=(10, 7))

# Graficar el potencial con un área sombreada para mejor look
plt.plot(v_data['x'], v_data['V'], color='#2c3e50', linewidth=3, label='Potencial $V(x)$')
plt.fill_between(v_data['x'], v_data['V'], v_data['V'].max(), color='#2c3e50', alpha=0.1)

# Graficar niveles de energía como líneas horizontales largas
for i, row in e_levels.iterrows():
    E = row['E']
    line = plt.axhline(y=E, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    plt.text(v_data['x'].min() + 0.2, E + 0.2, f'E{i}={E:.2f}', color='red', fontweight='bold')

# Ajustes de estilo
plt.title('Comparación de Niveles de Energía Cuántica', fontsize=15)
plt.xlabel('Posición (x)', fontsize=12)
plt.ylabel('Energía', fontsize=12)
plt.grid(True, which='both', linestyle=':', alpha=0.5)
plt.axvline(0, color='black', linewidth=1)

# Limitar vista para que se aprecie el fondo
plt.ylim(v_data['V'].min() - 2, e_levels['E'].max() + 5)
plt.legend(loc='upper right')

plt.tight_layout()
plt.show()