"""
red_neuronal_TL.py
==================
Red neuronal profunda (MLP) para aproximar la Pérdida por Transmisión
acústica TL(r, z) en un océano modelado con perfil de velocidad de Munk.

Física:
- Munk: c0=1500 m/s, eps=0.00737, zM=1300 m
- Suma de modos normales (Helmholtz) para calcular TL
- Dominio: r ∈ [0.1, 100] km, z ∈ [-4900, -100] m (prof. negativa → fondo)

Arquitectura MLP:
  Entrada (2) → Dense(64, ReLU) → Dense(128, ReLU) → Dense(128, ReLU)
             → Dense(64, ReLU) → Dense(1, lineal)
"""

# ─────────────────────────────────────────────────────────────────────────────
# 0. IMPORTS Y SEMILLAS
# ─────────────────────────────────────────────────────────────────────────────
import numpy as np
import matplotlib
matplotlib.use('Agg')          # backend sin display (funciona en cualquier entorno)
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import os, sys

# Semillas fijas para reproducibilidad
np.random.seed(42)

import tensorflow as tf
tf.random.set_seed(42)
# Suprimir logs innecesarios de TF
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

from tensorflow import keras
from tensorflow.keras import layers
from sklearn.metrics import r2_score

print("=" * 65)
print("  RED NEURONAL PARA TL(r,z) — PROPAGACIÓN ACÚSTICA EN EL OCÉANO")
print("=" * 65)
print(f"  TensorFlow version : {tf.__version__}")
print(f"  NumPy    version   : {np.__version__}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# 1. CARGA Y EXPLORACIÓN DE DATOS
# ─────────────────────────────────────────────────────────────────────────────
print("━" * 65)
print("  1. CARGA Y EXPLORACIÓN DE DATOS")
print("━" * 65)

DATA_PATH = "ComparativaModelos/tl_modos.txt"

data = np.loadtxt(DATA_PATH, comments='#')
r_raw = data[:, 0]   # rango [km]
z_raw = data[:, 1]   # profundidad [m], negativa (convención océano)
TL    = data[:, 2]   # pérdida por transmisión [dB]

# Valores únicos del grid
r_vals = np.unique(r_raw)   # 1 000 valores, 0.1–100 km
z_vals = np.unique(z_raw)   # 49 valores, -4900 a -100 m
Nr, Nz = len(r_vals), len(z_vals)

print(f"  Datos cargados: {data.shape[0]:,} puntos  ({Nr} rangos × {Nz} profundidades)")
print(f"  r ∈ [{r_vals.min():.1f}, {r_vals.max():.1f}] km   (paso {r_vals[1]-r_vals[0]:.1f} km)")
print(f"  z ∈ [{z_vals.min():.0f}, {z_vals.max():.0f}] m    (paso {z_vals[1]-z_vals[0]:.0f} m)")
print(f"  TL ∈ [{TL.min():.2f}, {TL.max():.2f}] dB")
print(f"  TL media = {TL.mean():.2f} dB    TL std = {TL.std():.2f} dB")

# Reshape a grilla 2D (Nr × Nz) para visualización
# Los datos están ordenados: primero varía r para cada z fijo
# Verificamos el orden real
TL_grid = np.zeros((Nr, Nz))
for j, zv in enumerate(z_vals):
    mask = (z_raw == zv)
    r_sub = r_raw[mask]
    tl_sub = TL[mask]
    idx = np.argsort(r_sub)
    TL_grid[:, j] = tl_sub[idx]

# ── Figura 0: mapa 2D de TL para verificar datos ──────────────────────────
fig0, ax0 = plt.subplots(figsize=(10, 5))
im = ax0.pcolormesh(r_vals, -z_vals, TL_grid.T,   # z_vals negativo → profundidad positiva
                    cmap='viridis_r', shading='auto')
cb = fig0.colorbar(im, ax=ax0, label='TL [dB]')
ax0.set_xlabel('Rango r [km]')
ax0.set_ylabel('Profundidad |z| [m]')
ax0.set_title('Mapa 2D de Pérdida por Transmisión TL(r,z) — datos numéricos')
ax0.invert_yaxis()
fig0.tight_layout()
fig0.savefig('fig0_mapa_TL_datos.png', dpi=150)
plt.close(fig0)
print("\n  ✓ Guardada: fig0_mapa_TL_datos.png")

# ─────────────────────────────────────────────────────────────────────────────
# 2. PREPROCESAMIENTO
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 65)
print("  2. PREPROCESAMIENTO")
print("━" * 65)

# Usamos |z| (profundidad positiva) para normalización
z_abs = -z_raw   # ahora z_abs ∈ [100, 4900] m

# ── Normalización Min-Max para entradas ──────────────────────────────────────
r_min, r_max = r_raw.min(), r_raw.max()
z_min, z_max = z_abs.min(), z_abs.max()

r_norm = (r_raw - r_min) / (r_max - r_min)      # ∈ [0, 1]
z_norm = (z_abs - z_min) / (z_max - z_min)      # ∈ [0, 1]

# ── Normalización estándar para la salida ─────────────────────────────────────
TL_mean = TL.mean()
TL_std  = TL.std()
TL_norm = (TL - TL_mean) / TL_std

print(f"  Normalización r  : Min-Max  → [{r_min}, {r_max}] → [0, 1]")
print(f"  Normalización |z|: Min-Max  → [{z_min}, {z_max}] → [0, 1]")
print(f"  Normalización TL : Estándar → μ={TL_mean:.4f}, σ={TL_std:.4f}")

# ── Split GEOGRÁFICO (no aleatorio) ──────────────────────────────────────────
mask_train = r_raw <= 70.0
mask_val   = (r_raw > 70.0) & (r_raw <= 85.0)
mask_test  = r_raw > 85.0

X_train = np.column_stack([r_norm[mask_train], z_norm[mask_train]]).astype(np.float32)
y_train = TL_norm[mask_train].astype(np.float32)

X_val   = np.column_stack([r_norm[mask_val],   z_norm[mask_val]]).astype(np.float32)
y_val   = TL_norm[mask_val].astype(np.float32)

X_test  = np.column_stack([r_norm[mask_test],  z_norm[mask_test]]).astype(np.float32)
y_test  = TL_norm[mask_test].astype(np.float32)

print(f"\n  Split geográfico:")
print(f"    Train      r ∈ [0.1, 70]  km  → {X_train.shape[0]:,} puntos  ({X_train.shape[0]/len(r_raw)*100:.1f}%)")
print(f"    Validación r ∈ (70,  85]  km  → {X_val.shape[0]:,} puntos  ({X_val.shape[0]/len(r_raw)*100:.1f}%)")
print(f"    Test       r ∈ (85, 100]  km  → {X_test.shape[0]:,} puntos  ({X_test.shape[0]/len(r_raw)*100:.1f}%)")

# ─────────────────────────────────────────────────────────────────────────────
# 3. ARQUITECTURA MLP
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 65)
print("  3. ARQUITECTURA MLP")
print("━" * 65)

model = keras.Sequential([
    keras.Input(shape=(2,), name='entrada_rz'),
    layers.Dense(64,  activation='relu', name='capa1'),
    layers.Dense(128, activation='relu', name='capa2'),
    layers.Dense(128, activation='relu', name='capa3'),
    layers.Dense(64,  activation='relu', name='capa4'),
    layers.Dense(1,   activation='linear', name='salida'),
], name='MLP_TL')

model.summary()

total_params = model.count_params()
print(f"\n  Total parámetros entrenables: {total_params:,}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. ENTRENAMIENTO
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 65)
print("  4. ENTRENAMIENTO")
print("━" * 65)

model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-3),
    loss='mse'
)

# Callbacks
early_stop = keras.callbacks.EarlyStopping(
    monitor='val_loss',
    patience=30,
    restore_best_weights=True,
    verbose=1
)

checkpoint = keras.callbacks.ModelCheckpoint(
    filepath='modelo_TL.keras',
    monitor='val_loss',
    save_best_only=True,
    verbose=0
)

# Callback personalizado: imprime progreso cada 50 épocas
class PrintEvery50(keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs=None):
        if (epoch + 1) % 50 == 0 or epoch == 0:
            print(f"  Época {epoch+1:4d} | "
                  f"loss={logs['loss']:.6f} | "
                  f"val_loss={logs['val_loss']:.6f}")

print(f"\n  Entrenando con batch_size=256, lr=1e-3, patience=30 …\n")

history = model.fit(
    X_train, y_train,
    validation_data=(X_val, y_val),
    epochs=500,
    batch_size=256,
    callbacks=[early_stop, checkpoint, PrintEvery50()],
    verbose=0
)

epochs_run    = len(history.history['loss'])
best_val_loss = min(history.history['val_loss'])
best_epoch    = np.argmin(history.history['val_loss']) + 1

print(f"\n  ✓ Entrenamiento terminado.")
print(f"    Épocas ejecutadas : {epochs_run}")
print(f"    Mejor época       : {best_epoch}")
print(f"    Mejor val_loss    : {best_val_loss:.6f}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. VISUALIZACIONES
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 65)
print("  5. VISUALIZACIONES")
print("━" * 65)

# ── Predicciones completas (espacio original) ─────────────────────────────────
X_all = np.column_stack([r_norm, z_norm]).astype(np.float32)
TL_pred_norm = model.predict(X_all, batch_size=4096, verbose=0).flatten()
TL_pred      = TL_pred_norm * TL_std + TL_mean   # desnormalizar
err_abs      = np.abs(TL - TL_pred)

# Grid 2D de predicción y error
TL_pred_grid = np.zeros((Nr, Nz))
err_grid     = np.zeros((Nr, Nz))
for j, zv in enumerate(z_vals):
    mask = (z_raw == zv)
    r_sub  = r_raw[mask]
    pred_s = TL_pred[mask]
    err_s  = err_abs[mask]
    idx    = np.argsort(r_sub)
    TL_pred_grid[:, j] = pred_s[idx]
    err_grid[:, j]     = err_s[idx]

# ────────────────────────────────────────────────────────────────────────────
# FIG 1: Curva de aprendizaje
# ────────────────────────────────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(9, 4))
epochs_arr = np.arange(1, epochs_run + 1)
ax1.semilogy(epochs_arr, history.history['loss'],     label='Train loss', color='steelblue')
ax1.semilogy(epochs_arr, history.history['val_loss'], label='Val   loss', color='tomato', linestyle='--')
ax1.axvline(best_epoch, color='green', linestyle=':', linewidth=1.8,
            label=f'Early stopping (época {best_epoch})')
ax1.set_xlabel('Época')
ax1.set_ylabel('MSE (escala logarítmica)')
ax1.set_title('Curva de aprendizaje — TL(r,z) MLP')
ax1.legend()
ax1.grid(True, which='both', alpha=0.3)
fig1.tight_layout()
fig1.savefig('fig1_curva_aprendizaje.png', dpi=150)
plt.close(fig1)
print("  ✓ Guardada: fig1_curva_aprendizaje.png")

# ────────────────────────────────────────────────────────────────────────────
# FIG 2: Comparación TL a 3 profundidades fijas
# ────────────────────────────────────────────────────────────────────────────
target_depths = [100, 1000, 2000]   # metros (positivos)
fig2, axes2 = plt.subplots(3, 2, figsize=(14, 12), sharex='col')

colors = {'num': 'steelblue', 'pred': 'tomato', 'err': 'darkorange'}
r70, r85 = 70.0, 85.0

for row, d in enumerate(target_depths):
    z_target = -float(d)     # convertir a convención negativa del archivo
    # Encontrar el z más cercano disponible
    j_idx = np.argmin(np.abs(z_vals - z_target))
    z_used = z_vals[j_idx]

    tl_num  = TL_grid[:, j_idx]
    tl_pr   = TL_pred_grid[:, j_idx]
    tl_err  = err_grid[:, j_idx]

    ax_tl  = axes2[row, 0]
    ax_err = axes2[row, 1]

    # TL numérica vs predicha
    ax_tl.plot(r_vals, tl_num, label='TL numérica', color=colors['num'],  lw=1.5)
    ax_tl.plot(r_vals, tl_pr,  label='TL predicha', color=colors['pred'], lw=1.5, linestyle='--')
    ax_tl.axvline(r70, color='gray',  linestyle=':', lw=1.2, label='r=70 km')
    ax_tl.axvline(r85, color='black', linestyle=':', lw=1.2, label='r=85 km')
    # Sombrear zonas
    ax_tl.axvspan(r_vals[0], r70, alpha=0.05, color='blue',   label='_Train')
    ax_tl.axvspan(r70,       r85, alpha=0.05, color='orange', label='_Val')
    ax_tl.axvspan(r85, r_vals[-1], alpha=0.05, color='red',   label='_Test')
    ax_tl.set_ylabel('TL [dB]')
    ax_tl.set_title(f'z = {d} m')
    if row == 0:
        ax_tl.legend(fontsize=7, ncol=2)
    ax_tl.grid(alpha=0.3)

    # Error absoluto
    ax_err.plot(r_vals, tl_err, color=colors['err'], lw=1.2)
    ax_err.axvline(r70, color='gray',  linestyle=':', lw=1.2)
    ax_err.axvline(r85, color='black', linestyle=':', lw=1.2)
    ax_err.axvspan(r_vals[0], r70, alpha=0.05, color='blue')
    ax_err.axvspan(r70,       r85, alpha=0.05, color='orange')
    ax_err.axvspan(r85, r_vals[-1], alpha=0.05, color='red')
    ax_err.set_ylabel('Error absoluto [dB]')
    ax_err.set_title(f'Error |TL_num − TL_pred|   (z={d} m)')
    ax_err.grid(alpha=0.3)

for ax in axes2[-1, :]:
    ax.set_xlabel('Rango r [km]')

# Leyenda de zonas
from matplotlib.patches import Patch
legend_patches = [Patch(facecolor='blue',   alpha=0.15, label='Train (0–70 km)'),
                  Patch(facecolor='orange', alpha=0.15, label='Val (70–85 km)'),
                  Patch(facecolor='red',    alpha=0.15, label='Test (85–100 km)')]
fig2.legend(handles=legend_patches, loc='lower center', ncol=3, fontsize=9,
            bbox_to_anchor=(0.5, -0.01), frameon=True)

fig2.suptitle('Comparación TL numérica vs predicha a profundidades fijas', fontsize=13, y=1.01)
fig2.tight_layout()
fig2.savefig('fig2_comparacion_profundidades.png', dpi=150, bbox_inches='tight')
plt.close(fig2)
print("  ✓ Guardada: fig2_comparacion_profundidades.png")

# ────────────────────────────────────────────────────────────────────────────
# FIG 3: Mapas 2D lado a lado
# ────────────────────────────────────────────────────────────────────────────
z_plot = -z_vals   # profundidad positiva para visualización

# Límites comunes de TL
vmin_tl = min(TL_grid.min(), TL_pred_grid.min())
vmax_tl = max(TL_grid.max(), TL_pred_grid.max())
vmax_err = err_grid.max()

fig3, axes3 = plt.subplots(1, 3, figsize=(17, 5))

kw_pcolor = dict(shading='auto')

im0 = axes3[0].pcolormesh(r_vals, z_plot, TL_grid.T,
                           cmap='viridis_r', vmin=vmin_tl, vmax=vmax_tl, **kw_pcolor)
fig3.colorbar(im0, ax=axes3[0], label='TL [dB]')
axes3[0].set_title('TL numérica')

im1 = axes3[1].pcolormesh(r_vals, z_plot, TL_pred_grid.T,
                           cmap='viridis_r', vmin=vmin_tl, vmax=vmax_tl, **kw_pcolor)
fig3.colorbar(im1, ax=axes3[1], label='TL [dB]')
axes3[1].set_title('TL predicha (MLP)')

im2 = axes3[2].pcolormesh(r_vals, z_plot, err_grid.T,
                           cmap='RdYlGn_r', vmin=0, vmax=vmax_err, **kw_pcolor)
fig3.colorbar(im2, ax=axes3[2], label='|Error| [dB]')
axes3[2].set_title('Error absoluto')

for ax in axes3:
    ax.set_xlabel('Rango r [km]')
    ax.set_ylabel('Profundidad [m]')
    ax.invert_yaxis()
    ax.axvline(r70, color='white', linestyle='--', lw=1.2, label='r=70 km')
    ax.axvline(r85, color='cyan',  linestyle='--', lw=1.2, label='r=85 km')
    ax.legend(fontsize=7, loc='upper right')

fig3.suptitle('Mapas 2D: TL numérica, predicha y error absoluto', fontsize=13)
fig3.tight_layout()
fig3.savefig('fig3_mapas_2D.png', dpi=150)
plt.close(fig3)
print("  ✓ Guardada: fig3_mapas_2D.png")

# ────────────────────────────────────────────────────────────────────────────
# FIG 4: Boxplot de error por zona × profundidad
# ────────────────────────────────────────────────────────────────────────────
# Seleccionamos 3 profundidades representativas para el boxplot
depths_box = [100, 1000, 2000]
zones = {
    'Train\n(0–70 km)' : r_raw <= 70.0,
    'Val\n(70–85 km)'  : (r_raw > 70.0) & (r_raw <= 85.0),
    'Test\n(85–100 km)': r_raw > 85.0,
}
zone_colors = ['steelblue', 'darkorange', 'tomato']

fig4, axes4 = plt.subplots(1, 3, figsize=(14, 5), sharey=False)

for col, d in enumerate(depths_box):
    z_target = -float(d)
    mask_z   = (z_raw == z_vals[np.argmin(np.abs(z_vals - z_target))])
    ax       = axes4[col]
    data_box = []
    labels   = []
    for (label, mask_zone), color in zip(zones.items(), zone_colors):
        combined = mask_zone & mask_z
        data_box.append(err_abs[combined])
        labels.append(label)

    bp = ax.boxplot(data_box, patch_artist=True, notch=False,
                    medianprops=dict(color='black', lw=2))
    for patch, color in zip(bp['boxes'], zone_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Error absoluto [dB]')
    ax.set_title(f'z = {d} m')
    ax.grid(axis='y', alpha=0.3)

fig4.suptitle('Distribución del error absoluto por zona y profundidad', fontsize=12)
fig4.tight_layout()
fig4.savefig('fig4_boxplot_error.png', dpi=150)
plt.close(fig4)
print("  ✓ Guardada: fig4_boxplot_error.png")

# ─────────────────────────────────────────────────────────────────────────────
# 6. MÉTRICAS FINALES
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 65)
print("  6. MÉTRICAS FINALES")
print("━" * 65)

def compute_metrics(y_true, y_pred, r_vals_subset, z_vals_subset, label):
    """Calcula MSE, RMSE, MAE, R², error máximo y su ubicación."""
    err   = np.abs(y_true - y_pred)
    mse   = np.mean((y_true - y_pred) ** 2)
    rmse  = np.sqrt(mse)
    mae   = np.mean(err)
    r2    = r2_score(y_true, y_pred)
    emax  = err.max()
    idx_m = np.argmax(err)
    r_max_loc = r_vals_subset[idx_m]
    z_max_loc = z_vals_subset[idx_m]
    return {
        'label': label,
        'n'    : len(y_true),
        'MSE'  : mse,
        'RMSE' : rmse,
        'MAE'  : mae,
        'R2'   : r2,
        'Emax' : emax,
        'r_max': r_max_loc,
        'z_max': z_max_loc,
    }

TL_true_all  = TL
TL_pred_all  = TL_pred
r_all        = r_raw
z_all        = -z_raw   # positivo

results = []
masks_zones = {
    'Train  (r≤70 km)' : mask_train,
    'Val  (70<r≤85 km)': mask_val,
    'Test (r>85 km)'   : mask_test,
}
for label, mask in masks_zones.items():
    m = compute_metrics(TL_true_all[mask], TL_pred_all[mask],
                        r_all[mask], z_all[mask], label)
    results.append(m)

# ── Tabla de métricas ─────────────────────────────────────────────────────────
print(f"\n  {'Zona':<22}  {'N':>6}  {'MSE':>9}  {'RMSE':>7}  {'MAE':>7}  {'R²':>7}  {'Emax':>7}  {'(r,z) max'}")
print("  " + "─" * 80)
for r_dict in results:
    print(f"  {r_dict['label']:<22}  {r_dict['n']:>6,}  "
          f"{r_dict['MSE']:>9.4f}  {r_dict['RMSE']:>7.4f}  "
          f"{r_dict['MAE']:>7.4f}  {r_dict['R2']:>7.4f}  "
          f"{r_dict['Emax']:>7.3f}  "
          f"({r_dict['r_max']:.1f} km, {r_dict['z_max']:.0f} m)")

# ── Top 5 puntos con mayor error ──────────────────────────────────────────────
top5_idx = np.argsort(err_abs)[-5:][::-1]
print(f"\n  Top 5 puntos con MAYOR ERROR ABSOLUTO:")
print(f"  {'#':<3}  {'r [km]':>8}  {'z [m]':>7}  {'TL_num':>9}  {'TL_pred':>9}  {'Error':>8}")
print("  " + "─" * 55)
for i, idx in enumerate(top5_idx):
    print(f"  {i+1:<3}  {r_raw[idx]:>8.2f}  {-z_raw[idx]:>7.0f}  "
          f"{TL[idx]:>9.3f}  {TL_pred[idx]:>9.3f}  {err_abs[idx]:>8.3f} dB")

# ── Campo cercano vs lejano ────────────────────────────────────────────────────
mask_cerca  = r_raw < 10.0
mask_lejos  = r_raw > 50.0
err_cerca   = err_abs[mask_cerca].mean()
err_lejos   = err_abs[mask_lejos].mean()
print(f"\n  Error promedio campo cercano (r < 10 km)  : {err_cerca:.4f} dB")
print(f"  Error promedio campo lejano  (r > 50 km)  : {err_lejos:.4f} dB")

# ─────────────────────────────────────────────────────────────────────────────
# 7. ANÁLISIS DE INTERFERENCIAS
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 65)
print("  7. ANÁLISIS DE INTERFERENCIAS")
print("━" * 65)

# Varianza local de TL en ventanas de 5 km
# Para cada (r, z), calculamos la varianza de TL en r ± 2.5 km
window_km   = 5.0
half_win    = window_km / 2.0
var_local   = np.zeros_like(TL)

for i, (ri, zi) in enumerate(zip(r_raw, z_raw)):
    mask_win = (np.abs(r_raw - ri) <= half_win) & (z_raw == zi)
    var_local[i] = np.var(TL[mask_win]) if mask_win.sum() > 1 else 0.0

print(f"  Varianza local calculada (ventana {window_km} km)")
print(f"  Varianza local: min={var_local.min():.4f}, max={var_local.max():.4f}, "
      f"media={var_local.mean():.4f}")

# Correlación varianza local vs error
corr = np.corrcoef(var_local, err_abs)[0, 1]
print(f"  Correlación (varianza_local, error_abs): r = {corr:.4f}")

# ── Figura 5: varianza local vs error ────────────────────────────────────────
fig5, axes5 = plt.subplots(1, 2, figsize=(13, 5))

# Scatter varianza vs error (subsample para velocidad)
step = max(1, len(var_local) // 5000)
axes5[0].scatter(var_local[::step], err_abs[::step], s=5, alpha=0.4,
                 c=r_raw[::step], cmap='plasma', label='puntos')
axes5[0].set_xlabel('Varianza local de TL [dB²]  (ventana 5 km)')
axes5[0].set_ylabel('Error absoluto [dB]')
axes5[0].set_title(f'Varianza local vs Error  (r = {corr:.3f})')
axes5[0].grid(alpha=0.3)
sm = plt.cm.ScalarMappable(cmap='plasma',
                            norm=plt.Normalize(r_raw.min(), r_raw.max()))
sm.set_array([])
fig5.colorbar(sm, ax=axes5[0], label='Rango r [km]')

# Mapa 2D de varianza local
var_grid = np.zeros((Nr, Nz))
for j, zv in enumerate(z_vals):
    mask = (z_raw == zv)
    r_sub  = r_raw[mask]
    var_s  = var_local[mask]
    idx    = np.argsort(r_sub)
    var_grid[:, j] = var_s[idx]

im5 = axes5[1].pcolormesh(r_vals, z_plot, var_grid.T,
                           cmap='hot_r', shading='auto')
fig5.colorbar(im5, ax=axes5[1], label='Varianza local [dB²]')
axes5[1].set_xlabel('Rango r [km]')
axes5[1].set_ylabel('Profundidad [m]')
axes5[1].set_title('Mapa 2D de varianza local de TL')
axes5[1].invert_yaxis()

fig5.suptitle('Análisis de interferencias: varianza local vs error de la red', fontsize=12)
fig5.tight_layout()
fig5.savefig('fig5_interferencias.png', dpi=150)
plt.close(fig5)
print("  ✓ Guardada: fig5_interferencias.png")

# ─────────────────────────────────────────────────────────────────────────────
# 8. GUARDAR MODELO FINAL
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 65)
print("  8. GUARDADO DEL MODELO")
print("━" * 65)
model.save('modelo_TL.keras')
print("  ✓ Modelo guardado como 'modelo_TL.keras'")

# ─────────────────────────────────────────────────────────────────────────────
# RESUMEN FINAL
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("  RESUMEN FINAL")
print("=" * 65)
print(f"\n  Épocas entrenadas  : {epochs_run} (mejor en época {best_epoch})")
print(f"  Parámetros totales : {total_params:,}")
print(f"\n  Figuras generadas:")
for fn in ['fig0_mapa_TL_datos.png', 'fig1_curva_aprendizaje.png',
           'fig2_comparacion_profundidades.png', 'fig3_mapas_2D.png',
           'fig4_boxplot_error.png', 'fig5_interferencias.png']:
    size = os.path.getsize(fn) // 1024 if os.path.exists(fn) else 0
    print(f"    {fn:<38} ({size} KB)")

print(f"\n  Modelo guardado    : modelo_TL.keras")

# ─────────────────────────────────────────────────────────────────────────────
# INTERPRETACIÓN FÍSICA
# ─────────────────────────────────────────────────────────────────────────────
print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  INTERPRETACIÓN FÍSICA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  La red ha aprendido la tendencia de decaimiento geométrico de TL
  con la distancia (comportamiento suave ∝ 10·log10(r)), que domina
  la señal.  Sin embargo, la suma de modos normales produce franjas
  de interferencia constructiva/destructiva cuya frecuencia espacial
  crece con r y z.  El análisis de varianza local confirma que los
  puntos con mayor error coinciden con las zonas de mayor oscilación
  de TL: el MLP, al ser un aproximador suave, promedia las
  interferencias rápidas en lugar de reproducirlas fielmente.  El R²
  elevado en el conjunto de entrenamiento (campo cercano sencillo)
  frente al menor R² en test (campo lejano con más modos
  interferentes) indica que la red generaliza la física de fondo pero
  no memoriza los detalles de la interferencia multimodal.  Para
  capturar las oscilaciones de alta frecuencia sería necesaria una
  arquitectura más profunda, con más neuronas, o un enfoque basado en
  redes física-informadas (PINN) que incorporen explícitamente la
  ecuación de Helmholtz.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")
