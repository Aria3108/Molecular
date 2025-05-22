import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np

# --- Funciones de energía ---
def energy1(a, b):
    pair = a + b
    if pair in ["CG", "GC"]:
        return -5
    elif pair in ["AU", "UA"]:
        return -4
    elif pair in ["GU", "UG"]:
        return -1
    else:
        return 0

def energy2(a, b):
    pair = a + b
    if pair in ["CG", "GC", "AU", "UA"]:
        return -1
    else:
        return 0

# --- Nussinov Algorithm ---
def nussinov(seq, energy_fn, min_loop_length=1):
    n = len(seq)
    dp = [[0]*n for _ in range(n)]
    bt = [[None]*n for _ in range(n)]

    for l in range(1, n):
        for i in range(n - l):
            j = i + l
            if j - i < min_loop_length:
                continue

            best = dp[i+1][j]
            bt[i][j] = (i+1, j)

            if dp[i][j-1] > best:
                best = dp[i][j-1]
                bt[i][j] = (i, j-1)

            if energy_fn(seq[i], seq[j]) < 0 and dp[i+1][j-1] + 1 > best:
                best = dp[i+1][j-1] + 1
                bt[i][j] = (i+1, j-1, True)

            for k in range(i+1, j):
                if dp[i][k] + dp[k+1][j] > best:
                    best = dp[i][k] + dp[k+1][j]
                    bt[i][j] = (i, k, k+1, j)

            dp[i][j] = best

    return dp, bt

# --- Backtracking ---
def get_pairs(i, j, bt, seq):
    pairs = []

    if i >= j:
        return pairs

    move = bt[i][j]
    if move is None:
        return pairs

    if len(move) == 2:
        return get_pairs(move[0], move[1], bt, seq)

    if len(move) == 3 and move[2]:
        pairs.append((i, j))
        pairs += get_pairs(move[0], move[1], bt, seq)
    elif len(move) == 4:
        pairs += get_pairs(move[0], move[1], bt, seq)
        pairs += get_pairs(move[2], move[3], bt, seq)

    return pairs

# --- Visualización en notación de paréntesis ---
def dot_bracket(seq, pairs):
    n = len(seq)
    structure = ['.'] * n
    for i, j in pairs:
        structure[i] = '('
        structure[j] = ')'
    return ''.join(structure)

# --- Visualización como grafo ---
def draw_structure(seq, pairs):
    n = len(seq)
    x = np.arange(n)
    y = np.zeros(n)

    fig, ax = plt.subplots(figsize=(n, 2))

    # Dibujar los nucleótidos como puntos con colores
    nucleotide_colors = {
    'A': '#e74c3c',  # rojo
    'U': '#e67e22',  # naranja
    'G': '#3498db',  # azul
    'C': '#2ecc71',  # verde
    }

    for i, base in enumerate(seq):
        ax.plot(x[i], y[i], 'o', color= nucleotide_colors[base], markersize=12)
        ax.text(x[i], y[i] - 0.1, base, ha='center', va='top', fontsize=12, fontweight='bold')

    # Dibujar el backbone
    ax.plot(x, y, '-', color='gray', alpha=0.3)

    # Dibujar arcos de emparejamiento
    for i, j in pairs:
        xm = (x[i] + x[j]) / 2
        r = abs(x[j] - x[i]) / 2
        theta = np.linspace(0, np.pi, 100)
        arc_x = xm + r * np.cos(theta)
        arc_y = r * np.sin(theta)
        ax.plot(arc_x, arc_y, color='purple', linewidth=1.5)

    # Ajustar estilo del gráfico
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-1, n)
    ax.set_ylim(-0.5, n/2)
    ax.axis('off')

    # Leyenda
    patches = [mpatches.Patch(color=color, label=base) for base, color in nucleotide_colors.items()]
    plt.legend(handles=patches, loc='upper right', frameon=False, title="Bases")

    plt.title("Estructura Secundaria del ARN", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()

# --- Ejemplo de uso ---
if __name__ == "__main__":
    seq = "GGAAAUCC"
    print("Secuencia:", seq)

    for idx, energy_fn in enumerate([energy1, energy2], start=1):
        print(f"\n--- Usando función de energía {idx} ---")
        dp, bt = nussinov(seq, energy_fn)
        pairs = get_pairs(0, len(seq)-1, bt, seq)

        # Calcular energía total liberada
        energia_total = sum(energy_fn(seq[i], seq[j]) for i, j in pairs)

        print("Pares emparejados:", pairs)
        print("Notación de paréntesis:", dot_bracket(seq, pairs))
        print("Energía total liberada:", energia_total)

        draw_structure(seq, pairs)
