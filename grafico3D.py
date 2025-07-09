import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def punto_en_plano(p, a, b, c, d, tol=1e-1):
    # Distance point-plane: |a·x + b·y + c·z - d| / √(a² + b² + c²)
    num = abs(a * p[0] + b * p[1] + c * p[2] - d)
    denom = np.linalg.norm([a, b, c])
    return num / denom < tol

def proyectar_a_plano(puntos, normal):
    """Projects 3D points to a 2D plane to calculate the ConvexHull."""
    normal = normal / np.linalg.norm(normal)
    # Create orthonormal basis
    if abs(normal[0]) > 1e-3:
        v = np.array([-(normal[1] + normal[2]) / normal[0], 1, 1])
    else:
        v = np.array([1, -(normal[0] + normal[2]) / normal[1], 1])
    v = v / np.linalg.norm(v)
    w = np.cross(normal, v)
    v = np.cross(w, normal)
    v, w = v / np.linalg.norm(v), w / np.linalg.norm(w)

    return np.array([[np.dot(p, v), np.dot(p, w)] for p in puntos])

# Main execution block when the script is run
if __name__ == '__main__':
    # Read JSON file path from command line arguments
    if len(sys.argv) < 2:
        print("Usage: python grafico3D.py <path_to_json_file> [winId]", file=sys.stderr)
        sys.exit(1)

    ruta_json = sys.argv[1]
    # winId_str = sys.argv[2] if len(sys.argv) > 2 else None # The winId is received but not used for embedding here

    try:
        with open(ruta_json, 'r') as f:
            data = json.load(f)

        # Load vertices
        vertices_data = data["vertices"]
        vertices = np.array([[v["x"], v["y"], v["z"]] for v in vertices_data])

        # Optimal point
        optimo = data.get("optimo", {})
        opt = np.array([optimo["x"], optimo["y"], optimo["z"]])

        # Constraints
        restricciones = data["restricciones"]

        # Colors for regions
        colores = ['#ffff99', '#ff9999', '#99ffcc', '#99ccff', '#ffcc99', '#ccccff']

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Analyze each constraint to find vertices on that plane
        for i, restriccion in enumerate(restricciones):
            a = restriccion['a']
            b = restriccion['b']
            c = restriccion['c']
            d = restriccion['rhs']
            normal = np.array([a, b, c])

            # Find vertices that are on the plane
            puntos_plano = [p for p in vertices if punto_en_plano(p, a, b, c, d)]

            if len(puntos_plano) >= 3:
                puntos_plano = np.array(puntos_plano)
                # Project to 2D
                puntos_2d = proyectar_a_plano(puntos_plano, normal)

                try:
                    hull = ConvexHull(puntos_2d)
                    cara = puntos_plano[hull.vertices]
                    poly = Poly3DCollection([cara], alpha=0.6, facecolor=colores[i % len(colores)], edgecolor='black')
                    ax.add_collection3d(poly)

                    # Label in the center
                    centro = np.mean(cara, axis=0)
                    ax.text(*centro, f"S{i+1}", fontsize=10, color='red')

                except Exception as e:
                    print(f"Could not build face for constraint S{i+1}: {e}", file=sys.stderr)

        # Draw vertices
        ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c='red', s=50, label='Vertices')

        # Optimal point
        ax.scatter(*opt, c='blue', s=100, label='Optimal', zorder=10)

        # Aesthetics
        ax.set_xlabel("x1")
        ax.set_ylabel("x2")
        ax.set_zlabel("x3")
        ax.set_title("Feasible Region 3D")
        ax.legend()
        ax.grid(True)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

        plt.tight_layout()
        plt.show() # This will open a separate Matplotlib window

    except FileNotFoundError:
        print(f"Error: JSON file not found at {ruta_json}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {ruta_json}. Check file format.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred in Python script: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr) # Print full traceback for debugging
        sys.exit(1)