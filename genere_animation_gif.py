from PIL import Image
import glob


directory_images = 'dossier_images_7.7.26'
gif_file = "dykes_animation_github7.7.26_2.gif"

# --- Récupérer toutes les images ---
# Suppose que tes fichiers s’appellent "dyke_step_60.png", "dyke_step_100.png", etc.
frames = [Image.open(img) for img in sorted(glob.glob(f"{directory_images}/position_pas_*.png"))]

# --- Créer le GIF ---
frames[0].save(
    gif_file,
    save_all=True,
    append_images=frames[1:],
    duration=300,   # durée entre frames en ms (ici 300 ms = 0.3 s)
    loop=1          # 0 = boucle infinie
)