import os
from PIL import Image

def alternative_white_background(image_path, max_size_kb=500, target_dim=(600, 600)):
    if not os.path.exists(image_path):
        print(f"Hata: {image_path} bulunamadı.")
        return

    # Orijinal resmi aç
    img = Image.open(image_path).convert("RGB")
    
    # Fond clair -> Convertit le fond gris/jaune/vert proche du blanc en blanc pur
    # (Seuil de tolérance fixé à 150 pour blanchir le fond sans toucher aux cheveux sombres)
    pixels = img.load()
    for y in range(img.size[1]):
        for x in range(img.size[0]):
            r, g, b = pixels[x, y]
            # Si le pixel est globalement clair (fond), on le force en blanc pur
            if r > 150 and g > 150 and b > 130:
                pixels[x, y] = (255, 255, 255)

    directory, filename = os.path.split(image_path)
    name, ext = os.path.splitext(filename)
    output_path = os.path.join(directory, f"{name}_final_biometric{ext}")

    # Dimensionnement et compression
    img.thumbnail(target_dim, Image.Resampling.LANCZOS)
    quality = 95
    img.save(output_path, quality=quality, optimize=True)
    
    while os.path.getsize(output_path) > max_size_kb * 1024 and quality > 10:
        quality -= 5
        img.save(output_path, quality=quality, optimize=True)
        
    print(f"\nİşlem tamamlandı: {output_path}")

if __name__ == "__main__":
    target_image = r"C:\Downloads\New folder (4)\muhammed.JPG"
    alternative_white_background(target_image)
