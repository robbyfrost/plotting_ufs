# --------------------------------
# Name: create_gif.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 10 July 2023
# Purpose: Create gifs using figures
# output from plotting scripts
# --------------------------------
from PIL import Image

image_path = []

for hr in range(7,19):
    # File path of the images to combine
    path = f'/scratch2/BMC/fv3lam/Robby.Frost/figures/2023041912/refc/refc_sidebyside_f{hr}_OK.png'
    # append to list
    image_path.append(path)

# Output GIF file path
output_path = '/scratch2/BMC/fv3lam/Robby.Frost/figures/2023041912/refc/refc_sidebyside_OK.gif'

# Create a list to store the image objects
images = []

# Open and append each image to the list
for file in image_path:
    image = Image.open(file)
    images.append(image)

# Save the images as an animated GIF
images[0].save(output_path, save_all=True, append_images=images[1:], optimize=False, duration=1000, loop=0)

print(f"GIF saved successfully at {output_path}")
