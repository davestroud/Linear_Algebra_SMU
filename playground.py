from skimage import data
camera = data.camera()
type(camera)
camera.shape
camera.size
camera.ndim
camera.min(), camera.max()
camera.mean()
camera[10, 20]
camera[3, 10] = 0
camera[3, 10]

import numpy as np
from skimage.util import img_as_float
image = np.arange(0, 50, 10, dtype=np.uint8)
print(image.astype(float))
print(img_as_float(image))