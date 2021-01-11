An R implementation of the palette extraction part of Tam, Echevarria and
Gingold's paper ”[Efficient palette-based decomposition and recoloring of
images via RGBXY-space geometry](https://cragl.cs.gmu.edu/fastlayers/)”
(2018).

Needs `qhull` and `glpk` installed locally.

### Examples
Extracting palettes using k-means often misses perceptually important colors:

![k-means](./images/kmeans.png)

Whereas `colorhull` includes these but might also produce additional colors
without much perceptual importance:

![colorhull](./images/colorhull.png)

(Example photo by author)
