# MXB324 Aquifer

## Vertices

Number | Description | Coordinates
------ | ----------- | -----------
1 | Bottom left corner | x = 0 && z = 0
2 | Bottom horizontal boundary | 0 < x < 500 && z = 0
3 | Bottom right corner | x = 500 && z = 0
4 | Left sandstone boundary | x = 0 && 0 < z < 30
5 | Pump | x = 450 && z = 10
6 | Sandstone evapotranspiration zone | 350 < x < 500 && 76 < z < 80
7 | Normal sandstone interior point | 0 < x < 500 && 0 < z < 30 || 350 < x < 500 && 30 < z < 80
8 | Right sandstone boundary | x = 500 && 0 < z < 76
9 | Right sandstone boundary evapotranspiration | x = 500 && 76 < z < 80
10 | Sandstone/Alluvium interface and left boundary | x = 0 && z = 30
11 | Sandstone/Alluvium interface | 0 < x < 50 && z = 30
12 | Sandstone/Alluvium/Confining interface | x = 50 && z = 30
13 | Sandstone/Confining horizontal interface | 50 < x < 350 && z = 30
14 | Sandstone/Confining corner | x = 350 && z = 30
15 | Alluvium left boundary | x = 0 && 30 < z < 78
16 | Alluvium left boundary with evapotranspiration | x = 0 && 78 < z < 80
17 | Pump | x = 100 && z = 50
18 | Evapotranspiration in alluvium | 0 < x < 350 && 78 < z < 80
19 | Normal alluvium interior point | 0 < x < 50 && 30 < z < 80 || 0 < x < 350 && 40 < z < 80
20 | Alluvium/Confining vertical interface | x = 50 && 30 < z < 40
21 | Normal confining interior point | 50 < x < 350 && 30 < z < 40
22 | Sandstone/Confining vertical interface | x = 350 && 30 < z < 40
23 | Confining/Alluvium corner point | x = 50 && z = 40
24 | Confining/Alluvium horizontal interface | 50 < x < 350 && z = 40
25 | Confining/Alluvium/Sandstone upper right corner | x = 350 && z = 40
26 | Alluvium/Sandstone vertical interface | x = 350 && 40 < z < 76
27 | Alluvium/Sandstone vertical interface with evapotranspiration | x = 350 && 76 < z < 78
28 | Alluvium/Sandstone vertical interface with evapotranspiration | x = 350 && 78 < z < 80
29 | Top left node | x = 0 && z = 80
30 | Alluvium/Rain boundary | 0 < x < 350 && z = 80
31 | Alluvium/Sandstone/Rain boundary | x = 350 && z = 80
32 | Sandstone/Rain boundary | 350 < x < 500 && z = 80
33 | Top right node | x = 500 && z = 80

