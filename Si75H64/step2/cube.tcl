##################

mol load cube sample.cube
mol modstyle 0 0 CPK 1.000000 0.300000 12.000000 12.000000
mol representation CPK 1.000000 0.300000 12.000000 12.000000

mol addrep 0
mol modstyle 1 0 Isosurface 0.04 0 0 0 1 1
mol modmaterial 1 0 Opaque
mol modcolor 1 0 ColorID 0

mol addrep 0
mol modcolor 2 0 ColorID 1
mol modstyle 2 0 Isosurface -0.04 0 0 0 1 1

color Display Background white

#rotate x by 225.0
scale by 0.5
axes location off

render TachyonInternal sample.tga
