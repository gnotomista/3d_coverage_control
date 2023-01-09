# 3d_coverage_control
Implementation of Lloyd's algorithm for multi-agent coverage control in 3D.

Videos of the simulations obtained by running the examples in the `src` folder:

* [10 agents moving in a cube with uniform density](https://github.com/gnotomista/3d_coverage_control/blob/main/src/coverage_3d_cube_uniform.m) (plotting agents, centroid and boundary of the Voronoi cells)

![example_1_cube_uniform_centroids_cells](https://user-images.githubusercontent.com/23127701/211418404-79ac2e4a-894f-42cf-897e-f7e66c01be47.gif)

* [10 agents moving in a dodecahedron with Gaussian density](https://github.com/gnotomista/3d_coverage_control/blob/main/src/coverage_3d_dodecahedron_gaussian.m) (plotting only agents)

![example_2_dodecahedron_gaussian](https://user-images.githubusercontent.com/23127701/211419165-36aaeb02-5148-4f1e-acfb-30d92ae307e3.gif)

---

Thanks to [Alessia Benevento](https://scholar.google.com/citations?user=us-LVKQAAAAJ&hl=en) for the first version of this code, based on the 2D implementation available in the [swarm_sim](https://github.com/gnotomista/swarm_sim) repository.

The code to compute bounded Voronoi cells in 3D is from [this](https://github.com/hyongju/Polytope-bounded-Voronoi-diagram) repository.

Additional utils from MATLAB Central File Exchange are used to draw [Platonic solids](https://www.mathworks.com/matlabcentral/fileexchange/77446-platonic-solids) and compute [permutations](https://www.mathworks.com/matlabcentral/fileexchange/7147-permn).
