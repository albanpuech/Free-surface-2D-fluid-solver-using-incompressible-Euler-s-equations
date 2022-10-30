# Project instructions :

Your second project grade will be based on labs 6-7-8, which should have resulted in a free-surface 2D fluid solver using incompressible Euler's equations. The report's repository should be sent to me by June ZZth, 23:59.

Notably you will need to do:

- A Voronoï diagram using Voronoï Parallel Linear Enumeration (Sec. 4.3.3, lab 6) with Sutherland-Hodgman polygon clipping algorithm (Sec. 4.2, lab 6). It is not required to accelerate the clipping using the kd-tree nearest neighbor search.
- Extending this to a Power diagram (Sec. 4.4.3, lab 7)
- Optimizing the weights of the power diagram using LBFGS (Sec 4.4.4, lab 7)
- Implement de Gallouet-Mérigot incompressible Euler scheme that prescribes each fluid cell area and the sum of air cell areas to given values, and add a spring force from each fluid particle to their Laguerre's cell centroid. (Sec. 5.4, lab 8).

# Compile :

```
g++ -Ofast -o project project.cpp lbfgs.c -fopenmp -I.  
```

# File description 

- report 2.pdf: Report with feature description
- project_2.cpp: Code of the free-surface 2D fluid solver
- gif_slow.gif/ gif_fast.gif:  fluid simulation videos
- frames.rar: frames outputted by the simulator
