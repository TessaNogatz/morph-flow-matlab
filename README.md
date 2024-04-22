# morph-flow-matlab

Displacement estimation based on the morphological coarse-to-fine scheme

### Usage

First add all subfolders to the Matlab search path. From download/clone root folder, call 
```matlab
addpath(genpath(pwd))
```
Note that also the morphological wavelet transform needs to be available in the search path. It can be found [here](https://github.com/TessaNogatz/morphological-wavelets-matlab)
```matlab
addpath(genpath('path\to\morphological-wavelets-matlab'));
```

Then, calculate the morphological wavelet transform in moving and fixed image for a pre-defined number of levels (must be a multiple of 3, here e.g. 6)
```matlab

[C0, S0] = QLiftDec3MinMaxMin(fixed,para.levels);
[C1, S1] = QLiftDec3MinMaxMin(moving,para.levels);
```

Lastly, calculate the displacement based on the decomposition. A set of example parameters can be found in ```run_morph.m```

```matlab    
[u, v, w] = morph_tvl1_of(C0, S0, C1, S1, para);
```

### License

This code is distributed under BSD-2-Clause license. If you want to usemorph-flow in a scientific context, please cite
<ul>
 <li> Nogatz, T., Redenbach, C., & Schladitz, K. (2023). MorphFlow: Estimating Motion in In Situ Tests of Concrete. arXiv preprint arXiv:2310.11109. </li>
</ul>
