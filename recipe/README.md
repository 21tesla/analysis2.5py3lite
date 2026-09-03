# Conda-forge recipe — ccpnmr (starting point)

`meta.yaml` in this directory is a **starting-point** conda-forge recipe for the
Python 3.14 port of ccpnmr 2.5.2. It is **not** yet a complete, CI-ready package.

## To actually publish to conda-forge

1. Create a fork of [conda-forge](https://github.com/conda-forge/staged-recipes)
   and a branch; copy this `recipe/` directory into it.
2. Fill in `meta.yaml`:
   - a public `source.url` (the sdist tarball) and its `sha256`, or a pinned
     VCS tag + commit;
   - the `extra.recipe-maintainers` handle(s).
3. Open a PR from `staged-recipes/` into `recipes/ccpnmr/`.
4. Make sure the conda-forge build environment can satisfy the **C-extension**
   build deps: a C/C++ compiler, GL, and Tk.

## Why this is optional

The canonical, primary distribution for this project is the **PyPI** sdist/wheel,
which already builds the C extensions at install time and is the "others can use
it" path. The conda recipe exists so the same build can be offered in the conda
world without diverging the packaging logic.
