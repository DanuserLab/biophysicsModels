This repository is moved to https://github.com/XWangLab-TU/BiophysicsModels.git. For inquiries, please contact former lab member, Xinxin Wang: xinxin-wang@utulsa.edu.

This is a project of a membrane model based on a recent manuscript (MS): 
[**Remeshing flexible membranes under the control of free energy**](https://doi.org/10.1371/journal.pcbi.1010766), *PLOS Computational Biology*, 2022, 18(12): e1010766, written by Xinxin Wang, [Gaudenz Danuser](https://www.danuserlab-utsw.org/).

The membrane is modularized as a matlab object @ModMembrane, and external control points as @ModSubstrate.
@ModMembrane follows the physics-based remeshing algorithm in the MS to be able to represent generic morphologies, e.g. red blood cell.
@ModSubstrate is mechanically coupled to @ModMembrane via @TypForce for simulating the morphologies of various case studies in the MS.
All the modular objects are organized by the object @model before running the simulations.
All the available simulations are written in @BiophysicsApp, and can be run from 'Scripts/runMS2022Examples.m'.
Graphics are adjusted to the exact formats in the MS, and can be run from  'Scripts/runMS2022Graphics.m'.
 
Author: Xinxin Wang, Danuser Lab
email: wangxinxin8627@gmail.com
date: 2022/07/25

### Danuser Lab Links
[Danuser Lab Website](https://www.danuserlab-utsw.org/)

[Software Links](https://github.com/DanuserLab/)
