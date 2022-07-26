This is a project of a membrane model based on a recent manuscript (MS): 
"Remeshing flexible membranes under the control of free enegy",
by Xinxin Wang and Gaudenz Danuser.
The membrane is modularized as a matlab object @ModMembrane, and external control points as @ModSubstrate.
@ModMembrane follows the physics-based remeshing algorithm in the MS to be able to represent generic morphologies, e.g. red blood cell.
@ModSubstrate is mechanically coupled to @ModMembrane via and @TypForce for simulating the morphologies of various case studies in the MS.
All the modular objects are organized by the organizing object @model for running the simulations.
All the available simulations are written in @BiophysicsApp, and can be run from 'Scripts/runMS2022Examples.m'.
Graphics are adjusted to the exact formats in the MS via 'Scripts/runMS2022Graphics.m'.
 
Author: Xinxin Wang, Danuser Lab
email: wangxinxin8627@gmail.com
date: 2022/07/25

