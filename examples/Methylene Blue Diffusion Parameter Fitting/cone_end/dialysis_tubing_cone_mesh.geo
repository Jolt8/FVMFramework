SetFactory("OpenCASCADE");
Merge "dialysis_tubing_model_cone.step";
Mesh.ScalingFactor = 0.001;
Coherence;
Physical Volume("dialysis_tubing_interior", 15) = {3, 2};
//+
Physical Volume("surrounding_fluid", 16) = {1};
//+
Physical Surface("dialysis_tubing_surface", 17) = {4, 5};
//+
Coherence;
Mesh 3;
Coherence Mesh;
Save "dialysis_tubing_cone_output.msh";
//+
