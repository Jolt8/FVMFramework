SetFactory("OpenCASCADE");
Merge "dehi_model.step";
Mesh.ScalingFactor = 0.001;
Coherence;
//+
Physical Volume("implant_interior", 19) = {2};
//+
Physical Volume("surrounding_tissue", 20) = {1};
//+
Physical Surface("implant_surface", 21) = {8, 7, 9};
//+
Coherence;
Mesh 3;
Coherence Mesh;
Save "dehi_model_output.msh";
//+


