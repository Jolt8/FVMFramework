SetFactory("OpenCASCADE");
Merge "packed bed 2D sketch.step";
Merge "packed_bed_2D_mesher.opt";
Mesh.ScalingFactor = 0.001;
Coherence;
Mesh 2;
Extrude {0, 3, 0} {
  Surface{2}; Surface{1}; Surface{3}; Layers {1}; Recombine;
}
Mesh 3;//+
Physical Surface("inlet_surface", 41) = {5};
//+
Physical Surface("outlet_surface", 42) = {18};
//+
Physical Volume("inlet_volume", 43) = {1};
//+
Physical Volume("outlet_volume", 44) = {3};
//+
Physical Volume("reacting_volume", 45) = {2};
//+
Coherence Mesh;
Save "output.msh";
