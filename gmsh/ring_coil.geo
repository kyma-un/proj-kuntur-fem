// ============================================================
// Una bobina: anillo filamentario en z = 0
// Malla volumétrica de aire + curva embebida como guía de refinamiento
// ============================================================
SetFactory("OpenCASCADE");

// --- Parámetros ---------------------------------------------
R     = 0.100;      // Radio del anillo [m]
R_dom = 3.0 * R;    // Radio del dominio de aire

lc_coil = R / 30;   // tamaño sobre el anillo
lc_far  = R_dom / 4;

// --- Dominio de aire (esfera) -------------------------------
Sphere(1) = {0, 0, 0, R_dom};

// --- Bobina: un círculo en el plano z = 0 -------------------
Circle(100) = {0, 0, 0, R, 0, 2*Pi};
Curve{100} In Volume{1};

// --- Grupos físicos -----------------------------------------
Physical Volume("air")       = {1};
Physical Curve("coil")       = {100};
Physical Surface("boundary") = {1};

// --- Control de tamaño de malla -----------------------------
Mesh.Algorithm3D = 10;   // HXT

Field[1] = Distance;
Field[1].CurvesList = {100};
Field[1].Sampling   = 200;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_coil;
Field[2].SizeMax  = lc_far;
Field[2].DistMin  = R / 12;
Field[2].DistMax  = R / 2;

Field[3] = Ball;
Field[3].Radius    = R / 3;
Field[3].Thickness = R / 4;
Field[3].VIn       = R / 20;
Field[3].VOut      = lc_far;
Field[3].XCenter = 0; Field[3].YCenter = 0; Field[3].ZCenter = 0;

Field[10] = Min;
Field[10].FieldsList = {2, 3};
Background Field = 10;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh 3;
Save "ring_coil.msh";
