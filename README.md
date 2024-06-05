# Knickpoint Propagation Model

These model codes simulate transient knickpoint propagation.

VectorExtractor: A function that takes in a TopoToolbox2 STREAMobj (i.e., georeferenced polyline) of a single stream channel and a TopoToolbox2 GRIDobj (i.e., georeferenced raster) of flow accumulation and outputs Matlab arrays of downstream length and upstream area values along the stream path defined by the input STREAMobj.

KPPropagation: A function that takes in the output arrays of VectorExtractor, plus user defined model parameters and a capture age, and outputs numerical values of the total knickpoint propagation distance and max, min, and mean knickpoint celerity (i.e., propagation velocity).

RunKPModel_OptimKT: A function that takes in a numeric array of observed knickpoint propagation distances (for multiple observed knickpoints), an array of STREAMobj corresponding to channels with observed knickpoints, a GRIDobj of flow accumulation, and user defined model parameters representing a sample space of possible parameter values, and runs KPPropagation to produce a arrays of modeled knickpoint propagation distances, optimal model parameters, and model misfits for all streams with observed knickpoints.

GenerateKPModelOutputFiles: A function that takes the output of RunKPModel_OptimKT and outputs a formatted .csv file with attributes of modeled knickpoints that can be imported into GIS software.
