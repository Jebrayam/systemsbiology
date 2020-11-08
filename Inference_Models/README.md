# Inference Models
- **Mean Cell**. Assumps that all cell can be simulated from the same set of parameters. This model works better to find estimate variability that comes from measure noise.
- **Mean Curve Fit**. Simplify version of Mean Cell model. Performs a fitting to the mean observations profile to estimate parameters based on a least square error function as a loss function.
- **Moment Based**. Computes the moments of the system upon the second order moment to perform estimation. This model is meant to approximate determinitically a biological system which expresses mainly intrinsic variability.
