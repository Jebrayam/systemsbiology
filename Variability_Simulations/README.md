# Simulation of different types of varibility used to simulate gene expression in cell populations.
Each variability source is based on different assumptions.

- **Measurement Noise**: Takes in account the noise that comes from the measure devices. Each species is affected differently accoirding to the way is measured. But mainly the noise is represented by two terms: aditive and multiplicative. Also, over the whole distribution of measurements is applied a Gaussian noise, that it is multiplied by the measurements.

- **Intrinsic Variability**: This variability source comes from the stochasticity inside the cell as a result of the randomness in the reactions that take place and the instant of time when occur.

- **Extrinsic Variability**: Represents all the fenotipic differences among a cell population as parameter distribution. So, that each cell has a unique set of parameters which are used to simulate its expression.
