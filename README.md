# Polynomial Approximation Optimizer

A high-precision polynomial approximation tool that uses linear programming to find optimal coefficients. This implementation specializes in handling high-degree polynomials with arbitrary precision arithmetic to minimize rounding errors.

## Features

- **Arbitrary Precision Arithmetic**: Uses `mpmath` library to perform calculations with user-defined precision, reducing rounding errors significantly
- **High-Degree Polynomials**: Capable of finding coefficients for high-order polynomial approximations
- **Weighted Approximation**: Supports custom weight functions to bound the approximation above and below
- **Visualization**: Includes plotting capabilities with LaTeX rendering for publication-quality figures

## Current Limitations

- Performance degrades significantly for discrete intervals with more than 30 points
- Best suited for problems requiring high precision with relatively few evaluation points
- Memory usage increases with both polynomial degree and number of evaluation points

## Installation

### Prerequisites
```bash
pip install mpmath matplotlib
```

## Usage

Here's a basic example of how to use the optimizer:

```python
import mpmath as mp
from optimal_poly import OptimalPolynomial

# Set up the function to approximate
def f(x):
    return mp.exp(x) * mp.cos(2 * mp.pi * x)

# Create optimizer with desired precision (e.g., 50 decimal places)
opt = OptimalPolynomial(precision=50)

# Find optimal coefficients
coefs = opt.get_coefs(
    f=f,          # Function to approximate
    degree=6,     # Polynomial degree
    a=-1,         # Left endpoint
    b=1,          # Right endpoint
    num=20        # Number of discretization points
)

# Print results
opt.print_status()

# Plot the approximation
opt.plot_fig(f, coefs)
```

### Using Weight Functions

You can provide weight functions to bound the approximation:

```python
# Define weight functions
def omega_sup(x):
    return mp.mpf('1.2')  # Upper bound: 20% above

def omega_inf(x):
    return mp.mpf('0.8')  # Lower bound: 20% below

# Find coefficients with bounds
coefs = opt.get_coefs(
    f=f,
    degree=6,
    a=-1,
    b=1,
    num=20,
    omega_sup=omega_sup,
    omega_inf=omega_inf
)
```

## Key Advantages

1. **High Precision**
   - Uses arbitrary-precision arithmetic to minimize rounding errors
   - Particularly useful for high-degree polynomials where standard floating-point arithmetic might fail

2. **Bounded Approximation**
   - Supports weight functions to control the approximation behavior
   - Allows for asymmetric bounds above and below the target function

3. **Flexible Usage**
   - Can be used for various types of functions
   - Works well with small discrete intervals requiring high precision
   - Suitable for applications where accuracy is more important than speed

## Best Practices

1. **Number of Points**
   - Keep the number of discretization points below 30 for reasonable performance
   - Use fewer points with higher precision for initial experiments

2. **Precision Setting**
   - Start with moderate precision (e.g., 50 decimal places) and adjust as needed
   - Higher precision is more important for higher-degree polynomials

3. **Weight Functions**
   - Use weight functions when you need to control the approximation bounds
   - Keep weight functions simple for better solver performance

## Technical Details

The implementation uses:
- Linear programming via the Revised Simplex Method
- Arbitrary precision arithmetic through `mpmath`
- Matrix operations optimized for arbitrary precision

## Known Issues

1. Performance significantly degrades with:
   - More than 30 discretization points

2. Memory usage can be high when:
   - Using many discretization points
   - Working with high-degree polynomials
   - Setting very high precision

## Contributing

Contributions are welcome! Areas that need improvement:
- Performance optimization for larger discrete intervals
- Memory usage optimization
- Additional visualization options
- More example cases and documentation

## License

This project is licensed under the MIT License - see the LICENSE file for details.