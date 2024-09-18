"""
This module computes symbolic expressions for coefficients of a 3D Hermite interpolation polynomial
using a 5-point derivative stencil with points (-3, -2, -1, 0, 1, 2, 3) in all dimensions.
It uses the derivative computations as specified in the provided Maple code.
"""

import os
import time
import multiprocessing
from itertools import product

import sympy as sp
from sympy.printing.pycode import PythonCodePrinter


def DATA(i, j, k):
    """Return a sympy Symbol for DATA at position (i, j, k)."""
    return sp.Symbol(f"DATA_{i}_{j}_{k}")


def create_coefficients():
    """Create symbolic coefficients for the interpolation polynomial."""
    coeffs = {}
    for i in range(4):
        for j in range(4):
            for k in range(4):
                coeffs[(i, j, k)] = sp.Symbol(f"c{i}{j}{k}")
    return coeffs


def create_interpolation_polynomial(coeffs):
    """Create the interpolation polynomial p(x, y, z)."""
    p_expr = sum(
        coeffs[(i, j, k)] * x_var**i * y_var**j * z_var**k
        for i in range(4)
        for j in range(4)
        for k in range(4)
    )
    return p_expr


def get_positions():
    """Get grid positions and function positions."""
    # Grid positions from -3 to +3 (inclusive)
    grid_positions = list(product(range(-3, 4), repeat=3))
    # Positions where we have function values and derivatives (from 0 to 1)
    fn_positions = list(product([0, 1], repeat=3))
    return grid_positions, fn_positions


def build_function_value_equations(p_expr, fn_positions):
    """Build equations for function values at specified positions."""
    equations = []
    for i, j, k in fn_positions:
        eq = sp.Eq(p_expr.subs({x_var: i, y_var: j, z_var: k}), DATA(i, j, k))
        equations.append(eq)
    return equations


def compute_derivatives(p_expr):
    """Compute symbolic derivatives of the interpolation polynomial."""
    dpdx = sp.diff(p_expr, x_var)
    dpdy = sp.diff(p_expr, y_var)
    dpdz = sp.diff(p_expr, z_var)
    dpdxdy = sp.diff(p_expr, x_var, y_var)
    dpdxdz = sp.diff(p_expr, x_var, z_var)
    dpdydz = sp.diff(p_expr, y_var, z_var)
    dpdxdydz = sp.diff(p_expr, x_var, y_var, z_var)
    return dpdx, dpdy, dpdz, dpdxdy, dpdxdz, dpdydz, dpdxdydz


def finite_difference_approximations():
    """Define finite difference approximations for derivatives based on the Maple code."""
    h = 1  # Grid spacing

    # Define the dx_3point and dx_5point functions
    def dx_3point(f, x0):
        return sp.Rational(1, 2) * (-f(x0 - 1) + f(x0 + 1))

    def dx_5point(f, x0):
        return sp.Rational(1, 12) * (
            f(x0 - 2) - 8 * f(x0 - 1) + 8 * f(x0 + 1) - f(x0 + 2)
        )

    # First derivatives using 5-point stencil
    def deriv_3d_dx(i, j, k):
        def f(mi):
            return DATA(i + mi, j, k)

        return dx_5point(f, 0)

    def deriv_3d_dy(i, j, k):
        def f(mj):
            return DATA(i, j + mj, k)

        return dx_5point(f, 0)

    def deriv_3d_dz(i, j, k):
        def f(mk):
            return DATA(i, j, k + mk)

        return dx_5point(f, 0)

    # Mixed derivatives using nested 5-point stencils
    def deriv_3d_dxy(i, j, k):
        def f(mi):
            def g(mj):
                return DATA(i + mi, j + mj, k)

            return dx_5point(g, 0)

        return dx_5point(f, 0)

    def deriv_3d_dxz(i, j, k):
        def f(mi):
            def g(mk):
                return DATA(i + mi, j, k + mk)

            return dx_5point(g, 0)

        return dx_5point(f, 0)

    def deriv_3d_dyz(i, j, k):
        def f(mj):
            def g(mk):
                return DATA(i, j + mj, k + mk)

            return dx_5point(g, 0)

        return dx_5point(f, 0)

    def deriv_3d_dxyz(i, j, k):
        def f(mi):
            def g(mj):
                def h(mk):
                    return DATA(i + mi, j + mj, k + mk)

                return dx_5point(h, 0)

            return dx_5point(g, 0)

        return dx_5point(f, 0)

    return (
        deriv_3d_dx,
        deriv_3d_dy,
        deriv_3d_dz,
        deriv_3d_dxy,
        deriv_3d_dxz,
        deriv_3d_dyz,
        deriv_3d_dxyz,
    )


def build_derivative_matching_equations(derivatives, fd_approximations, fn_positions):
    """Build equations for derivative matching at specified positions."""
    equations = []
    (
        dpdx,
        dpdy,
        dpdz,
        dpdxdy,
        dpdxdz,
        dpdydz,
        dpdxdydz,
    ) = derivatives
    (
        deriv_3d_dx,
        deriv_3d_dy,
        deriv_3d_dz,
        deriv_3d_dxy,
        deriv_3d_dxz,
        deriv_3d_dyz,
        deriv_3d_dxyz,
    ) = fd_approximations

    for i, j, k in fn_positions:
        eq_dx = sp.Eq(dpdx.subs({x_var: i, y_var: j, z_var: k}), deriv_3d_dx(i, j, k))
        eq_dy = sp.Eq(dpdy.subs({x_var: i, y_var: j, z_var: k}), deriv_3d_dy(i, j, k))
        eq_dz = sp.Eq(dpdz.subs({x_var: i, y_var: j, z_var: k}), deriv_3d_dz(i, j, k))
        eq_dxy = sp.Eq(
            dpdxdy.subs({x_var: i, y_var: j, z_var: k}), deriv_3d_dxy(i, j, k)
        )
        eq_dxz = sp.Eq(
            dpdxdz.subs({x_var: i, y_var: j, z_var: k}), deriv_3d_dxz(i, j, k)
        )
        eq_dyz = sp.Eq(
            dpdydz.subs({x_var: i, y_var: j, z_var: k}), deriv_3d_dyz(i, j, k)
        )
        eq_dxyz = sp.Eq(
            dpdxdydz.subs({x_var: i, y_var: j, z_var: k}), deriv_3d_dxyz(i, j, k)
        )
        equations.extend([eq_dx, eq_dy, eq_dz, eq_dxy, eq_dxz, eq_dyz, eq_dxyz])
    return equations


def solve_equations(equations, unknown_coeffs):
    """Solve the system of equations for the unknown coefficients."""
    a_matrix, b_vector = sp.linear_eq_to_matrix(equations, unknown_coeffs)
    solutions = sp.linsolve((a_matrix, b_vector), unknown_coeffs)
    # linsolve returns a FiniteSet
    solution = list(solutions)[0]
    coeff_dict = dict(zip(unknown_coeffs, solution))
    return coeff_dict


def collect_data_variables(grid_positions):
    """Collect data variables based on grid positions."""
    data_vars_local = [DATA(i, j, k) for i, j, k in grid_positions]
    return data_vars_local


def format_index(i):
    """Format index for variable names."""
    i = int(i)
    if i == 0:
        return "0"
    if i > 0:
        return f"p{i}"
    return f"m{abs(i)}"


class CustomPythonCodePrinter(PythonCodePrinter):
    """Custom code printer to ensure sp.Rational is used in the output."""

    def _print_Rational(self, expr):
        p_num, q_den = expr.p, expr.q
        return f"sp.Rational({p_num}, {q_den})"

    def _print_Symbol(self, expr):
        # Ensure symbols are printed as is, without any prefix
        return str(expr)

    def _get_loop_opening_ending(self, *args, **kwargs):
        """Stub method."""
        return "", ""

    def _rate_index_position(self, *args, **kwargs):
        """Stub method."""
        return 0


def output_expressions_to_py(coeff_exprs, cse_results, filename, coeffs_name_prefix):
    """
    Output expressions to a Python file.

    Parameters:
        coeff_exprs: List of tuples (data_variable, coefficient_expression)
        cse_results: Tuple (replacements, reduced_exprs)
        filename: Output filename
        coeffs_name_prefix: Prefix for coefficient variable names
    """
    # cse_results is a tuple (replacements, reduced_exprs)
    replacements, reduced_exprs = cse_results

    # Collect all free symbols to define at the top
    free_symbols = set()
    for expr in reduced_exprs:
        free_symbols.update(expr.free_symbols)
    for repl in replacements:
        free_symbols.update(repl[1].free_symbols)

    # Exclude symbols that will be defined in the script
    predefined_symbols = set(repl[0] for repl in replacements)
    predefined_symbols.update([dv for dv, _ in coeff_exprs])
    free_symbols -= predefined_symbols

    # Open file
    with open(filename, "w", encoding="utf-8") as file_obj:
        # Write the imports and symbol definitions
        file_obj.write("import sympy as sp\n\n")
        # Define free symbols
        if free_symbols:
            file_obj.write("# Define free symbols\n")
            for sym in sorted(free_symbols, key=str):
                file_obj.write(f"{sym} = sp.Symbol('{sym}')\n")
            file_obj.write("\n")
        # Create an instance of the custom code printer
        printer = CustomPythonCodePrinter()

        # Write the temporary variables
        for _, (var, expr) in enumerate(replacements):
            expr_code = printer.doprint(expr)
            line = f"{var} = {expr_code}\n"
            file_obj.write(line)
        file_obj.write("\n")
        # Write the coefficients
        for _, ((dv_sym, _), expr) in enumerate(zip(coeff_exprs, reduced_exprs)):
            # Get the position from the data variable name
            dv_name = str(dv_sym)
            # dv_name is of the form 'DATA_i_j_k'
            # We need to extract i, j, k
            _, i_str, j_str, k_str = dv_name.split("_")
            i_idx = int(i_str)
            j_idx = int(j_str)
            k_idx = int(k_str)
            # Now, construct the coefficient name
            i_formatted = format_index(i_idx)
            j_formatted = format_index(j_idx)
            k_formatted = format_index(k_idx)
            coeff_name = (
                f"{coeffs_name_prefix}coeff_{i_formatted}_{j_formatted}_{k_formatted}"
            )
            # Now, expr may be an expression involving temporary variables
            expr_code = printer.doprint(expr)
            line = f"{coeff_name} = {expr_code}\n\n"
            file_obj.write(line)


def process_derivative(
    deriv_name,
    deriv_expr,
    coeff_dict_local,
    data_vars_local,
    output_dir_local,
):
    """
    Process a derivative and output expressions.

    Parameters:
        deriv_name: Name of the derivative
        deriv_expr: Expression of the derivative
        coeff_dict_local: Dictionary mapping coefficients to their solutions
        data_vars_local: List of data variables
        output_dir_local: Directory to output the files
    """
    print(f"Started processing derivative '{deriv_name}'")
    start_time_local = time.time()

    deriv_subs = deriv_expr.subs(coeff_dict_local)
    deriv_collected = sp.collect(sp.expand(deriv_subs), data_vars_local)
    coeff_exprs_deriv = []
    for dv_sym in data_vars_local:
        coeff_expr = deriv_collected.coeff(dv_sym)
        if coeff_expr != 0:
            coeff_exprs_deriv.append((dv_sym, coeff_expr))
    expressions_deriv = [expr for _, expr in coeff_exprs_deriv]
    cse_results_deriv = sp.cse(expressions_deriv)
    coeffs_filename = os.path.join(output_dir_local, f"coeffs_{deriv_name}.py")
    coeffs_name_prefix_local = f"coeffs_{deriv_name}__"
    os.makedirs(os.path.dirname(coeffs_filename), exist_ok=True)
    output_expressions_to_py(
        coeff_exprs_deriv, cse_results_deriv, coeffs_filename, coeffs_name_prefix_local
    )

    end_time_local = time.time()
    print(
        f"Finished processing derivative '{deriv_name}' in {end_time_local - start_time_local:.2f} seconds"
    )


# Define symbols globally to avoid redefining from outer scope
x_var, y_var, z_var = sp.symbols("x y z")


def main():
    """Main function to orchestrate the computation of symbolic coefficients."""
    # Create coefficients
    coeffs = create_coefficients()

    # Create interpolation polynomial
    p_expr = create_interpolation_polynomial(coeffs)

    # Get positions
    grid_positions, fn_positions = get_positions()

    # Build equations for function values
    equations = build_function_value_equations(p_expr, fn_positions)

    # Compute symbolic derivatives
    derivatives = compute_derivatives(p_expr)

    # Define finite difference approximations
    fd_approximations = finite_difference_approximations()

    # Build equations for derivative matching
    derivative_equations = build_derivative_matching_equations(
        derivatives, fd_approximations, fn_positions
    )
    equations.extend(derivative_equations)

    # List of unknown coefficients
    unknown_coeffs = list(coeffs.values())

    # Solve the system of equations
    coeff_dict_local = solve_equations(equations, unknown_coeffs)

    # Collect data variables
    data_vars_local = collect_data_variables(grid_positions)

    return p_expr, coeff_dict_local, data_vars_local, derivatives


if __name__ == "__main__":
    # Run the main function to get necessary variables
    p_expr, coeff_dict, data_vars, derivatives = main()

    # Process the interpolation polynomial (I)
    output_dir = "3d_coeffs/3d_cube_order3"
    print("Constructing interpolating polynomial coefficients...")
    start_time = time.time()

    os.makedirs(output_dir, exist_ok=True)
    p_subs = p_expr.subs(coeff_dict)
    p_collected = sp.collect(sp.expand(p_subs), data_vars)

    coeff_exprs_I = []
    for dv_sym in data_vars:
        coeff_expr = p_collected.coeff(dv_sym)
        if coeff_expr != 0:
            coeff_exprs_I.append((dv_sym, coeff_expr))
    expressions_I = [expr for _, expr in coeff_exprs_I]
    cse_results_I = sp.cse(expressions_I)
    coeffs_I_filename = os.path.join(output_dir, "coeffs_I.py")
    coeffs_name_prefix_I = "coeffs_I__"
    print("Outputting coeffs_I...")
    output_expressions_to_py(
        coeff_exprs_I, cse_results_I, coeffs_I_filename, coeffs_name_prefix_I
    )

    end_time = time.time()
    print(
        f"Finished constructing interpolating polynomial coefficients in {end_time - start_time:.2f} seconds.\n"
    )

    # List of derivatives to process (only first derivatives and mixed first derivatives)
    derivative_list = [
        ("dx", derivatives[0]),
        ("dy", derivatives[1]),
        ("dz", derivatives[2]),
        ("dxx", sp.diff(p_expr, x_var, x_var)),
        ("dxy", derivatives[3]),
        ("dxz", derivatives[4]),
        ("dyy", sp.diff(p_expr, y_var, y_var)),
        ("dyz", derivatives[5]),
        ("dzz", sp.diff(p_expr, z_var, z_var)),
    ]

    # Prepare arguments for multiprocessing
    args_list = [
        (deriv_name, deriv_expr, coeff_dict, data_vars, output_dir)
        for deriv_name, deriv_expr in derivative_list
    ]

    print(f"Processing {len(derivative_list)} derivatives in parallel...")
    start_time = time.time()
    with multiprocessing.Pool() as pool:
        pool.starmap(process_derivative, args_list)
    end_time = time.time()
    print(
        f"Finished processing all derivatives in {end_time - start_time:.2f} seconds."
    )
