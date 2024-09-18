import importlib
import sys

import nrpy.validate_expressions.validate_expressions as ve


def main():
    # List of base filenames without extensions
    bases = ["dx", "dxx", "dxy", "dxz", "dy", "dyy", "dyz", "dz", "dzz", "I"]

    for base in bases:
        coeffs_module_name = f"3d_coeffs.3d_cube_order2.coeffs_{base}"
        trusted_coeffs_module_name = f"AEILocalInterp_output.trusted_coeffs_{base}"

        try:
            # Dynamically import the coeffs_* module
            coeffs = importlib.import_module(coeffs_module_name)
        except ImportError:
            print(f"Error: Module '{coeffs_module_name}.py' not found.")
            sys.exit(1)

        try:
            # Dynamically import the trusted_coeffs_* module
            trusted_coeffs = importlib.import_module(trusted_coeffs_module_name)
        except ImportError:
            print(f"Error: Module '{trusted_coeffs_module_name}.py' not found.")
            sys.exit(1)

        # Retrieve all attributes from coeffs module that start with 'coeffs_[base]__'
        coeffs_vars = [
            var
            for var in dir(coeffs)
            if var.startswith(f"coeffs_{base}__") and not var.startswith("__")
        ]

        if not coeffs_vars:
            print(
                f"No variables found in '{coeffs_module_name}.py' matching pattern 'coeffs_{base}__*'."
            )
            sys.exit(1)

        for var in coeffs_vars:
            # Construct the corresponding trusted_coeffs_* variable name
            suffix = var[len(f"coeffs_{base}__") :]
            trusted_var = f"trusted_coeffs_{base}__{suffix}"

            # Check if the trusted_var exists in the trusted_coeffs module
            if not hasattr(trusted_coeffs, trusted_var):
                print(
                    f"Warning: '{trusted_var}' not found in '{trusted_coeffs_module_name}.py'. Skipping validation."
                )
                sys.exit(1)

            # Retrieve the variables from both modules
            coeffs_value = getattr(coeffs, var)
            trusted_coeffs_value = getattr(trusted_coeffs, trusted_var)

            # Print validation message
            print(f"Validating '{trusted_var}' against '{var}'")

            try:
                # Perform the equality assertion
                ve.assert_equal(trusted_coeffs_value, coeffs_value)
                print(f"Validation Passed: '{trusted_var}' == '{var}'\n")
            except AssertionError as e:
                print(f"Validation Failed: '{trusted_var}' != '{var}'\n")
                print(e)
                sys.exit(1)

    print("All validations completed successfully.")


if __name__ == "__main__":
    main()
