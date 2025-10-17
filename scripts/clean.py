import os
import glob
import re


def delete_generated_files(base_dir):
    # Find all .pyx files in the directory structure
    pyx_files = glob.glob(os.path.join(base_dir, "**", "*.pyx"), recursive=True)

    # Track files to delete
    files_to_delete = []

    # Regular expression for shared object files (e.g., .so)
    so_regex = re.compile(r".*\.cpython-\d{2}.*\.so")

    for pyx_file in pyx_files:
        # Remove the extension to get the base file name
        base_name = os.path.splitext(pyx_file)[0]

        # Possible generated files
        generated_files = [
            f"{base_name}.c",
            f"{base_name}.cpp",
        ]

        # Add the matching .so files using the base name
        so_files = glob.glob(f"{base_name}*.so")
        for so_file in so_files:
            if so_regex.match(so_file):
                generated_files.append(so_file)

        # Add existing generated files to the deletion list
        for gen_file in generated_files:
            if os.path.exists(gen_file):
                files_to_delete.append(gen_file)

    # Delete files and report
    for file in files_to_delete:
        try:
            os.remove(file)
            print(f"Deleted: {file}")
        except Exception as e:
            print(f"Error deleting {file}: {e}")


if __name__ == "__main__":
    # Set the base directory of your project
    project_dir = "."

    # Call the function
    delete_generated_files(project_dir)
