#!/bin/bash
#move files to data repostiory
mv src/workflow/* data/ && mv build/MeshBuilder data/ && mv src/metis_builds/METIS/build/programs/gpmetis data/
cd data

# Prompt the user to input information
echo "Enter the path to the .in file to be used by MeshBuilder: "
read file_path
echo "Enter the number of computer nodes for partitioning: "
read nn
echo "Enter partitioning method,
		1 --> SurfaceFlow (SF),
		2 --> Surface-Subsurface Flow (SSF),
        3 --> Surface-Subsurface Flow with Headwaters (SSF-H): "
read OPT_Part
echo "Enter simulation basename: "
read basename

# Check if the file exists
if [ -f "$file_path" ]; then
    # If the file exists, execute the programs with the provided inputs
    echo "Executing MeshBuilder with the file: $file_path"
    ./MeshBuilder "$file_path"  # Replace "MeshBuilder" with the actual program name

    echo "Executing run_metis.zsh with parameters: $nn, $OPT_Part, $basename"
    ./run_metis.zsh "$nn" "$OPT_Part" "$basename"
else
    # If the file does not exist, display an error message
    echo "Error: The specified file does not exist."
fi
