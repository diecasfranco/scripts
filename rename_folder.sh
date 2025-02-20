#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <folder_path> -l <list_file>"
    exit 1
}

# Parse command line arguments
while getopts "i:l:" opt; do
    case ${opt} in
        i )
            DIR_PATH=$OPTARG
            ;;
        l )
            LIST_FILE=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Check if both arguments were provided
if [ -z "$DIR_PATH" ] || [ -z "$LIST_FILE" ]; then
    usage
fi

# Read the rename list and process each line
while IFS= read -r line; do
    # Split the line into old name and new name
    old_name=$(echo $line | awk '{print $1}')
    new_name=$(echo $line | awk '{print $2}')
    
    # Rename the folder
    mv "${DIR_PATH}/${old_name}" "${DIR_PATH}/${new_name}"

    # Rename the FASTA file inside the folder and change the header
    for fasta_file in "${DIR_PATH}/${new_name}"/*.fna; do
        if [ -f "$fasta_file" ]; then
            # Rename the FASTA file to match the folder name
            mv "$fasta_file" "${DIR_PATH}/${new_name}/${new_name}.fna"
            
            # Change the FASTA header
            sed -i "1s/.*/>${new_name}/" "${DIR_PATH}/${new_name}/${new_name}.fna"
        fi
    done
done < "$LIST_FILE"
