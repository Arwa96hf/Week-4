 #!/bin/bash

# Loop over each genome file in the directory
for file in *.fna; do
    echo "Processing $file"
    python3 ORF.py "$file" > "${file%.fna}_output.txt"
done

