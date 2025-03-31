def remove_unicode(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as infile, open(output_file, 'w', encoding='utf-8') as outfile:
        for line in infile:
            # Remove non-ASCII characters
            cleaned_line = ''.join(char for char in line if ord(char) < 128)
            outfile.write(cleaned_line)

# Specify file paths
input_path = 'RP3_SK.inp'  # Replace with your file
output_path = 'RP3_SK.inp'

# Remove Unicode characters
remove_unicode(input_path, output_path)
print(f"Cleaned file saved to {output_path}")

