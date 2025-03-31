import PyPDF2

# File paths
pdf_path = 'ao3c05087_si_001.pdf'
txt_path = 'mechanism.inp'

# Open and extract text
with open(pdf_path, 'rb') as pdf_file, open(txt_path, 'w') as txt_file:
    reader = PyPDF2.PdfReader(pdf_file)
    for page in reader.pages:
        text = page.extract_text()
        if text:  # Write only if there's text
            txt_file.write(text + '\n')

print(f'Text extracted and saved to {txt_path}')

