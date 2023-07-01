#Ejercicio 2.1#
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

def source(id):
    # Conexión al NCBI
    Entrez.email = "carloz1palacios@gmail.com"  # Reemplaza con tu correo electrónico registrado en NCBI
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    
    # Extracción del nombre del organismo fuente
    source_name = record.annotations["source"]
    
    # Cálculo de la frecuencia de cada especie
    species_frequency = {}
    for feature in record.features:
        if feature.type == "source":
            organism = feature.qualifiers.get("organism", ["Unknown Organism"])[0]
            species_frequency[organism] = species_frequency.get(organism, 0) + 1
            print(species_frequency)
    
    # Guardar el resultado en un archivo CSV
    df = pd.DataFrame.from_dict(species_frequency, orient="index", columns=["Frequency"])
    df.index.name = "Species"
    df.to_csv("results/source.csv")
    
    # Imprimir el nombre del organismo fuente
    print("Source organism:", source_name)


##Ejercicio 2.2##

import pandas as pd
from Bio import Entrez, SeqIO, SeqUtils

def sequences(id):
    # Conexión al NCBI
    Entrez.email = "emily.cherrez@est.ikiam.edu.ec" 
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    
    # Extracción del nombre del organismo fuente
    sequence_name = record.description
    
    # Extracción de la secuencia de ADN
    dna_sequence = record.seq
    
    # Verificar que la longitud de la secuencia sea divisible por tres
    if len(dna_sequence) % 3 != 0:
        print("La longitud de la secuencia no es divisible por tres. No se puede realizar la traducción completa.")
        return
    
    # Traducción de la secuencia de ADN a péptidos
    peptide_sequence = dna_sequence.translate()
    
    # Separación de péptidos
    peptides = str(peptide_sequence).split('*')
    peptides = [p for p in peptides if len(p) > 0]
    
    # Cálculo del peso molecular y frecuencia de cada péptido
    peptide_data = []
    for peptide in peptides:
        mw = SeqUtils.molecular_weight(peptide)
        peptide_data.append({'Peptide': peptide, 'Molecular Weight': mw})
    
    # Guardar el resultado en un archivo CSV
    df = pd.DataFrame(peptide_data)
    df.to_csv("results/peptides.csv", index=False)
    
    # Imprimir el nombre de la secuencia
    print("Sequence name:", sequence_name)

# Número de acceso de ejemplo
accession_number = "NC_064536.1"  # Reemplaza con tu número de acceso real

# Ejecutar la función para el número de acceso proporcionado
sequences(accession_number)