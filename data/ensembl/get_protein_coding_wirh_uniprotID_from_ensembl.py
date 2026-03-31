import requests
import pandas as pd
from io import StringIO

xml_query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
				
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "upstream_flank" value = "500"/>
		<Filter name = "with_uniprotswissprot" excluded = "0"/>
		<Attribute name = "coding" />
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "uniprotswissprot" />
		<Attribute name = "strand" />
		<Attribute name = "transcript_start" />
		<Attribute name = "transcript_end" />
		<Attribute name = "ensembl_transcript_id" />
	</Dataset>
</Query>cd 
"""

url = "https://www.ensembl.org/biomart/martservice"
response = requests.get(url, params={"query": xml_query}, timeout=120)
response.raise_for_status()

text = response.text

# Αν υπάρχει completion stamp, αφαίρεσέ το πριν το διάβασμα ως TSV
text = text.replace("[success]", "").strip()

df = pd.read_csv(StringIO(text), sep="\t")
print(df.head())
print(df.shape)

# αποθήκευση
df.to_csv("ensembl_biomart_human_swissprot_coding_upstream500.tsv", sep="\t", index=False)