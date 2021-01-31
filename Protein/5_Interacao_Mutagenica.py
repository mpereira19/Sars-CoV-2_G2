import vcf
from tabulate import tabulate
from gtfparse import read_gtf

print('Variante Portuguesa do SARS-CoV-2:')
vcf_reader = vcf.Reader(open("ERR4157959.sars_cov.raw.vcf", "r"))
for record in vcf_reader:
    print(record)

print('--------------------------------------------------------------------')

print('Variante Espanhola do SARS-CoV-2:')
vcf_reader = vcf.Reader(open("ERR4395294.sars_cov.raw.vcf", "r"))
for record in vcf_reader:
    print(record)


print()
df = read_gtf("Sars_cov_2.ASM985889v3.101.gtf")
#print(tabulate(df, headers='keys', tablefmt='github'))

df_genes = df[df["gene_name"] == "ORF8"]
print(tabulate(df_genes, headers='keys', tablefmt='github'))
df_genes1 = df[df["gene_name"] == "N"]
print(tabulate(df_genes1, headers='keys', tablefmt='github'))
