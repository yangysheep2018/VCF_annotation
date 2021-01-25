import argparse
import vcf
import json
import urllib.request
import jsonpath
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", required=True, type=str, help="input .vcf file")


args = vars(parser.parse_args())
vcf_file = args["i"]


base_url = "http://exac.hms.harvard.edu/rest/variant/"


##http://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html to get the sorted 
##the listed from most harmful to least
ll = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost", 
      "start_lost","transcript_amplification","inframe_insertion","inframe_deletion","missense_variant", 
      "protein_altering_variant","splice_region_variant","incomplete_terminal_codon_variant","start_retained_variant",
      "stop_retained_variant","synonymous_variant","initiator_codon_variant","coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant",
      "3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","NMD_transcript_variant",
      "non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification",
      "TF_binding_site_variant","regulatory_region_ablation","regulatory_region_amplification","feature_elongation",
      "regulatory_region_variant","feature_truncation","intergenic_variant"]

##create a dictionary with key as effect and value as number from 0 to 36,
##the value is smaller, the consequence is more harmful
dict_effect_importance = { ll[i] : i for i in range(0, len(ll))}


df = pd.DataFrame(columns=['Chrom', 'Pos', 'Ref', 'Alt', 'Variant_Type', 'Variant_Effect', 'DP', 'AO', 'Ratio', 'AF_ExAC'])

j = 0


##load input vcf file
vcf_reader = vcf.Reader(filename=vcf_file)


for record in vcf_reader:

    snp = (str(record.CHROM) + "-" + str(record.POS) + "-" + str(record.REF[0]) + "-" + str(record.ALT[0]))
#     print(snp)

    contents = urllib.request.urlopen(base_url + snp).read()
    result_object = json.loads(contents)

    effect_pool = jsonpath.jsonpath(result_object,'$..major_consequence')

    ##by referring the dict_effect_importance rank, get the effect pool containing all the potential
    ##consequences the most deleterious possibility
    Importance = 37 # biggest number of effect
    Effect = 'NA'
    if effect_pool != False:
        for i in effect_pool:
            if dict_effect_importance[i] < Importance:
                Importance = dict_effect_importance[i]
                Effect = i

    ##initialize the ratio between variants and reference reads depth to -1,
    ##if reference read depth equals to 0, then the ratio equals to -1
    
    ratio = -1
    if record.INFO["RO"]>0:
        ratio = float((str(record.INFO["AO"]).strip("[]"))[0])/float(record.INFO["RO"])

    df.loc[j]=[record.CHROM,record.POS,record.REF[0],record.ALT[0],record.INFO["TYPE"][0],
               Effect,record.INFO['DP'],str(record.INFO['AO']).strip("[]"),
               ratio,str(jsonpath.jsonpath(result_object,'$..allele_freq')).strip("[]")]
    j = j+1
    # print(j)

df.to_csv("vcf_output.csv", index = False)










