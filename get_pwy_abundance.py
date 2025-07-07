#!/usr/bin/env python

import pandas as pd
import math
from subprocess import PIPE, run
import argparse

parser = argparse.ArgumentParser(
    description="Run fungal functional profiling, Input is the diamond blastx output.\n"
)

parser.add_argument(
    "-i",
    "--input",
    dest="input_path",
    action="store",
    default=None,
    help="Input blastx output path.\n",
    required=True,
)
parser.add_argument(
    "-s",
    "--script",
    dest="script_path",
    action="store",
    default=None,
    help="folder storing needed tables and scripts.\n",
    required=True,
)


parser.add_argument(
    "-o1",
    "--output1",
    dest="output_path1",
    action="store",
    default="",
    help="Output path of the fungal pathway profiling.\n",
    required=True,
)
parser.add_argument(
    "-o2",
    "--output2",
    dest="output_path2",
    action="store",
    default="",
    help="Output path of the fungal pathway_class profiling.\n",
    required=True,
)
parser.add_argument(
    "-o3",
    "--output3",
    dest="output_path3",
    action="store",
    default="",
    help="Output path of the fungal pathway_type profiling.\n",
    required=True,
)
parser.add_argument(
    "-o4",
    "--output4",
    dest="output_path4",
    action="store",
    default="",
    help="Output path of the fungal full annotation.\n",
    required=True,
)

option = parser.parse_args()

sample_to_process = option.input_path
scriptdir = option.script_path
output1 = option.output_path1
output2 = option.output_path2
output3 = option.output_path3
output4 = option.output_path4

print("LOADING DATABASES")
ann_path = scriptdir + "/jgi_ann_05-2022_reordered.tab"
uniprot_catalog = scriptdir + "/ncbi_uniprot_species.txt"
taxa_path = scriptdir + "/taxonomy_for_function.csv"

column_names = [
    "proteinId",
    "ecNum",
    "definition",
    "uniprotID",
    "KO",
    "pathway",
    "pathway_class",
    "pathway_type",
    "speciesID",
    "catalyticActivity",
    "cofactors",
    "associatedDiseases",
]
jgi_ann = pd.read_csv(ann_path, delimiter="\t", names=column_names, index_col=False)
blastx_out = pd.read_csv(
    sample_to_process,
    index_col=False,
    delimiter="\t",
    names=[
        "qid",
        "rid",
        "id",
        "lenght",
        "mismatch",
        "gapopne",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bscore",
    ],
)
rid_to_uniprot_species = pd.read_csv(
    uniprot_catalog,
    delimiter="\t",
    names=["rid", "uniprotId", "speciesID"],
    dtype={"first_column": "str"},
    index_col=False,
)
taxa = pd.read_csv(
    taxa_path, delimiter="\t", names=["speciesID", "taxa"], index_col=False
)


print("PARSING BLASTX RESULT")
if len(blastx_out) >= 1:
    # raw_counts=blastx_out['rid'].value_counts()
    raw_counts = (
        blastx_out["rid"].value_counts().rename_axis("rid").reset_index(name="counts")
    )
    # raw_counts_clean=pd.merge(raw_counts,id_catalog)
    column_names = [
        "proteinId",
        "ecNum",
        "definition",
        "uniprotID",
        "KO",
        "pathway",
        "pathway_class",
        "pathway_type",
        "speciesID",
        "catalyticActivity",
        "cofactors",
        "associatedDiseases",
    ]
    pwy_counts = pd.DataFrame(columns=column_names, dtype=object)

    for i in range(len(raw_counts)):
        # print(f"CHECKING {raw_counts.iloc[i]['rid']}")
        inter = jgi_ann.loc[jgi_ann["proteinId"] == raw_counts.iloc[i]["rid"]]
        # print(f'INTER JGI: {inter}')
        ## if len(inter)==0, if starts with JGI, assign NAs to all fields of inter; else if starts with NCBI, map to uniprot catalog, then use the R script to get KEGG path; else if starts with neither, use the R script to get KEGG path
        if len(inter) == 0:
            speciesID = rid_to_uniprot_species.loc[
                rid_to_uniprot_species["rid"] == raw_counts.iloc[i]["rid"]
            ]["speciesID"]
            # print(f'Species ID: {speciesID}')
            if len(speciesID) == 0:
                speciesID = "unclassified"
            else:
                speciesID = list(speciesID)[0]  # sometimes it is duplicated
            if raw_counts.iloc[i]["rid"].startswith("jgi"):
                # if starts with JGI and was not found in jgi_ann file: set to unidentified
                inter = [
                    raw_counts.iloc[i]["rid"],
                    "NA",
                    "NA",
                    "NA",
                    "unidentified",
                    "unidentified",
                    "unidentified",
                    "unidentified",
                    speciesID,
                    "NA",
                    "NA",
                    "NA",
                ]
                # print(f'INTER JGI UNIDENT {inter}')
                inter = pd.DataFrame([inter], columns=column_names)
            else:
                # else try mapping with uniprot
                uptIDs = rid_to_uniprot_species.loc[
                    rid_to_uniprot_species["rid"] == raw_counts.iloc[i]["rid"]
                ]["uniprotId"]
                # print(f'UNIPROT ID: {uptID}')
                if len(uptIDs) == 0:
                    uptIDs = "unidentified"
                    inter = [
                        raw_counts.iloc[i]["rid"],
                        "NA",
                        "NA",
                        "NA",
                        "unidentified",
                        "unidentified",
                        "unidentified",
                        "unidentified",
                        speciesID,
                        "NA",
                        "NA",
                        "NA",
                    ]
                    # print(f'INTER UNIPROT EMPTY {inter}')
                    inter = pd.DataFrame([inter], columns=column_names)
                else:
                    inter = []
                    for uptID in uptIDs:
                        # merge raw_counts.iloc[i]['rid']] with NCBI_UNIPROT_CATALOG (link to uniprot accession)
                        # uptID=uptID.to_string(index=False).strip()
                        uptID = str(uptID).strip()
                        if uptID == "nan":  # UNIPROT accession is NA
                            row = [
                                raw_counts.iloc[i]["rid"],
                                "NA",
                                "NA",
                                "NA",
                                "unidentified",
                                "unidentified",
                                "unidentified",
                                "unidentified",
                                speciesID,
                                "NA",
                                "NA",
                                "NA",
                            ]
                            # print(f'INTER UNIPROT NAN {inter}')
                            inter.append(row)
                        else:  # invoke the keggConv script to get pwy info, return pwy, pwy class
                            # print(f'Runing keggConv for {uptID}')
                            command = ["keggConv.R", uptID]
                            result = run(
                                command,
                                stdout=PIPE,
                                stderr=PIPE,
                                universal_newlines=True,
                            )
                            # print(f'Keggconv result: {result}')
                            if result.stderr:
                                print(result.stderr)
                            conv = result.stdout.split("|")
                            # print(f'Keggconv split result: {conv}')
                            for c in range(0, len(conv), 3):
                                # print(f'Keggconv result rowset {c}')
                                try:
                                    rows = [
                                        raw_counts.iloc[i]["rid"],
                                        "NA",
                                        "NA",
                                        uptID,
                                        conv[c].strip("\n"),
                                        conv[c + 1].strip("\n"),
                                        conv[c + 2].strip("\n"),
                                        "unidentified",
                                        speciesID,
                                        "NA",
                                        "NA",
                                        "NA",
                                    ]
                                except IndexError as e:
                                    print(f"ID {raw_counts.iloc[i]['rid']}")
                                    print(f"UniprotID {uptID}")
                                    print(f"Keggconv result: {result}")
                                    print(f"Keggconv split result: {conv}")
                                    raise (e)
                                # print(f'{rows}')
                                inter.append(rows)

                    inter = pd.DataFrame(inter, columns=column_names)

        inter = inter.assign(counts=raw_counts.iloc[i]["counts"] / len(inter))
        # print(f"FINAL INTER FOR {raw_counts.iloc[i]['rid']}: {inter}")
        pwy_counts = pwy_counts.append(inter, ignore_index=True, sort=False)

        # input('Press enter to continue ...') #debugging purposes
    # raise()

    def get_stratified_taxa(ann_type):
        annotation = (
            pwy_counts.reset_index()
            .groupby([ann_type, "speciesID"])
            .agg({"counts": "sum"})
            .reset_index()
        )
        annotation.columns = [ann_type, "speciesID", "counts"]
        if len(annotation) >= 1:
            column_names = [ann_type, "counts"]
            ann_abd = pd.DataFrame(columns=column_names, dtype=object)
            for p in annotation[ann_type].unique():
                ann_subset = annotation[annotation[ann_type] == p]
                ann_taxa = pd.merge(taxa, ann_subset)
                # taxa_list=pd.DataFrame(ann_taxa.taxa.str.split('|').tolist(), dtype=object)
                # column_names=['genus_species']
                # taxa_table=pd.DataFrame(taxa_list[6].astype(str),columns=column_names, dtype=object)
                # stratified_pwy=pd.concat([ann_taxa,taxa_table],axis=1)
                ann_taxa[ann_type] = (
                    ann_taxa[ann_type].astype(str) + "|" + ann_taxa["taxa"].astype(str)
                )
                ann_taxa = ann_taxa[[ann_type, "counts"]]
                all_counts = ann_taxa["counts"].sum()
                ann_taxa.loc[-1] = [p, all_counts]  # adding a row
                ann_taxa.index = ann_taxa.index + 1  # shifting index
                ann_taxa.sort_index(inplace=True)
                stratified_pwy = (
                    ann_taxa.groupby([ann_type]).agg({"counts": "sum"}).reset_index()
                )
                ann_abd = ann_abd.append(stratified_pwy, ignore_index=True)
            return ann_abd
        else:
            print("no pathway found\n")
            pass

    pwy_abd = get_stratified_taxa("pathway")
    pwyCls_abd = get_stratified_taxa("pathway_class")
    pwyTyp_abd = get_stratified_taxa("pathway_type")

    pwy_abd.to_csv(output1, index=False, sep="\t")
    pwyCls_abd.to_csv(output2, index=False, sep="\t")
    pwyTyp_abd.to_csv(output3, index=False, sep="\t")
    pwy_counts.to_csv(output4, index=False, sep="\t")
else:
    print("no protein found\n")
    pass
