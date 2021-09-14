import os
import yaml
import argparse
import pandas as pd



def get_settings(config_file):
    """
    Load in settings from config file.
    """
    stream = open(config_file, 'r')
    settings = yaml.load(stream, Loader=yaml.Loader)
    return settings


def get_args():
    """
    Get command line arguments. 
    Make a flag for every gene name within the config file.
    Validate and process input arguments.
    """
    parser = argparse.ArgumentParser(
        description='Makes an expanded list of variants for checking in suspected mosiac cases',
        epilog='See https://github.com/erikwaskiewicz/filter_potential_mosiacs/blob/master/README.md for more details'
    )
    parser.add_argument('run_folder', action='store', help='Path to run folder in the format /data/results/<run_id>/<panel>/')
    parser.add_argument('sample_id', action='store', help='Sample ID e.g. 19M123456')

    # make flag for each gene in config
    for gene in GENES:
        parser.add_argument(('--'+gene), action='store_true', default=False, help=f'Include {gene} variants in the output')

    args = parser.parse_args()

    # check run path and sample ids are valid and parse run id from path
    run_folder = args.run_folder
    sample_id = args.sample_id
    
    if os.path.isdir(run_folder):
        run_path = os.path.abspath(run_folder)
        run_id = run_path.split('/')[-2]
    else:
        raise ValueError('The folder you entered isnt a folder')

    if not os.path.isdir(os.path.join(run_folder, sample_id)):
        raise ValueError(f'There isnt a sample {sample_id} on this run')
    
    # get a list of all genes selected
    genes = []
    for gene in GENES:
        if vars(args)[gene]:
            genes.append(gene)
    if len(genes) == 0:
        raise ValueError('No genes selected')

    # package all arguments into dictionary
    args_dict = {'run_path': run_path, 'run_id': run_id, 'sample_id': sample_id, 'genes': genes}

    return args_dict

    
def load_vcf(bcftools_command, genes=None):
    """
    Take a bcftools command and sort the output into a dataframe. 
    If a dictionary of regions is provided then apply the panel.
    """
    # run subprocess and capture stdout
    data = os.popen(bcftools_command).read()

    # format bcftools output into a list of rows, split rows into lists of columns
    output_list = []
    for line in data.split('\n'):
        line_split = line.split(';')
        if len(line_split) > 1:
            # convert pos into an integer so that it isnt converted incorrectly into the dataframe
            line_split[1] = int(line_split[1])
            output_list.append(line_split)

    # add all output to a pandas dataframe
    df = pd.DataFrame(data=output_list, columns=['chr', 'pos', 'ref', 'alt', 'AD', 'DP'])

    # apply a bed file to the data- chr, start and end come from the list within the config file
    if genes:
        df_panel = pd.DataFrame(columns=['chr', 'pos', 'ref', 'alt', 'AD', 'DP'])
        for g in genes:
            chrom = GENES[g][0]
            start = GENES[g][1]
            end = GENES[g][2]
            temp_df = df.loc[lambda df: df.chr == chrom, :].loc[lambda df: df.pos > start, :].loc[lambda df: df.pos < end, :]
            df_panel = df_panel.append(temp_df, ignore_index=True)
        return df_panel
        
    return df


def main():
    # arguments
    args = get_args()
    run_path = args['run_path']
    run_id = args['run_id']
    sample_id = args['sample_id']
    genes = args['genes']
    
    # testing variables
    #run_id = '190227_M00766_0196_000000000-CCD8K'
    #sample_id = '19M01966'
    #genes = ['NF1', 'NF2']
    

    ### First VCF

    vcf_1 = f'{run_path}/post_processing/results/annotated_vcf/{run_id}_anno.vcf.gz'
    bcftools_command_1 =  f'{BCFTOOLS_PATH} view -Ou -s {sample_id} {vcf_1} | {BCFTOOLS_PATH} query -f "%CHROM;%POS;%REF;%ALT;[ %AD];[ %DP]\n"'
    df_vcf_1 = load_vcf(bcftools_command_1, genes=genes)


    ### Second VCF

    vcf_2 = f'{run_path}/post_processing/results/gvcf/{run_id}_merged.vcf.gz'
    bcftools_command_2 =  f'{BCFTOOLS_PATH} view -Ou -s {sample_id} {vcf_2} | {BCFTOOLS_PATH} query -f "%CHROM;%POS;%REF;%ALT;[ %AD];[ %DP]\n"'
    df_vcf_2 = load_vcf(bcftools_command_2, genes=genes)

    # make list of all variants that have been seen in vcf 1
    df_vcf_1['hgvs'] = df_vcf_1[['chr', 'pos', 'ref']].apply(lambda x: '-'.join(x.map(str)), axis=1)
    df_vcf_2['hgvs'] = df_vcf_2[['chr', 'pos', 'ref']].apply(lambda x: '-'.join(x.map(str)), axis=1)
    seen_list = df_vcf_1['hgvs']

    # filter out any variants that were seen in vcf 1
    df_vcf_2_filtered = df_vcf_2[~df_vcf_2['hgvs'].isin(seen_list)]


    ### Third VCF

    vcf_3 = f'{run_path}/post_processing/results/gvcf/{run_id}_{sample_id}.g.vcf.gz'
    bcftools_command_3 =  f'{BCFTOOLS_PATH} query -f "%CHROM;%POS;%REF;%ALT;[ %AD];[ %DP]\n" {vcf_3}'
    df_vcf_3 = load_vcf(bcftools_command_3, genes=genes)

    # remove all ref alleles and remove trailing <NON_REF> text
    df_vcf_3 = df_vcf_3.loc[lambda df: df.alt != '<NON_REF>', :]
    df_vcf_3['alt'] = df_vcf_3['alt'].apply(lambda x: x.rstrip(',<NON_REF>'))

    # update seen list with variants from vcf 2
    seen_list = seen_list.append(df_vcf_2_filtered['hgvs'])

    # filter out any variants that were in either vcf 1 or 2
    df_vcf_3['hgvs'] = df_vcf_3[['chr', 'pos', 'ref']].apply(lambda x: '-'.join(x.map(str)), axis=1)
    df_vcf_3_filtered = df_vcf_3[~df_vcf_3['hgvs'].isin(seen_list)]


    ### Write to file

    # Remove hgvs columns
    df_vcf_1 = df_vcf_1.drop('hgvs', axis=1)
    df_vcf_2_filtered = df_vcf_2_filtered.drop('hgvs', axis=1)
    df_vcf_3_filtered = df_vcf_3_filtered.drop('hgvs', axis=1)

    # Write to file
    output_path = os.path.join(OUTPUT_PATH, f'{run_id}_{sample_id}_filtered_variants.tsv')
    with open(output_path, 'w+') as f:
        f.write('WARNING! This has not been validated. These variants would usually be filtered out by the pipeline, make sure you validate any findings using another method.\n')
        f.write('Variants are shown in reverse order of when they were filtered. For more info, see https://github.com/erikwaskiewicz/filter_potential_mosiacs/blob/master/README.md\n')
        f.write(f'Run:\t{run_id}\n')
        f.write(f'Sample:\t{sample_id}\n')
        f.write(f'Gene(s):\t{", ".join(genes)}\n')
        f.write(f'\nVCF filter 1 - {vcf_1.split("/")[-1]}\n')
        df_vcf_1.to_csv(f, index=False, sep='\t')
        f.write(f'\nVCF filter 2 - {vcf_2.split("/")[-1]}\n')
        df_vcf_2_filtered.to_csv(f, index=False, sep='\t')
        f.write(f'\nVCF filter 3 - {vcf_3.split("/")[-1]}\n')
        df_vcf_3_filtered.to_csv(f, index=False, sep='\t')


if __name__ == '__main__':
    # global settings
    settings = get_settings('config.yaml')
    BCFTOOLS_PATH = os.path.abspath(settings['bcftools_path'])
    OUTPUT_PATH = os.path.abspath(settings['output_path'])
    GENES = settings['genes']
    
    # run script
    main()
