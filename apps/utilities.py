import pandas as pd
import numpy as np
import logging
import glob
import sys
import io

class EC19_project():
    '''
    This is the class to handle the outputs generated from EDGE-Covid19 workflow.
    '''
    samples     = []
    df_meta     = pd.DataFrame()
    df_snps     = pd.DataFrame()
    df_indels   = pd.DataFrame()
    df_alnstats = pd.DataFrame()
    df_pango    = pd.DataFrame()
    df_gaps     = pd.DataFrame()
    df_qc       = pd.DataFrame()
    
    def __init__(self, sample_path, store_in_class=True, error_on_duplication=True):
        # input files
        self.sample_path = sample_path # last project path
        self.config      = {}
        self.df_meta     = pd.DataFrame()
        self.df_snps     = pd.DataFrame()
        self.df_indels   = pd.DataFrame()
        self.df_alnstats = pd.DataFrame()
        self.df_pango    = pd.DataFrame()
        self.df_gaps     = pd.DataFrame()
        self.df_qc       = pd.DataFrame()
        
        self.store_in_class = store_in_class
        self.error_on_duplication = error_on_duplication
        
    def add_project(self, sample_name=None):
        '''
        add project
        '''
        # init sample name
        if not sample_name:
            # parsing config file
            file = f'{self.sample_path}/config.txt'
            self.config = self.parse_config(file)
            sample = self.get_project_name()
            logging.info(f'Project name found: {sample}')
        
        if sample in EC19_project.samples:
            # skip duplication
            logging.fatal(f'Duplication found: {sample}.')
            
            if self.error_on_duplication:
                sys.exit(1)
            else:
                logging.fatal(f'Skipping to store to class level.')
                self.store_in_class = False
        else:
            EC19_project.samples.append(sample)
        
        # parse project
        self.parse_project(sample_name)
        
        if self.store_in_class:
            EC19_project.df_snps = pd.concat([EC19_project.df_snps, self.df_snps]).reset_index(drop=True)
            EC19_project.df_indels = pd.concat([EC19_project.df_indels, self.df_indels]).reset_index(drop=True)
            EC19_project.df_alnstats = pd.concat([EC19_project.df_alnstats, self.df_alnstats]).reset_index(drop=True)
            EC19_project.df_gaps = pd.concat([EC19_project.df_gaps, self.df_gaps]).reset_index(drop=True)
            EC19_project.df_meta = pd.concat([EC19_project.df_meta, self.df_meta]).reset_index(drop=True)
            EC19_project.df_pango = pd.concat([EC19_project.df_pango, self.df_pango]).reset_index(drop=True)
            EC19_project.df_qc = pd.concat([EC19_project.df_qc, self.df_qc]).reset_index(drop=True)
        
    def parse_project(self, sample_name=None):
        '''
        parse project
        '''
        self.df_snps = self.parse_snps()
        self.df_indels = self.parse_indels()
        self.df_alnstats = self.parse_alnstats()
        self.df_gaps = self.parse_gaps()
        self.df_meta = self.parse_metadata()
        self.df_pango = self.parse_pango()
        self.df_qc = self.parse_qc_stats()

        
    def parse_config(self, file=None):
        '''
        parse EDGE config.txt
        '''
        import re
        
        config = {}
        section = 'main'
        
        with open(file) as config_fh:
            for line in config_fh:
                line = line.strip()
                if line and not line.startswith('#'):
                    if line.startswith('['):
                        section = re.search('^\[(.+)\]$', line).group(1)
                        config[section] = {}
                    else:
                        try:
                            (key, value) = line.split('=')
                            config[section][key] = value
                        except:
                            print(line)

        return config
        
    def parse_metadata(self, file=None):
        '''
        Adding EC19 metadata
        '''
        import glob
        
        files = glob.glob(f'{self.sample_path}/metadata_*.txt')
        if file:
            files = [file]

        metadata = {}
        dfs = pd.DataFrame()
        sample = self.get_project_name()
  
        try:
            for file in files:
                with open(file) as fh:
                    for line in fh:
                        line = line.strip()
                        (key, value) = line.split('=')
                        metadata[key] = value

                logging.info(f'{sample}: EC19 metadata file ({file}) loaded.')

            df = pd.DataFrame(metadata, index=[0])
            df['SAMPLE'] = sample
            
            dfs = pd.concat([dfs, df])
        except:
            logging.info(f'{sample}: EC19 metadata file not loaded.')
            
        return dfs

    def parse_qc_stats(self, file=None):
        '''
        Adding EC19 dataset info from QC stats
        '''
        import glob
        import re
        
        files = glob.glob(f'{self.sample_path}/QcReads/QC.stats.txt')

        qc_stats = {}
        sample = self.get_project_name()
        prefix = ''
        df = pd.DataFrame()

        try:
            for file in files:
                with open(file) as fh:
                    for line in fh:
                        line = line.strip()

                        if 'After Trimming' in line:
                            prefix = 'TRIMMED_'
                        if line.startswith('Reads #'):
                            val = re.search('Reads #: (\d+)', line).group(1)
                            qc_stats[f'{prefix}READ_NUM'] = int(val)
                        if line.startswith('Total bases'):
                            val = re.search('Total bases: (\d+)', line).group(1)
                            qc_stats[f'{prefix}TOL_BASES'] = int(val)
                        if line.startswith('Reads Length'):
                            val = re.search('Reads Length: (\S+)', line).group(1)
                            qc_stats[f'{prefix}READ_LEN'] = float(val)

                df = pd.DataFrame(qc_stats, index=[0])
                df['SAMPLE'] = sample
                logging.info(f'{sample}: EC19 QC stats file ({file}) loaded.')
        except Exception as e:
            logging.info(f'{sample}: Exception: {e}')
            
        return df

    
    def parse_indels(self, file=None):
        ''' 
        Parsing EC19 indel report (NC_045512.2_consensus.Indels_report.txt)
        '''
        
        files = []
        indel_files = glob.glob(f'{self.sample_path}/ReadsBasedAnalysis/readsMappingToRef/*.Indels_report.txt')

        # for file in indel_files:
        #     # keep readToRef.alnstats.txt
        #     if file.endswith('readsToRef.Indels_report.txt'):
        #         files.append(file)

        if not files:
            files = indel_files          

        dfs = pd.DataFrame()
        sample = self.get_project_name()
        
        try:
            for file in files:
                logging.debug(f'{sample}: EC19 Indel file ({file}) processing...')
                df = pd.read_csv(file, sep='\t')

                # skip empty file
                if len(df)==0:
                    logging.info(f'{sample}: Skip processing Indel file.')
                    # concat headers
                    if 'Chromosome' in df:
                        dfs = pd.concat([dfs, df])
                    continue
                
                # converting product names to abbreviations
                df.loc[:,'Product'] = df.loc[:,'Product'].str.lower()
                df.loc[:,'Product'] = df.loc[:,'Product'].replace(regex=[r'^\S+\d+:', ' \S*protein'], value='')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('nucleocapsid', 'N')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('surface', 'S')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('membrane', 'M')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('envelope', 'E')
                df.loc[:,'Product'] = df.loc[:,'Product'].str.replace('orf', 'ORF', case=True)

                idx = ~df['CDS_start'].isna()
                df.loc[idx, 'CDS_start']    = df.loc[idx, 'CDS_start'].astype(int)
                df.loc[idx, 'CDS_end']      = df.loc[idx, 'CDS_end'].astype(int)

                df['SAMPLE'] = sample
                
                dfs = pd.concat([dfs, df])     
                
                logging.info(f'{sample}: EC19 indels file ({file}) loaded.')
        except Exception as e:
            logging.info(f'{sample}: Exception: {e}')
        
        return dfs

    
    def parse_snps(self, file=None, intergenic=True, synonymous=True):
        ''' 
        Parsing EC19 SNP report (*_consensus.SNPs_report.txt)
        '''
        
        files = glob.glob(f'{self.sample_path}/ReadsBasedAnalysis/readsMappingToRef/readsToRef.SNPs_report.txt')
        
        if len(files)==0:
            files += glob.glob(f'{self.sample_path}/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.SNPs_report.txt')

        if file:
            files = [file]
        
        dfs = pd.DataFrame()
        sample = self.get_project_name()
        
        try:
            for file in files:
                logging.debug(f'{sample}: EC19 SNPs file ({file}) processing...')
                df = pd.read_csv(file, sep='\t')

                # skip empty file
                if len(df)==0:
                    logging.info(f'{sample}: Skip empty file.')
                    continue
                
                # skip duplicated references
                if len(dfs)>0 and len(df)>0:
                    if df.iloc[0,0] in dfs['Chromosome']:
                        logging.info(f'{sample}: Skipped duplicated reference.')
                        continue

                # add sample name
                df['SAMPLE'] = sample
                        
                # dropping merged snp
                df = df.dropna(subset=['Product'])
                        
                # loading dropping mutations in intergenic regions
                if not intergenic:
                    df = df.dropna(subset=['CDS_start'])    
                
                # remove synonymous changes
                if not synonymous:
                    df = df.drop(df[df.Synonymous=='Yes'].index)
                
                # converting product names to abbreviations
                df.loc[:,'Product'] = df.loc[:,'Product'].str.lower()
                df.loc[:,'Product'] = df.loc[:,'Product'].replace(regex=[r'^\S+\d+:', ' \S*protein'], value='')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('nucleocapsid', 'N')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('surface', 'S')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('membrane', 'M')
                df.loc[:,'Product'] = df.loc[:,'Product'].replace('envelope', 'E')
                df.loc[:,'Product'] = df.loc[:,'Product'].str.replace('orf', 'ORF', case=True)
                df['SNP_position']  = df['SNP_position'].astype(int)
                df['AA_pos']        = np.nan
                
                # Non-intergenic region
                idx = ~df.CDS_start.isna()                
                df.loc[idx, 'CDS_start'] = df.loc[idx, 'CDS_start'].astype(int)
                df.loc[idx, 'CDS_end']   = df.loc[idx, 'CDS_end'].astype(int)
                df.loc[idx, 'AA_pos']    = np.int16((df.loc[idx, 'SNP_position']-df.loc[idx, 'CDS_start'])/3)+1
                
                # Non-Synonymous only
                idx = df.Synonymous=='No'
                
                df.loc[idx, 'Mutation'] = df.loc[idx, 'Product'] + ':' + \
                                                df.loc[idx, 'aa_Ref'] + \
                                                df.loc[idx, 'AA_pos'].astype(str).str.replace(r'\.0$', '', regex=True) + \
                                                df.loc[idx, 'aa_Sub']

                # set SNP type
                df['Type'] = 'Intergenic region'
                df.loc[df.Synonymous=='Yes', 'Type'] = 'Synonymous'
                df.loc[df.Synonymous=='No', 'Type'] = 'Non-synonymous'
                
                dfs = pd.concat([dfs, df])            
                
                logging.info(f'{sample}: EC19 SNPs file ({file}) loaded.')

            # add missing columns if needed
            if not 'Count' in dfs:
                dfs['Count'] = np.nan
            if not 'F:R' in dfs:
                dfs['F:R'] = np.nan

        except Exception as e:
            logging.info(f'{sample}: Exception: {e}')
        
        return dfs
        
    def parse_alnstats(self, file=None):
        '''
        Parsing EC19 alnstats report (*.alnstats.txt)
        '''
        import re

        files = []
        alnstats_files = glob.glob(f'{self.sample_path}/ReadsBasedAnalysis/readsMappingToRef/*.alnstats.txt')
        # file = f'{self.sample_path}/ReadsBasedAnalysis/readsMappingToRef/readToRef.alnstats.txt'

        for file in alnstats_files:
            # keep readToRef.alnstats.txt
            if file.endswith('readsToRef.alnstats.txt'):
                files.append(file)

        if not files:
            files = alnstats_files          

        if file:
            files = [file]

        dfs = pd.DataFrame()
        sample = self.get_project_name()

        try:        
            for file in files:

                mapped_pct = 0
                content = ''
                content_flag = False

                logging.debug(f'{sample}: EC19 alnstats file ({file}) processing...')
                
                with open(file) as fh:
                    for line in fh:
                        if 'mapped (' in line:
                            mapped_pct = re.search('mapped \(([^%]+)% :', line).group(1)
                            mapped_pct = float(mapped_pct)

                        if line.startswith('Ref'):
                            content_flag = True

                        if content_flag:
                            content += line

                if len(content)>0:
                    df = pd.read_csv(io.StringIO(content), sep='\t', index_col=False)
                    
                    # skip empty file
                    if len(df)==0:
                        logging.info(f'{sample}: Skip empty file.')
                        continue
                    
                    df = df.rename(columns={'Ref_GC%':'Ref_GC_pct',  'Ref_recovery%': 'Ref_recovery_pct',  'Avg_fold(x)': 'Avg_fold_x'})
                    df = df.astype({'Ref_GC_pct': 'float', 
                                    'Mapped_reads': 'int64', 
                                    'Ref_recovery_pct': 'float', 
                                    'Avg_fold_x': 'float',
                                    'Fold_std': 'float', 
                                    'Num_of_Gap': 'int64',
                                    'Total_Gap_bases': 'int64'})
                    # only 
                    if 'Num_of_SNPs' in df:
                        df['Num_of_SNPs'] = df['Num_of_SNPs'].astype('int64')
                    elif len(self.df_snps):
                        df['Num_of_SNPs'] = len(self.df_snps)

                    if 'Num_of_INDELs' in df:
                        df['Num_of_INDELs'] = df['Num_of_INDELs'].astype('int64')
                    elif len(self.df_indels):
                        df['Num_of_INDELs'] = len(self.df_indels)

                    df['Mapped_reads_pct'] = mapped_pct
                    df['SAMPLE'] = sample

                    dfs = pd.concat([dfs, df])

                logging.info(f'{sample}: EC19 alnstats file ({file}) loaded.')
        except Exception as e:
            logging.info(f'{sample}: Exception: {e}')
        
        return dfs

    def parse_pango(self, file=None):
        ''' 
        Parsing EC19 pango lineage (NC_045512.2_consensus_lineage.txt)
        '''
        
        if not file:
            file = f'{self.sample_path}/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus_lineage.txt'
        
        df_pango = pd.DataFrame()
        sample = self.get_project_name()
        
        try:
            df_pango = pd.read_csv(file)
            df_pango.loc[:,'scorpio_call'] = df_pango.loc[:,'scorpio_call'].replace(regex=r' .*$', value='')
            df_pango['variant'] = df_pango['lineage']
            idx = df_pango['scorpio_call'].notnull()
            df_pango.loc[idx, 'variant'] = df_pango.loc[idx, 'scorpio_call']
            df_pango['SAMPLE'] = sample
            logging.info(f'{sample}: EC19 pango file ({file}) loaded.')
        except:
            logging.info(f'{sample}: EC19 pango file ({file}) not loaded.')
            
        return df_pango
              
    def parse_gaps(self, file=None):
        '''
        Parsing EC19 gaps report (*_consensus.gaps_report.txt)
        '''
        files = glob.glob(f'{self.sample_path}/ReadsBasedAnalysis/readsMappingToRef/*_consensus.gaps_report.txt')
        if file:
            files = [file]

        dfs = pd.DataFrame()
        sample = self.get_project_name()

        try:
            for file in files:
                logging.debug(f'{sample}: EC19 gaps file ({file}) processing...')
                df = pd.read_csv(file, sep='\t')
                df = df.replace(r'^\s*$', 0, regex=True)
                df = df.astype({'cds_start': 'int64', 
                                'cds_end': 'int64',
                                'gap_start': 'int64', 
                                'gap_end': 'int64'
                               })

                # skip empty file
                if len(df)==0:
                    logging.info(f'{sample}: Skip empty file.')
                    continue
                
                df['gap_aa_start'] = np.int64((df['gap_start']-df['cds_start'])/3)+1
                df['gap_aa_end'] = np.int64((df['gap_end']-df['cds_start'])/3)+1

                idx = df.cds_start==0
                df.loc[idx, 'gap_aa_start'] = np.nan
                df.loc[idx, 'gap_aa_end'] = np.nan

                df = df.drop(columns=['cds_start', 'cds_end', 'cds_product'])

                idx = df['gap_aa_start']<=0
                df.loc[idx, 'gap_aa_start'] = 1
                df['SAMPLE'] = sample

                dfs = pd.concat([dfs, df])

                logging.info(f'{sample}: EC19 gaps file ({file}) loaded.')
        except Exception as e:
            logging.info(f'{sample}: Exception: {e}')
        
        return dfs
    

    def get_project_name(self):
        ''' get project name '''
        try:
            return self.config['project']['projname']
        except:
            return ''

    def get_sample_df_snps(self, sample=None):
        ''' get df_snps by sample name'''
        sample = sample if sample else self.get_project_name()
        df = EC19_project.df_snps
        return df[df.SAMPLE==sample]
        
    def get_sample_df_alnstats(self, sample=None):
        ''' get df_alnstats by sample name'''
        sample = sample if sample else self.get_project_name()
        df = EC19_project.df_alnstats
        return df[df.SAMPLE==sample]
        
    def get_sample_df_pango(self, sample=None):
        ''' get df_snps by sample name'''
        sample = sample if sample else self.get_project_name()
        df = EC19_project.df_pango
        return df[df.SAMPLE==sample]

    def get_sample_df_gaps(self, sample=None):
        ''' get df_gaps by sample name'''
        sample = sample if sample else self.get_project_name()
        df = EC19_project.df_gaps
        return df[df.SAMPLE==sample]

    def get_sample_df_meta(self, sample=None):
        ''' get df_meta by sample name'''
        sample = sample if sample else self.get_project_name()
        df = EC19_project.df_meta
        return df[df.SAMPLE==sample]

    def get_sample_df_qc(self, sample=None):
        ''' get df_qc by sample name'''
        sample = sample if sample else self.get_project_name()
        df = EC19_project.df_qc
        return df[df.SAMPLE==sample]