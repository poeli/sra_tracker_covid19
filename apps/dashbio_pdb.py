import pandas as pd
import logging
from dash_bio_utils import pdb_parser, create_mol3d_style
from Bio.SeqUtils import seq1
from Bio import SearchIO

class SpikePdbData:
    """
    This is the class to process PDB file for spike protein
    """
    def __init__(self, pdb, blastxml):
        """
        blastxml: Query sequence is the reference spike protein sequence. Subject sequence is PDB sequences
        """
        self.pdb = pdb
        self.blastxml = blastxml
        self.pdb_data = None
        self.pdb_df = None
        self.blast_hsp = None
        self._process()

    def _process(self):
        """
        processing input files
        """
        # parsing pdb file for mol3d and pdb_style
        pdb_data = pdb_parser.PdbParser(self.pdb).mol3d_data()
        
        # create datatable
        df = pd.DataFrame(pdb_data['atoms'])
        df = df.drop_duplicates(subset=['residue_index'])
        df['residue_index_chain'] = df['residue_index']%len(df[df.residue_index=='A'])
        df['residue'] = df['residue_name'].str[0:3:1].map(seq1)
        
        # alignment
        blast_qresult = SearchIO.read(self.blastxml, "blast-xml")
        blast_hsp = blast_qresult[0][0]

        # rebuild residue index per chain for protein trimer
        chain_length = int(len(''.join(df.residue).rstrip('X'))/3)
        df['residue_index_chain'] = -1
        idx = df['residue']!='X'
        df.loc[idx, 'residue_index_chain'] = df.loc[idx, 'residue_index']%chain_length

        self.pdb_df = df
        self.blast_hsp = blast_hsp
        self.pdb_data = pdb_data
    
    def pdb_seq(self):
        """
        return all protein sequences (concat all chains) of the pdb file
        """
        return ''.join(self.pdb_df.residue.tolist())
    
    def pdb_style(self, residue_indexes=[], residue_indexes_color_map={}, highlight_bg_indexes=[], highlight_bg_indexes_color_map={}):
        styles = []

        for x in self.pdb_data['atoms']:
            rsd_idx = x['residue_index']

            if rsd_idx in residue_indexes_color_map:
                styles.append({'visualization_type': 'sphere', 'color': residue_indexes_color_map[rsd_idx]}),                        

            elif rsd_idx in residue_indexes:
                styles.append({'visualization_type': 'sphere', 'color': '#459DF8'})

            elif rsd_idx in highlight_bg_indexes_color_map:
                styles.append({'visualization_type': 'stick', 'color': highlight_bg_indexes_color_map[rsd_idx]})

            elif rsd_idx in highlight_bg_indexes:
                styles.append({'visualization_type': 'stick', 'color': '#E45E68'})

            else:
                styles.append({'visualization_type': 'line', 'color': '#AAAAAA'})
                    
        return styles

    
    def residue_labels(self, poss=[], label_text=[], font_size=12, background_opacity=0.8, poss_color_map={}):
        labels = []
        
        for pos, text in zip(poss, label_text):
            
            # setting background colors
            backgroundColor = poss_color_map[pos] if pos in poss_color_map else "#459DF8"
            fontColor = "black"
            
            # convert spike position to pdb_residue_index_chain
            residue_index = self.spike2pdb_residue_index_chain(pos)

            # retrieve atoms by residue_index_chain
            df = self.pdb_df.query(f'residue_index_chain=={residue_index}')

            for index, row in df.iterrows():
                label = {
                    "text": text,
                    "fontSize": font_size,
                    "fontColor": fontColor,
                    "backgroundColor": backgroundColor,
                    "backgroundOpacity": background_opacity,
                    "position": {
                        "x": row["positions"][0],
                        "y": row["positions"][1],
                        "z": row["positions"][2],
                    }
                }
                labels.append(label)
        
        return labels
    
    def spike2pdb_residue_index_chain(self, pos: int):
        (qs, qe) = self.blast_hsp.query_range
        (hs, he) = self.blast_hsp.hit_range
        
        aln_pos = pos - qs

        if aln_pos < 0:
            logging.info("position out of range.")
            return None
        else:
            # adjust alignment position if gaps are in the query seuqneces
            query_gaps = str(self.blast_hsp.query.seq)[:aln_pos].count('-')
            padding = 0
            while query_gaps != padding:
                padding = query_gaps
                aln_pos += padding
                query_gaps = str(self.blast_hsp.query.seq)[:aln_pos].count('-')
        
        # get corresponding target position
        target_base = str(self.blast_hsp.hit[aln_pos-1:aln_pos].seq)

        if target_base=='-':
            # corresponding base is missing from the pdb sequences
            # return None
            return None
        else:
            gaps = str(self.blast_hsp.hit[:aln_pos].seq).count('-')
            # The position need to -1 because pdb residue_index position is 0-basis
            return aln_pos-gaps+hs-1

    
    def spike2pdb_residue_indexes(self, poss=[]):
        """
        Pass in spike positions and get pdb_residue_indexes
        
        : param poss: list 
        : return residue_indexes: list 
        """
        residue_indexes = []
        
        for pos in poss:
            # get residue_index_per chain
            residue_index = self.spike2pdb_residue_index_chain(pos)
            # get residue_index series
            idxes = self.pdb_df.query(f'residue_index_chain=={residue_index}').residue_index.tolist()
            residue_indexes += idxes
        
        return residue_indexes
