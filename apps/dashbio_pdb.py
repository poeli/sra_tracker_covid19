import pandas as pd
import logging
from dash_bio_utils import pdb_parser
from Bio.SeqUtils import seq1
from Bio import SearchIO

class SpikePdbData:
    """
    This is the class to process PDB file for spike protein
    """
    def __init__(self, pdb:str, blastxml:str):
        """
        : param pdb : A .pdb file
        : param blastxml: Query sequence is the reference spike protein sequence. Subject sequence is PDB sequences
        """
        self.pdb = pdb
        self.blastxml = blastxml
        self.pdb_data = None
        self.pdb_df = None
        self.blast_qresult = None
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
        blast_qresult = SearchIO.read(self.blastxml, "blast-xml")[0]
        # blast_hsp = blast_qresult

        self.pdb_df = df
        self.blast_qresult = blast_qresult
        self.pdb_data = pdb_data
    

    def pdb_seq(self) -> str:
        """
        return all protein sequences (concat all chains) of the pdb file
        """
        return ''.join(self.pdb_df.residue.tolist())
    

    def pdb_style(self, residue_indexes: list=[], 
                        residue_indexes_color_map: dict={}, 
                        highlight_bg_indexes: list=[], 
                        highlight_bg_indexes_color_map: dict={},
                        highlight_bg_chain: list=[],
                        highlight_bg_chain_color_map: dict={}) -> list:
        """
        return a list of dict for 3DMol
        """
        styles = []

        for x in self.pdb_data['atoms']:
            # get residue index
            ridx = x['residue_index']
            chain = x['chain']

            if ridx in residue_indexes_color_map:
                styles.append({'visualization_type': 'sphere', 'color': residue_indexes_color_map[ridx]}),                        
            elif ridx in residue_indexes:
                # coloring blue
                styles.append({'visualization_type': 'sphere', 'color': '#459DF8'})
            elif ridx in highlight_bg_indexes_color_map:
                styles.append({'visualization_type': 'stick', 'color': highlight_bg_indexes_color_map[ridx]})
            elif ridx in highlight_bg_indexes:
                # coloring red
                styles.append({'visualization_type': 'stick', 'color': '#E68E96'})
            elif chain in highlight_bg_chain_color_map:
                styles.append({'visualization_type': 'stick', 'color': highlight_bg_chain_color_map[chain]})
            elif chain in highlight_bg_chain:
                # coloring green
                styles.append({'visualization_type': 'stick', 'color': '#73C991'})
            else:
                # coloring grey
                styles.append({'visualization_type': 'line', 'color': '#AAAAAA'})
                    
        return styles

    
    def residue_labels(self, poss: list=[], 
                             label_text: list=[], 
                             font_size: int=12, 
                             background_opacity: float=0.8, 
                             poss_color_map: dict={}) -> list:
        labels = []
        
        for pos, text in zip(poss, label_text):
            
            # setting background colors
            backgroundColor = poss_color_map[pos] if pos in poss_color_map else "#459DF8"
            fontColor = "black"
            
            # convert spike position to pdb_residue_index_chain
            residue_indexes = self.spike_pos2pdb_residue_index([pos])

            # retrieve atoms by residue_index_chain
            dfidx = self.pdb_df.residue_index.isin(residue_indexes)
            df = self.pdb_df[dfidx]

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
    

    def spike_pos2pdb_residue_index(self, poss:list=[]) -> list:

        residue_indexes = []

        for pos in poss:
            # iterate all blast hsps, it will convert spike position to the positions of all possible chains
            for idx, blast_hsp in enumerate(self.blast_qresult):
                (qs, qe) = blast_hsp.query_range
                (hs, he) = blast_hsp.hit_range
                
                aln_pos = pos - qs

                if aln_pos < 0:
                    logging.debug(f"S position {aln_pos} is out of range on HSP{idx}.")
                    # residue_indexes.append(None)
                    continue
                else:
                    # adjust alignment position if gaps are in the query seuqneces
                    query_gaps = str(blast_hsp.query.seq)[:aln_pos].count('-')
                    padding = 0
                    while query_gaps != padding:
                        padding = query_gaps
                        aln_pos += padding
                        query_gaps = str(blast_hsp.query.seq)[:aln_pos].count('-')
                
                # get corresponding target position
                target_base = str(blast_hsp.hit[aln_pos-1:aln_pos].seq)

                if target_base=='-':
                    # corresponding base is missing from the pdb sequences
                    logging.debug(f"S position {aln_pos} is a gap on HSP{idx}.")
                    # residue_indexes.append(None)
                    continue
                else:
                    gaps = str(blast_hsp.hit[:aln_pos].seq).count('-')
                    # The position need to -1 because pdb residue_index position is 0-basis
                    residue_indexes.append(aln_pos-gaps+hs-1)

        return residue_indexes