use std::u8;

/* Shared types and functions for parse_ensembl_gtf */
// Enum that records which strand a gene is on
#[derive(PartialEq)]
pub enum Strand{
    Forward,
    Reverse,
    Unknown,
}

#[derive(PartialEq)]
pub enum CodonPosition{
    First,
    Second,
    Third,
}
pub struct Exon{
    pub start: usize,
    pub end: usize,
    pub start_offset: Option<CodonPosition>,
}
// For introns, need to know where they start and end in the contig.
// and whether they are complementary to the DNA sequence
pub struct Protein{
    pub name: Option<String>,
    pub mrna_start: usize, // Start of the mRNA within the contig
    pub mrna_end: usize, // End of the mRNA within the contig
    pub coding_start: usize, //start position of the coding sequence within the mRNA
    pub coding_end: usize, // end position of the coding sequence within the mRNA
    pub exons: Vec<Exon>,
    pub codon_start: i8, // How many positions (0-2) is the first codon start offset from the
    // DNA sequence start
    pub dna_sequence: Vec<char>,
    pub ncbi_mrna: Vec<char>,
    pub computed_coding: Vec<char>,
    pub computed_protein: Vec<char>,
    pub ncbi_protein: Vec<char>,  //DNA sequence of the protein
    pub translation_table: i8,
    pub alternative_start_codon: bool,
    pub start_codon: Vec<usize>,
    pub start_positions_found: usize,
    pub stop_codon: Vec<usize>,
    pub stop_positions_found: usize,
    pub wrong_transcript: bool, // record whether the protein's note says it does not match the transcript
    pub selenocystine: bool,
    pub isoform_label: String, // text extracted describing the isoform of the protein.  
    pub transcript_label: String,
    // This gets gathered early, and then used later, once we have information about all the proteins in a gene.
}

impl Default for Protein{
    fn default() -> Self{
        Self {
            name: None, // name == None is how we detect uninitialized proteins
            mrna_start: std::usize::MAX, 
            mrna_end: 0, 
            coding_start: std::usize::MAX, 
            coding_end: 0, 
            exons: Vec::new(), 
            codon_start: 0, 
            dna_sequence: Vec::new(), 
            ncbi_mrna: Vec::new(),
            computed_coding: Vec::new(),
            computed_protein: Vec::new(),
            ncbi_protein: Vec::new(),
            translation_table: 0 ,
            alternative_start_codon: false,
            start_codon: vec!(0,0,0),
            start_positions_found: 0,
            stop_codon: vec!(0,0,0),
            stop_positions_found: 0,
            wrong_transcript: false, 
            selenocystine: false,
            isoform_label: "".to_string(),
            transcript_label: "".to_string(),
        }
    }
}

pub struct Gene{
    pub gene_id: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    pub proteins: Vec<Protein>,
}

// Data structure to hold either a protein or an mRNA sequence
// from the NCBI datasets
#[derive(PartialEq)]
pub struct NcbiSequence{
    pub id: String,
    pub mrna_sequence: Vec<char>,
    pub protein_sequence: Vec<char>,
    pub already_found: bool,
}

impl Default for NcbiSequence{
    fn default() -> Self{
        Self{
            id: "".to_string(),
            mrna_sequence: Vec::new(),
            protein_sequence: Vec::new(),
            already_found: false,
        }
    }
}

pub fn translate_codon(codon: Vec<char>, selenocystine: bool) -> char{
    if codon.len() != 3 { panic!("Codon must be three characters long"); }

    let amino_acid:char;
    let codon_string = codon.iter().collect::<String>();
    match codon_string.as_str(){
        "AAA" => amino_acid = 'K',
        "AAC" => amino_acid = 'N',
        "AAG" => amino_acid = 'K',
        "AAT" => amino_acid = 'N',
        "ACA" => amino_acid = 'T',
        "ACC" => amino_acid = 'T',
        "ACG" => amino_acid = 'T',
        "ACT" => amino_acid = 'T',
        "AGA" => amino_acid = 'R',
        "AGC" => amino_acid = 'S',
        "AGG" => amino_acid = 'R',
        "AGT" => amino_acid = 'S',
        "ATA" => amino_acid = 'I',
        "ATC" => amino_acid = 'I',
        "ATG" => amino_acid = 'M',
        "ATT" => amino_acid = 'I',
        "CAA" => amino_acid = 'Q',
        "CAC" => amino_acid = 'H',
        "CAG" => amino_acid = 'Q',
        "CAT" => amino_acid = 'H',
        "CCA" => amino_acid = 'P',
        "CCC" => amino_acid = 'P',
        "CCG" => amino_acid = 'P',
        "CCT" => amino_acid = 'P',
        "CGA" => amino_acid = 'R',
        "CGC" => amino_acid = 'R',
        "CGG" => amino_acid = 'R',
        "CGT" => amino_acid = 'R',
        "CTA" => amino_acid = 'L',
        "CTC" => amino_acid = 'L',
        "CTG" => amino_acid = 'L',
        "CTT" => amino_acid = 'L',
        "GAA" => amino_acid = 'E',
        "GAC" => amino_acid = 'D',
        "GAG" => amino_acid = 'E',
        "GAT" => amino_acid = 'D',
        "GCA" => amino_acid = 'A',
        "GCC" => amino_acid = 'A',
        "GCG" => amino_acid = 'A',
        "GCT" => amino_acid = 'A',
        "GGA" => amino_acid = 'G',
        "GGC" => amino_acid = 'G',
        "GGG" => amino_acid = 'G',
        "GGT" => amino_acid = 'G',
        "GTA" => amino_acid = 'V',
        "GTC" => amino_acid = 'V',
        "GTG" => amino_acid = 'V',
        "GTT" => amino_acid = 'V',
        "TAA" => amino_acid = '*',
        "TAC" => amino_acid = 'Y',
        "TAG" => amino_acid = '*',
        "TAT" => amino_acid = 'Y',
        "TCA" => amino_acid = 'S',
        "TCC" => amino_acid = 'S',
        "TCG" => amino_acid = 'S',
        "TCT" => amino_acid = 'S',
        "TGA" => amino_acid = match selenocystine{
            true => 'U',
            false => '*',
        },
        "TGC" => amino_acid = 'C',
        "TGG" => amino_acid = 'W',
        "TGT" => amino_acid = 'C',
        "TTA" => amino_acid = 'L',
        "TTC" => amino_acid = 'F',
        "TTG" => amino_acid = 'L',
        "TTT" => amino_acid = 'F',
        _ => amino_acid = 'X',
    }

    amino_acid
}

//FSM states
#[derive(Clone, PartialEq, Debug)]
pub enum FsmState{
    Intergenic, 
    Start,
    LastStart, //special case for when the first coding nucleotide is the last one in an exon
    Stop,
    FirstStop, // special case to handle situation where the first nucleotide in an exon is 
    // the last nucleotide in a protein.  Doesn't affect the training data, just the 
    // internal bookkeeping and self-checks
    Exon0,
    Exon1,
    Exon2,
    Intron0,
    Intron1,
    Intron2,
    Ass0,
    Ass1,
    Ass2,
    Dss0,
    Dss1,
    Dss2,
}

//Function to generate numeric representation for an FSM State
pub fn state_to_u8(state: FsmState) -> u8{
    match state{
        FsmState::Intergenic=> {'0' as u8}, 
        FsmState::Start | FsmState::LastStart => {'1' as u8},
        FsmState::Stop | FsmState::FirstStop=> {'2' as u8}, 
        // all stop variants generate the same training state, the differences are just
        // for sanity checks
        FsmState::Exon0=> {  '3' as u8},
        FsmState::Exon1=> {'4' as u8},
        FsmState::Exon2=> {'5' as u8},
        FsmState::Intron0=> {'6' as u8},
        FsmState::Intron1=> {'7' as u8},
        FsmState::Intron2=> {'8' as u8},
        FsmState::Ass0=>{'9' as u8},
        FsmState::Ass1=>{'A' as u8},
        FsmState::Ass2=>{'B' as u8},
        FsmState::Dss0=>{'C' as u8},
        FsmState::Dss1=>{'D' as u8},
        FsmState::Dss2=>{'E' as u8},
    }
}

pub const EXONSTATES: &[FsmState] = &[FsmState::Exon0, FsmState::Exon1, FsmState::Exon2];

pub const INTRONSTATES: &[FsmState] = &[FsmState::Intron0, FsmState::Intron1, FsmState::Intron2];

// These match the exon states because Dss replaces the last coding position in an exon
pub const DSSSTATES: &[FsmState] = &[FsmState::Dss0, FsmState::Dss1, FsmState::Dss2];

// These match the exon states because Ass replaces the first coding position in an exon
pub const ASSSTATES: &[FsmState] = &[FsmState::Ass0, FsmState::Ass1, FsmState::Ass2];