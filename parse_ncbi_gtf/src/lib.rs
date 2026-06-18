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
    pub ncbi_protein: Vec<char>,  //DNA sequence of the protein
    pub translation_table: i8,
    pub alternative_start_codon: bool,
    pub start_codon: Vec<usize>,
    pub start_positions_found: usize,
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
            ncbi_protein: Vec::new(),
            translation_table: 0 ,
            alternative_start_codon: false,
            start_codon: vec!(0,0,0),
            start_positions_found: 0}
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
}

impl Default for NcbiSequence{
    fn default() -> Self{
        Self{
            id: "".to_string(),
            mrna_sequence: Vec::new(),
            protein_sequence: Vec::new(),
        }
    }
}

