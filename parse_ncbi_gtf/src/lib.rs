/* Shared types and functions for parse_ensembl_gtf */
use serde::Deserialize;
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
// and whether they are complementary to the DNA sequence.
pub struct Intron{
    pub start: i64,
    pub end: i64,
}


pub struct Protein{
    pub name: Option<String>,
    pub uniprot_name: String,
    pub mrna_start: usize, // Start of the mRNA within the contig
    pub mrna_end: usize, // End of the mRNA within the contig
    pub coding_start: usize, //start position of the coding sequence within the mRNA
    pub coding_end: usize, // end position of the coding sequence within the mRNA
    pub introns: Vec<Intron>,
    pub exons: Vec<Exon>,
    pub codon_start: i8, // How many positions (0-2) is the first codon start offset from the
    // DNA sequence start
    pub dna_sequence: Vec<u8>,  //DNA sequence of the protein
    pub translation_table: i8,
}

impl Default for Protein{
    fn default() -> Self{
        Self { name: None, // name == None is how we detect uninitialized proteins
            uniprot_name: "".to_string(),
            mrna_start: 0, 
            mrna_end: 0, 
            coding_start: 0, 
            coding_end: 0, 
            introns: Vec::new(), 
            exons: Vec::new(), 
            codon_start: 0, 
            dna_sequence: Vec::new(), 
            translation_table: 0 }
    }
}

pub struct Gene{
    pub gene_id: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    pub proteins: Vec<Protein>,
}

#[derive(Deserialize, Debug)]
struct EnsemblSequence {
    seq: String,
    id: String,
    molecule: String,
}