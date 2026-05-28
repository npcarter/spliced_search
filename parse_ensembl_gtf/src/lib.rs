/* Shared types and functions for parse_ensembl_gtf */



pub struct Exon{
    pub start: i64,
    pub end: i64,
    pub complement: bool,
}

// For introns, need to know where they start and end in the contig.
// and whether they are complementary to the DNA sequence.
pub struct Intron{
    pub start: i64,
    pub end: i64,
}


pub struct Protein{
    pub name: String,
    pub uniprot_name: String,
    pub mrna_start: i64, // Start of the mRNA within the contig
    pub mrna_end: i64, // End of the mRNA within the contig
    pub coding_start: i64, //start position of the coding sequence within the mRNA
    pub coding_end: i64, // end position of the coding sequence within the mRNA
    pub introns: Vec<Intron>,
    pub exons: Vec<Exon>,
    pub codon_start: i8, // How many positions (0-2) is the first codon start offset from the
    // DNA sequence start
    pub dna_sequence: Vec<char>,  //DNA sequence of the protein
    pub protein_sequence: Vec<char>,
    pub translation_table: i8,
}

pub struct Gene{
    pub 
    pub proteins: Vec<Protein>,
}
