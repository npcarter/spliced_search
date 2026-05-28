use std::fs::File;
use std::fs::exists;
use std::io::{BufRead,BufReader, Lines};
use glob::glob;
use std::iter::Peekable;
use clap::Parser;
use noodles_gtf as gtf;
use noodles_fasta as fasta;
use parse_ensembl_gtf::*;


// Configure command-line arguments
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    gtffile: String,
    fastafile: String,
}


fn main() {
    let args = Cli::parse();
    let gtffile = args.gtffile;
    let fastafile = args.fastafile;

    // Open the GTF file and create a reader
    let file = File::open(gtffile.clone());
    let mut gtf_reader : gtf::io::Reader<BufReader<File>>;
    match file {
        Ok(file) => {
            gtf_reader = gtf::io::Reader::new(BufReader::new(file));
        }
        Err(e) => {
            eprintln!("Error opening GTF file {} : {}", gtffile, e);
            std::process::exit(1);
        }
    }

    // Open the FASTA file and create a reader
    let file = File::open(fastafile.clone());
    let mut fasta_reader : fasta::io::Reader<BufReader<File>>;
    match file {
        Ok(file) => {
            fasta_reader = fasta::io::Reader::new(BufReader::new(file));
        }
    
        Err(e) => {
            eprintln!("Error opening FASTA file {} : {}", fastafile, e);
            std::process::exit(1);
        }
    }

    let mut fasta_result = fasta_reader.records().next();
    let mut current_fasta_record : fasta::Record;
    match fasta_result {
        Some(record) => {
            //println!("Fasta record: {}", String::from_utf8(record.expect("Failed to get record").name().to_vec()).expect("Failed to convert record name to string"));
            current_fasta_record = record.expect("Failed to get record");
        }
        None => {
            eprintln!("Error reading FASTA file: No record found at start of file");
            std::process::exit(1)   
        }
    }
   

    // This loop is a bit convoluted because we don't know when we've found the end of 
    // a GTF gene until we hit the start of the next one.  So, we keep a state variable
    // (in_protein_coding) that records whether we're in GTF records for a protein coding gene.

    // Iterate through GTF records until we find the start of a gene.
    // If we're in a protein_coding gene, build a vector of the records that describe the gene.

    // When we find the start of a gene, if we're currently in a protein-coding gene, process that
    // into a gene record and handle it.  Otherwise, check whether the new gene is protein-coding
    // and make sure we have the FASTA entry that corresponds to it.

    let mut in_protein_coding = false;
    let mut gene_gtf_records: Vec<gtf::record> = Vec::new();

    for record_result in gtf_reader.record_bufs() {
        match record_result {
            Err(e) => {
                eprintln!("Error reading GTF file: {}", e);
                std::process::exit(1);
            }
            Ok(record) => { // skip until we find a gene record for a protein coding gene
                if record.ty() != "gene" {
                    if in_protein_coding{
                        gene_gtf_records.push(record.clone());
                    }
                    continue; // Skip non-gene records
                }
                else{  // This is a gene record, make sure it's for a protein coding gene
                    if in_protein_coding{ // We just finished reading the records for a protein coding gene
                        // so process it
                        let the_gene = parse_gtf_gene(gene_gtf_records, current_fasta_record);
                    }
                    gene_gtf_records = Vec::new(); // clear this, since we're starting a new gene.
                    let attrib = record.attributes();
                    match attrib.get(b"gene_biotype"){
                        None =>{
                            in_protein_coding = false;
                            continue; 
                        }
                        Some(val) =>{
                            match val.as_string(){
                                None => {
                                    eprintln!("GTF gene entry gene_biotype attribute wasn't a string");
                                    std::process::exit(1);
                                }
                                Some(the_str) => {
                                    if the_str != b"protein_coding"{
                                        in_protein_coding = false;
                                        continue;
                                    } 
                                    else{
                                        in_protein_coding = true;
                                        gene_gtf_records.push(record.clone());
                                    }
                                }
                            }
                        }
                    }
                
                    // If we make it here without hitting a continue, the GTF record is the start of a gene with the protein_coding type,
                    // so we want to parse it
                
                    // First, though, make sure we have the FASTA record that goes with the source of this gene
                    while(current_fasta_record.name() != record.reference_sequence_name()) {
                        fasta_result = fasta_reader.records().next();
                        match fasta_result {
                            Some(record) => {
                                current_fasta_record = record.expect("Failed to get record");
                            }
                            None => {
                                eprintln!("End of FASTA file reached before finding sequence for gene {}", record.reference_sequence_name());
                                std::process::exit(1);
                            }
                        }
                    }

                    
                }   
            }
        }
    }   
}

fn parse_gtf_gene(gtf_records:Vec<gtf::record>, fasta_record:fasta::record) => Gene{
    
}