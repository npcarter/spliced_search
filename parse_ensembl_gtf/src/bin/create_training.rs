use std::alloc::handle_alloc_error;
use std::fs::File;
use std::fs::exists;
use std::io::{BufRead,BufReader, Lines};
use std::result;
use glob::glob;
use noodles_fasta::record;
use std::iter::Peekable;
use clap::Parser;
use noodles_gtf as gtf;
use noodles_fasta as fasta;
use noodles_gff as gff;
use noodles_core;
use parse_ensembl_gtf::*;
use bstr::ByteSlice;

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
    let mut gene_gtf_records = Vec::new();

    // If we get here, we've opened both the GTF and FASTA files successfully
    for record_result in gtf_reader.record_bufs() {             
        let the_recbuf = record_result.expect("Unable to get GTF record buffer from file");

        if the_recbuf.ty() != "gene" {
            if in_protein_coding{
                gene_gtf_records.push(the_recbuf.clone());
                }
            continue; // Skip non-gene records
            }
        else{  // This is a gene record, make sure it's for a protein coding gene
            if in_protein_coding{ // We just finished reading the records for a protein coding gene
                        // so process it
                let the_gene = parse_gtf_gene(&gene_gtf_records, &current_fasta_record);
            }
            gene_gtf_records = Vec::new(); // clear this, since we're starting a new gene.
            let attrib = the_recbuf.attributes();
            match attrib.get(b"gene_biotype"){
                None =>{
                    in_protein_coding = false;
                    continue; 
                }
                Some(val) =>{
                    let the_value = val;
                    if the_value.as_string().expect("Couldn't get gene_biotype as string") != "protein_coding"{
                        in_protein_coding = false;
                        continue;
                    }
                    else{
                        in_protein_coding = true;
                        gene_gtf_records.push(the_recbuf.clone());
                    }
                }
            }
                
            // If we make it here without hitting a continue, the GTF record is the start of a gene with the protein_coding type,
            // so start building the vector of records that represent the gene
                
            // First, though, make sure we have the FASTA record that goes with the source of this gene
            while(current_fasta_record.name() != the_recbuf.reference_sequence_name()) {
                fasta_result = fasta_reader.records().next();
                match fasta_result {
                    Some(record) => {
                        current_fasta_record = record.expect("Failed to get record");
                    }
                    None => {
                        eprintln!("End of FASTA file reached before finding sequence for gene {}", the_recbuf.reference_sequence_name());
                        std::process::exit(1);
                    }
                }
            }

                    
        }   
    }
}   



// parse the GTF records that describe a gene into a structure
// requires that the gene is a protein coding gene
fn parse_gtf_gene(gtf_records:&Vec<gff::feature::RecordBuf> , fasta_record:&fasta::Record) -> Gene{
    let mut the_gene = Gene{
        gene_id: String::new(),
        start: 0,
        end: 0,
        strand: Strand::Unknown,
        proteins: Vec::new(),
    };

    let mut the_protein = Protein::default();

    for record in gtf_records {
        match record.ty().to_str() {
            Ok("gene")  => {
                let attrib = record.attributes();
                match attrib.get(b"gene_id") {
                    None => {
                        eprintln!("GTF gene entry gene_id attribute wasn't found");
                        the_gene.gene_id = "NONE".to_string();
                    }
                    Some(val) => {
                        assert!(attrib.get(b"gene_biotype").expect("Parse_gtf_gene couldn't get gene_biotype of gene record").as_string().expect("Couldn't get gene_biotype as string") == "protein_coding");
                        the_gene.gene_id = val.as_string().expect("GTF gene entry gene_id attribute wasn't a string").to_string();
                    }
                }
                match record.strand(){ // This is a bit ugly because of duelling Strand types between this code and noodles_gff
                    noodles_gff::feature::record::Strand::None => {                 
                        eprintln!("Gene record found with no strand");
                        std::process::exit(1);
                    }
                    noodles_gff::feature::record::Strand::Unknown => {                 
                        eprintln!("Gene record found with unknown strand");
                        std::process::exit(1);
                    }
                    noodles_gff::feature::record::Strand::Forward =>{
                        assert!(the_gene.strand!= Strand::Reverse); //Complain if gene seems to switch strands
                        the_gene.strand = Strand::Forward;
                    }
                    noodles_gff::feature::record::Strand::Reverse =>{
                        assert!(the_gene.strand != Strand::Forward); // Complain if gene seems to switch strands
                        the_gene.strand = Strand::Reverse;
                    }
                }
            }
            Ok("transcript") => {
            // start of new protein, so push old one onto gene list if there is one.
                if the_protein.name != None{ // There was a previous protein.  Push it onto the gene's list and create a new one
                    the_gene.proteins.push(the_protein);
                    the_protein = Protein::default();
                };
                // Grab what information we can from this record
                the_protein.mrna_start = record.start().get() ;
                the_protein.mrna_end = record.end().get();
                let start = noodles_core::Position::try_from(the_protein.mrna_start).expect("Couldn't generate Position from mrna_start");
                let end = noodles_core::Position::try_from(the_protein.mrna_end).expect("Couldn't generate position from mrna_end");


                // Fetch the DNA sequence from the FASTA record
                match the_gene.strand{
                    Strand::Unknown => {                 
                        eprintln!("Gene record had unknown strand after transcript line.");
                        std::process::exit(1);
                    }
                    Strand::Forward => {
                        the_protein.dna_sequence = Vec::from(fasta_record.sequence().get(start..=end).expect("Couldn't extract DNA sequence from forward strand transcript"));
                    }
                    Strand::Reverse => {
                        the_protein.dna_sequence = Vec::new(); // make sure this starts empty
                        for residue in fasta_record.sequence().slice(start..=end).expect("Couldn't extract DNA sequence from reverse strand transcript").complement(){
                            match residue{
                                Err(e) => {
                                    eprintln!("Error encountered complementing DNA sequence: {}", e);
                                    std::process::exit(1);
                                }
                                Ok(val) => {the_protein.dna_sequence.push(val);}
                            }
                        }
                    }
                }
            }
            Ok("exon") => {
                assert!(record.start().get() <= record.end().get()); // check to make sure our expectations about start and end are preserved
                // convert from 1-indexed positions within the FASTA record 
                // to 0-indexed offsets within the transcript
                // if reverse strand, handle the fact that we've reversed the transcript as part of complementing it
                let start:usize = match the_gene.strand{  
                    Strand::Forward => {
                        record.start().get() - the_protein.mrna_start
                    }
                    Strand::Reverse => {
                        the_protein.mrna_end -record.end().get()
                    }
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing exon start");
                    }
                };
                let end:usize = match the_gene.strand{
                    Strand::Forward => {
                        record.end().get() - the_protein.mrna_start
                    }
                    Strand::Reverse => {
                        the_protein.mrna_end -record.start().get()
                    }
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing exon start");
                    }
                };
                    record.start().get();

                let the_exon = Exon{
                    start: start,
                    end: end,
                };
                the_protein.exons.push(the_exon);
            }

            Err(e)=> {
                eprintln!("Unable to parse record type as string: {}", e);
                std::process::exit(1);
            }
            Ok(val) => {
                eprintln!("Unexpected record type field found: {}", val);
                std::process::exit(1);
            } 
        }
    }
    return the_gene;
}