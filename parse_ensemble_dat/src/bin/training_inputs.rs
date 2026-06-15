// parse_ensemble_dat: takes an ensemble dat file that describes a genome and
// parses it into DNA chunks that contain proteins along with vectors that
// describe their introns, exons, start points, and end points.


use std::fs::File;
use std::fs::exists;
use std::io::{Write};
use std::io::{BufRead,BufReader, Lines};
use glob::glob;
use std::iter::Peekable;
use clap::Parser;
use std::cmp::{min, max};

use parse_ensemble_dat::*;

// define constants for the HMM states
const ILLEGAL_STATE: char = 'Z';
const INTERGENIC: char='0';
const START: char = '1';
const STOP: char = '2';
const INTRON0: char = '3';
const INTRON1: char = '4';
const INTRON2: char = '5';
const EXON0: char = '6';
const EXON1: char = '7';
const EXON2: char = '8';
const ASS0: char = '9';
const ASS1: char = 'A';
const ASS2: char = 'B';
const DSS0: char = 'C';
const DSS1: char = 'D';
const DSS2: char = 'E';




#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    sourcedir: String,
    infile: String,
    outdir: String,
}

fn write_contig_training_data(contig: &Contig, base_outdir: &String){
    let mut outdir = base_outdir.clone();
    let base_outdir_exists = exists(&outdir);
    match base_outdir_exists {
        Err(e) => panic!("Unable to verify if output directory {} exists: {}", outdir, e),
        Ok(false) => {
            std::fs::create_dir(&outdir).expect(format!("Unable to create output directory {}", outdir).as_str());
        },
        Ok(true) => {},    
    }
    for taxon in contig.taxonomy.iter() {
        if !outdir.ends_with("/") {
            outdir.push_str("/");
        }
    
        outdir.push_str(taxon.to_lowercase().as_str());
        let outdir_exists = exists(&outdir);
        match outdir_exists {
            Err(e) => panic!("Unable to verify if output directory {} exists: {}", outdir, e),
            Ok(false) => {
                std::fs::create_dir(&outdir).expect(format!("Unable to create output directory {}", outdir).as_str());
            },
            Ok(true) => {},    
        }
    }
    for gene in contig.genes.iter(){
        let mut prot_number = 1;
        for protein in gene.proteins.iter(){
            let mut protein_outfile = outdir.clone();
            if !protein_outfile.ends_with("/") {
                protein_outfile.push_str("/"); 
            }
            if protein.uniprot_name != "" {
                protein_outfile.push_str(format!("{}_{}_p{}", contig.name,protein.uniprot_name, prot_number).replace(" ", "_").replace(":", "_").to_lowercase().as_str());
            }
            else{
                protein_outfile.push_str(format!("{}_{}_p{}", contig.name, protein.name, prot_number).replace(" ", "_").to_lowercase().as_str());
            }
            // build the set of hmm states we want to train the NN to model and the HMM to output
            // remember here that the dna sequence is always ordered from 5' to 3' regardless of source strand
            let mut hmm_states:Vec<char> = vec![ILLEGAL_STATE; protein.dna_sequence.len()];
            let start_codon_start : usize = match protein.exons[0].complement{
                true => {max(protein.coding_start, protein.exons[0].end) as usize}
                false => {max(protein.coding_start, protein.exons[0].start) as usize}
            };
            let subset: &[char] = &protein.dna_sequence[(max(10, start_codon_start) - 10.. min(start_codon_start +10, protein.dna_sequence.len()-1))];
            let substring: String = subset.iter().collect();
            let start_setby = match start_codon_start == protein.coding_start as usize{
                true => {"coding_start"},
                false => {"exon"},
            };
 
            // First: label the start and stop codons
            assert!(start_codon_start < (protein.coding_end -2) as usize);
            hmm_states[start_codon_start+2]= START;
            if protein.dna_sequence[start_codon_start] == 'G' || 
            protein.dna_sequence[start_codon_start +1] != 'T' ||
            protein.dna_sequence[start_codon_start +2] != 'G'{
                println!("{} {}", protein.name, substring);
                eprintln!("Protein {} had unusual start codon of {}{}{}, complement was {}, first AA was {}", protein.name, protein.dna_sequence[start_codon_start],
                protein.dna_sequence[start_codon_start +1],protein.dna_sequence[start_codon_start+2], protein.exons[0].complement, protein.protein_sequence[0] );
            }
/*        else{

                eprintln!("Protein {} had start codon in expected position, complement was {}, start_setby = {}", protein.name, protein.exons[0].complement, start_setby);
            } */
            //everything before the start codon is intergenic
            for i in 0..start_codon_start +2{
                hmm_states[i as usize] = INTERGENIC;
                }


                
    
            let mut outfile =File::create(protein_outfile).expect("Unable to create output file");
            outfile.write_all(protein.dna_sequence.iter().collect::<String>().as_bytes()).expect("Unable to write protein sequence to file");
            outfile.write_all("\n".as_bytes()).expect("Unable to write newline to file");
            outfile.write_all(hmm_states.iter().collect::<String>().as_bytes()).expect("Unable to write boundaries to file");   
            prot_number += 1;
        }
    }
}
    


fn main() {
    let args = Cli::parse();

    let mut base_path = args.sourcedir.clone();
    if !base_path.ends_with("/") {
        base_path.push_str("/");
    }

    for entry in glob(format!("{}{}*", base_path, args.infile).as_str()).expect("Failed to read glob pattern") {
        match entry {
            Ok(path) => process_ensembl_file(path.to_str().unwrap(), &args),
            Err(e) => println!("Error processing file: {:?}", e),
        }
    }

}

fn process_ensembl_file(filepath: &str, args: &Cli) {
    println!("Processing file {}", filepath);
    let file = File::open(filepath).expect(format!("Unable to open file {}", filepath).as_str());
    let reader = BufReader::new(file);
    let mut contig = Contig{dna_sequence: Vec::new(),
        genes: Vec::new(),
        length: 0,
        organism: "".to_string(),
        taxonomy: Vec::new(),
        name: "".to_string()};
    let mut lines:Peekable<Lines<BufReader<File>>> = reader.lines().peekable();
    while let Some(line_result) = lines.next() { // there's still data in the file
        let line = line_result.unwrap();
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        match fields[0] {
            "ID" => {
                assert!(contig.name == ""); // make sure we've processed any previous contig
                contig.name = fields[1].to_string();
                contig.length = fields[5].parse::<i64>().unwrap();
                println!("Starting to process contig {}", contig.name);
                //           println!("contig length is {}", contig.length);
            },
            "OS" => {
                contig.organism = fields[1..].join(" ").to_string();
            },
            "OC" => {
                let old_taxonomy = contig.taxonomy.clone();
                contig.taxonomy = parse_taxonomy(&line, &mut lines);
                if old_taxonomy.len() > 0 {
                    if old_taxonomy != contig.taxonomy {
                        panic!("Contig {}: contains multiple taxonomys old taxonomy was {:?}, new taxonomy is {:?}", contig.name, old_taxonomy, contig.taxonomy);
                    }
                }           
            },
            "FT" => {
                if fields[1] == "gene" {
                    let temp_gene = parse_gene(&mut lines);
                    if temp_gene.proteins.len() > 0 {
                        contig.genes.push(temp_gene);
                    }
                } else {};
            }
            "SQ" => {
                contig.dna_sequence = parse_sequence(fields, &mut lines);
            },
            "//" => {
                parse_contig(&mut contig);
                write_contig_training_data(&contig, &args.outdir);
                println!("Finished processing contig {}", contig.name);
                contig = Contig{dna_sequence: Vec::new(),
                    genes: Vec::new(),
                    length: 0,
                    organism: "".to_string(),
                    taxonomy: Vec::new(),
                    name: "".to_string()};
            },
            "XX" => {}, // spacer line
            "AC" => {}, // accession line, we don't care about this
            "SV" => {}, // sequence version line, we don't care about this
            "DT" => {}, // date line, we don't care about this
            "DE" => {}, // description line, we don't care about this
            "KW" => {}, // keywords line, we don't care about this
            "CC" => {}, // comment line, we don't care about this
            "FH" => {}, // feature header line, we don't care about this
            _ => println!("Unexpected line: {}", line),
        }
    }

   assert!(contig.name == ""); // make sure we've processed all the contigs

}
