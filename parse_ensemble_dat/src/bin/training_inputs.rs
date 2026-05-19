// parse_ensemble_dat: takes an ensemble dat file that describes a genome and
// parses it into DNA chunks that contain proteins along with vectors that
// describe their introns, exons, start points, and end points.


use std::fs::File;
use std::fs::exists;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::io::{BufRead,BufReader, Lines};

use std::iter::Peekable;
use clap::Parser;

use parse_ensemble_dat::*;

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
            println!("writing training data for protein {} to file {}", protein.name, protein_outfile);
            prot_number += 1;
            }
        }
    }
}

fn main() {
    let args = Cli::parse();

    let file = File::open(args.infile.clone()).expect(format!("Unable to open {}", args.infile).as_str());

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
                if contig.name != "" {
                    parse_contig(&mut contig);
                    write_contig_training_data(&contig, &args.outdir);
                }
                contig = Contig { dna_sequence: Vec::new(), genes: Vec::new(), length: 0, organism: "".to_string(), taxonomy: Vec::new(), name: "".to_string() };
                contig.name = fields[1].to_string();
                contig.length = fields[5].parse::<i64>().unwrap();
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

    // ok, we've processed the data file, now extract what we need from it.

}
