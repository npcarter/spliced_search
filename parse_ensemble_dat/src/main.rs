use std::fs::File;
use std::io::{self, BufRead,BufReader, Lines};
use std::path::Path;
use std::iter::Peekable;
use clap::Parser;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    infile: String,
}

struct exon{
    start: i64,
    end: i64,
}

struct intron{
    start: i64,
    end: i64,
}

struct gene{
    start: i64, // start position of the gene in the contig
    end: i64, // end position of the gene in the contig.
    // gene runs from start to end, inclusive
    dna_sequence: Vec<char>,
    protein_sequence: Vec<char>,
    exons: Vec<exon>,
    introns: Vec<intron>,
}
struct contig{
    length: i64,
    name: String,
    organism: String,
    dna_sequence: Vec<char>,
    genes: Vec<gene>,
}

fn process_contig(contig: &contig){
    println!("Contig {} has length {} and organism {}", contig.name, contig.length, contig.organism);
}

fn parse_sequence(fields: Vec<&str>, lines: &mut Peekable<Lines<BufReader<File>>>) -> Vec<char>{
    let mut sequence:Vec<char> = Vec::new();
    let length:i64 = fields[1].parse::<i64>().unwrap();
    let mut residue_count:i64 = 0;
    while (&mut *lines).peek().unwrap().unwrap() != "//".to_string(){
        let line = (&mut *lines).next().unwrap().unwrap();
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        let end:i64 = fields[-1].parse::<i64>().unwrap();
        for residues in fields[0..-1].iter(){
            for residue in residues.chars() {
                sequence.push(residue);
                residue_count += 1;
            }
        }
        if residue_count != end{println!("Error: residue count {} does not match end position {}", residue_count, end);}
    }

    sequence
}

fn parse_gene(fields: Vec<&str>, lines: &mut Peekable<Lines<BufReader<File>>>) -> gene{
    let mut the_gene = gene{start:0, end:0, dna_sequence: Vec::new(), protein_sequence: Vec::new(), exons: Vec::new(),
    introns: Vec::new()};
    let mut next_line:String = (&mut *lines).peek().unwrap().unwrap();
    let mut next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    while next_fields[0] != "XX" && next_fields[0] != "SQ" && next_fields[0] != "//" &&
        !(next_fields[0] == "FT" && next_fields[1] == "gene"){
        // We're still in the gene


        next_line = (&mut *lines).peek().unwrap().unwrap();
        next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    }
    the_gene
}
fn main() {
    let args = Cli::parse();

    let file = File::open(args.infile.clone()).expect(format!("Unable to open {}", args.infile).as_str());
    let reader = io::BufReader::new(file);
    let mut contig = contig{dna_sequence: Vec::new(),
        genes: Vec::new(),
        length: 0,
        organism: "".to_string(),
        name: "".to_string()};
    let mut lines:Peekable<Lines<BufReader<File>>> = reader.lines().peekable();
    while let Some(Ok(line)) = lines.next() { // there's still data in the file
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        match fields[0]{
            "ID" => {
                contig.name = fields[1].to_string();
                contig.length = fields[5].parse::<i64>().unwrap();
                break;
            },
            "OS" => {contig.organism = fields[1..].join(" ").to_string(); break;},
            "FT" => {
                if fields[1] == "gene"{contig.genes.push(parse_gene(fields[2..], &lines))}
                else{println!("Unexpected FT line: {}", fields.join(" "))}
                break;
            }
            "SQ" =>{contig.dna_sequence = parse_sequence(fields, &mut lines); break;},
            "//" => {process_contig(&contig); break;},
            _ => println!("Unexpected line: {}", line)
        }


    }

}
